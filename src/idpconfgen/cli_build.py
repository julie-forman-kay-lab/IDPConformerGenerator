"""
Builds IDP conformers.

Build from a database of torsion angles and secondary structure
information. Database is as created by `idpconfgen torsions` CLI.

USAGE:
    $ idpconfgen build -db torsions.json -seq MMMMMMM...

"""
import argparse
import re
import sys
from collections import Counter, namedtuple
from functools import partial
from itertools import cycle
from multiprocessing import Pool, Queue
# from numbers import Number
from random import choice as randchoice
from random import randint
from time import time

import numpy as np
from numba import njit

import idpcpp
#from idpconfgen.cpp.faspr import faspr_sidechains as fsc
from idpconfgen import log, Path
from idpconfgen import assert_type as AT
from idpconfgen.core.build_definitions import (
    amber_pdbs,
    atom_labels_amber,
    backbone_atoms,
    bonds_equal_3_inter,
    bonds_le_2_inter,
    build_bend_H_N_C,
    build_bend_angles_CA_C_Np1,
    build_bend_angles_Cm1_N_CA,
    build_bend_angles_N_CA_C,
    distance_H_N,
    distances_CA_C,
    distances_C_Np1,
    distances_N_CA,
    expand_topology_bonds_apart,
    generate_residue_template_topology,
    inter_residue_connectivities,
    n_terminal_h_coords_at_origin,
    read_ff14SB_params,
    sidechain_templates,
    topology_3_bonds_apart,
    )
from idpconfgen.core.definitions import aa1to3  # , vdW_radii_dict
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs import libcli
from idpconfgen.libs.libcalc import (
    # calc_all_vs_all_dists_square,
    calc_all_vs_all_dists,
    calc_residue_num_from_index,
    calc_torsion_angles,
    make_coord_Q,
    make_coord_Q_CO,
    make_coord_Q_COO,
    place_sidechain_template,
    rotate_coordinates_Q_njit,
    )
from idpconfgen.libs.libfilter import aligndb, regex_search
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libpdb import atom_line_formatter
from idpconfgen.libs.libtimer import timeme, ProgressCounter


faspr_sc = idpcpp.faspr_sidechains
dun2010bbdep_path = Path(
    Path(__file__).myparents(),
    'core',
    'data',
    'dun2010bbdep.bin',
    ).str()
_name = 'build'
_help = 'Builds conformers from database.'


_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    '-db',
    '--database',
    help='The IDPConfGen database.',
    required=True,
    )

ap.add_argument(
    '-seq',
    '--input_seq',
    help='The Conformer residue sequence.',
    required=True
    )

ap.add_argument(
    '-nc',
    '--nconfs',
    help='Number of conformers to build.',
    default=1,
    type=int,
    )

ap.add_argument(
    '-dr',
    '--dssp-regexes',
    help='Regexes used to search in DSSP',
    default='(?=(L{2,6}))',
    nargs='+',
    )

ap.add_argument(
    '-dsd',
    '--disable-sidechains',
    help='Whether or not to compute sidechais. Defaults to True.',
    action='store_true',
    )

# TODO: these three parameters must be discontinued
# libcli.add_argument_vdWb(ap)
# libcli.add_argument_vdWr(ap)
# libcli.add_argument_vdWt(ap)
libcli.add_argument_ncores(ap)


SLICES = []
ANGLES = None
CONF_NUMBER = Queue()

ConfMasks = namedtuple(
    'ConfMaks',
    [
        'bb3',
        'bb4',
        'NHs',
        'COs',
        'Hterm',
        'OXT1',
        'OXT2',
        'non_Hs',
        'non_Hs_non_OXT',
        'H1_N_CA_HA',
        ]
    )


# Other functions should have the same APIs:
# parameters = input_seq
def init_faspr_sidechains(input_seq):
    """."""
    dun2010db = dun2010bbdep_path
    faspr_func = faspr_sc
    def compute_faspr_sidechains(coords):
        """Does calculation."""
        return faspr_func(coords, input_seq, dun2010db)

    return compute_faspr_sidechains


compute_sidechains = {
    'faspr': init_faspr_sidechains,
    }


def init_confmasks(atom_labels):
    """
    Create a ConfMask object (namedtuple).

    ConfMask is a named tuple which attributes are integer masks for the
    respective groups.

    Parameters
    ----------
    atom_labels : array-like
        The atom names of the protein.

    Returns
    -------
    namedtuple

    Notes
    -----
    ConfMask attributes map to the following atom groups:

    bb3 : N, CA, C
    bb4 : N, CA, C, O
    NHs : amide protons
    Hterm : N-terminal protons
    OXT1 : O atom of C-terminal carboxyl group
    OXT2 : OXT atom of the C-terminal carboxyl group
    non_Hs : all but hydrogens
    non_Hs_non_OXT : all but hydrogens and the only OXT atom
    H1_N_CA_HA : these four atoms from the first residue
    """
    bb3 = np.where(np.isin(atom_labels, ('N', 'CA', 'C')))[0]
    assert len(bb3) % 3 == 0

    bb4 = np.where(np.isin(atom_labels, ('N', 'CA', 'C', 'O')))[0]
    assert len(bb4) % 4 == 0

    NHs = np.where(atom_labels == 'H')[0]

    # last O belongs to the C-term carboxyl, we don't want it in the carbonyl
    # mask
    COs = np.where(atom_labels == 'O')[0][:-1]

    OXT1 = np.where(atom_labels == 'O')[0][-1]
    # the actual OXT atom, this is only one O of the C-term carboxyl pair
    OXT2 = np.where(atom_labels == 'OXT')[0][0]  # int instead of list

    rr = re.compile(r'H+')
    hs_match = np.vectorize(lambda x: bool(rr.match(x)))
    non_Hs = np.where(np.logical_not(hs_match(atom_labels)))[0]
    non_Hs_non_OXT = non_Hs[:-1]

    # used to rotate the N-terminal Hs to -60 degrees to  HA during
    # the building process
    _H1_idx = np.where(atom_labels == 'H1')[0]
    assert len(_H1_idx) == 1

    # of the first residue
    _N_CA_idx = np.where(np.isin(atom_labels, ('N', 'CA')))[0][:2]
    assert len(_N_CA_idx) == 2, _N_CA_idx

    _HA_HA2_idx = np.where(np.isin(atom_labels, ('HA', 'HA2')))[0][0:1]
    assert len(_HA_HA2_idx) == 1

    H1_N_CA_HA = list(_H1_idx) + list(_N_CA_idx) + list(_HA_HA2_idx)
    assert len(H1_N_CA_HA) == 4

    Hterm = np.where(np.isin(atom_labels, ('H1', 'H2', 'H3')))[0]
    assert len(Hterm) == 3

    conf_mask = ConfMasks(
        bb3=bb3,
        bb4=bb4,
        NHs=NHs,
        COs=COs,
        Hterm=Hterm,
        OXT1=OXT1,
        OXT2=OXT2,
        non_Hs=non_Hs,
        non_Hs_non_OXT=non_Hs_non_OXT,
        H1_N_CA_HA=H1_N_CA_HA,
        )

    return conf_mask


def get_cycle_distances_backbone():
    """
    Return an inifinite iterator of backbone atom distances.
    """
    return cycle((
            distances_N_CA,  # used for OMEGA
            distances_CA_C,  # used for PHI
            distances_C_Np1,  # used for PSI
            ))


def get_cycle_bend_angles():
    """
    Return an infinite iterator of the bend angles.
    """
    return cycle((
        build_bend_angles_Cm1_N_CA,  # used for OMEGA
        build_bend_angles_N_CA_C,  # used for PHI
        build_bend_angles_CA_C_Np1,  # used for PSI
        ))


def main(
        input_seq,
        database,
        dssp_regexes=r'(?=(L{2,6}))',
        func=None,
        nconfs=1,
        ncores=1,
        **kwargs,
        ):
    """
    Execute main client logic.

    Distributes over processors.
    """
    global ANGLES
    # Calculates how many conformers are built per core
    if nconfs < ncores:
        ncores = 1
        core_chunks = nconfs
        remaining_chunks = 0
    else:
        core_chunks = nconfs // ncores
        # in case nconfs is not multiple of ncores, builds the remaining confs
        # at the end
        remaining_chunks = nconfs % ncores

    _slices, ANGLES = read_db_to_slices(database, dssp_regexes, ncores=ncores)
    SLICES.extend(_slices)
    # creates a list of reversed numbers to name the conformers
    for i in range(1, nconfs + 1):
        CONF_NUMBER.put(i)

    # prepars execution function
    execute = partial(
        _build_conformers,
        input_seq=input_seq,  # string
        nconfs=core_chunks,  # int
        **kwargs,
        )

    start = time()
    with Pool(ncores) as pool:
        imap = pool.imap(execute, range(ncores))
        for _ in imap:
            pass

    if remaining_chunks:
        execute(core_chunks * ncores, nconfs=remaining_chunks)

    log.info(f'{nconfs} conformers built in {time() - start:.3f} seconds')


def read_db_to_slices(database, dssp_regexes, ncores=1):
    """Create database base of slice and angles."""
    # reads db dictionary from disk
    db = read_dictionary_from_disk(database)
    log.info(f'Read DB with {len(db)} entries')

    # reads and prepares IDPConfGen data base
    timed = partial(timeme, aligndb)
    pdbs, angles, dssp, resseq = timed(db)

    # searchs for slices in secondary structure, according to user requests
    timed = partial(timeme, regex_search, ncores=ncores)
    dssp_regexes = \
        [dssp_regexes] if isinstance(dssp_regexes, str) else dssp_regexes
    slices = []
    for dssp_regex_string in dssp_regexes:
        slices.extend(timed(dssp, dssp_regex_string))
    log.info(f'Found {len(slices)} indexes for {dssp_regexes}')

    return slices, angles


def create_conformer_labels(input_seq, input_seq_3_letters):
    """
    Create all atom labels model based on an input sequence.

    The use of `input_seq_3_letters` is still experimental. Is likely
    to be refactored.
    """
    # /
    # prepares data based on the input sequence
    # considers sidechain all-atoms
    atom_labels = np.array(make_list_atom_labels(input_seq, atom_labels_amber))
    num_atoms = len(atom_labels)

    # /
    # per atom labels
    residue_numbers = np.empty(num_atoms, dtype=np.int)
    residue_labels = np.empty(num_atoms, dtype='<U3')

    # generators
    _res_nums_gen = gen_residue_number_per_atom(atom_labels, start=1)
    _res_labels_gen = \
        gen_3l_residue_labels_per_atom(input_seq_3_letters, atom_labels)

    # fills empty arrays from generators
    _zipit = zip(range(num_atoms), _res_nums_gen, _res_labels_gen)
    for _i, _num, _label in _zipit:
        residue_numbers[_i] = _num
        residue_labels[_i] = _label

    # maniatic cleaning from pre-function isolation
    del _res_labels_gen, _res_nums_gen, _zipit

    # ensure
    assert len(residue_numbers) == num_atoms
    assert len(residue_labels) == num_atoms, (len(residue_labels), num_atoms)
    # ?
    return atom_labels, residue_numbers, residue_labels


# private function because it depends on the global `CONF_NUMBER`
# which is assembled in `main()`
def _build_conformers(
        *args,
        input_seq=None,
        conformer_name='conformer',
        nconfs=1,
        **kwargs):
    """
    Arrange building of conformers and saves them to PDB files.
    """
    ROUND = np.round

    # TODO: this has to be parametrized for the different HIS types
    input_seq_3_letters = [
        'HIP' if _res == 'H' else aa1to3[_res]
        for _res in input_seq
        ]

    builder = conformer_generator(input_seq=input_seq, **kwargs)

    atom_labels, residue_numbers, residue_labels = next(builder)

    for conf_n in range(nconfs):

        coords = next(builder)

        pdb_string = gen_PDB_from_conformer(
            input_seq_3_letters,
            atom_labels,
            residue_numbers,
            ROUND(coords, decimals=3),
            )

        fname = f'{conformer_name}_{CONF_NUMBER.get()}.pdb'

        with open(fname, 'w') as fout:
            fout.write(pdb_string)

    del builder
    return



# the name of this function is likely to change in the future
def conformer_generator(
        *,
        input_seq=None,
        generative_function=None,
        disable_sidechains=True,
        sidechain_method='faspr',
        # TODO: these parameters must be discontinued
        # vdW_bonds_apart=3,
        # vdW_tolerance=0.4,
        # vdW_radii='tsai1999',
        ):
    """
    Build conformers.

    `conformer_generator` is actually a Python generator. Examples on
    how it works:

    Note that all arguments are **named** arguments.

    >>> builder = conformer_generator(
    >>>    input_seq='MGAETTWSCAAA'  # the primary sequence of the protein
    >>>    )

    `conformer_generator` is a generator, you can instantiate it simply
    providing the residue sequence of your protein of interest.

    The **very first** iteration will return the labels of the protein
    being built. Labels are sorted by all atom models. Likewise,
    `residue_number` and `residue_labels` sample **all atoms**. These
    three are numpy arrays and can be used to index the actual coordinates.

    >>> atom_labels, residue_numbers, residue_labels = next(builder)

    After this point, each iteraction `next(builder)` yields the coordinates
    for a new conformer. There is no limit in the generator.

    >>> new_coords = next(builder)

    `new_coords` is a (N, 3) np.float64 array where N is the number of
    atoms. As expected, atom coordinates are aligned with the labels
    previously generated.

    When no longer needed,

    >>> del builder

    Should delete the builder generator.

    You can gather the coordinates of several conformers in a single
    multi dimensional array with the following:

    >>> builder = conformer_generator(
    >>>     input_seq='MGGGGG...',
    >>>     generative_function=your_function)
    >>>
    >>> atoms, res3letter, resnums = next(builder)
    >>>
    >>> num_of_conformers = 10_000
    >>> shape = (num_of_conformers, len(atoms), 3)
    >>> all_coords = np.empty(shape, dtype=float64)
    >>>
    >>> for i in range(num_of_conformers):
    >>>     all_coords[i, :, :] = next(builder)
    >>>

    Parameters
    ----------
    input_seq : str, mandatory
        The primary sequence of the protein being built in FASTA format.
        `input_seq` will be used to generate the whole conformers' and
        labels arrangement.
        Example: "MAGERDDAPL".

    generative_function : callable, optional
        The generative function used by the builder to retrieve torsion
        angles during the building process.

        The builder expects this function to receive two parameters:
            - `nres`, the residue chunk size to get angles from
            - `cres`, the next residue being built. For example,
                with cres=10, the builder will expect a minimum of three
                torsion angles (phi, psi, omega) for residue 10.

        Depending on the nature of the `generative function` the two
        pameters may be ignored by the function itself (use **kwargs
        for that purpose).

        If `None` provided, the builder will use the internal `SLIDES`
        and `ANGLES` variables and will assume the `cli_build.main` was
        executed priorly, or that ANGLES and SLICES were populated
        properly.

    disable_sidechains : bool
        Disables sidechain creation. Defaults to `False`, computes
        sidechains.

    nconfs : int
        The number of conformers to build.
    """
    assert input_seq, f'`input_seq` must be given! {input_seq}'
    BUILD_BEND_H_N_C = build_bend_H_N_C
    # CALC_DISTS = calc_all_vs_all_dists_square
    # TODO
    CALC_DISTS = calc_all_vs_all_dists
    CALC_TORSION_ANGLES = calc_torsion_angles
    DISTANCE_NH = distance_H_N
    MAKE_COORD_Q_COO_LOCAL = make_coord_Q_COO
    MAKE_COORD_Q_CO_LOCAL = make_coord_Q_CO
    MAKE_COORD_Q_LOCAL = make_coord_Q
    NAN = np.nan
    NANSUM = np.nansum
    N_TERMINAL_H = n_terminal_h_coords_at_origin
    PI2 = np.pi * 2
    PLACE_SIDECHAIN_TEMPLATE = place_sidechain_template
    RC = randchoice
    RAD_60 = np.radians(60)
    RINT = randint
    ROT_COORDINATES = rotate_coordinates_Q_njit
    ROUND = np.round
    SIDECHAIN_TEMPLATES = sidechain_templates
    VEC_1_X_AXIS = np.array([1, 0, 0])
    angles = ANGLES
    slices = SLICES
    # bellow an optimization for the builder loop

    # TODO: correct for HIS/HIE/HID/HIP
    def translate_seq_to_3l(input_seq):
        return [
            'HIP' if _res == 'H' else aa1to3[_res]
            for _res in input_seq
            ]

    r_input_seq = input_seq
    r_input_seq_3_letters = translate_seq_to_3l(input_seq)
    del input_seq

    # semantic exchange for speed al readibility
    with_sidechains = not(disable_sidechains)

    if with_sidechains:
        calc_sidechains = compute_sidechains[sidechain_method](r_input_seq)


    # tests generative function complies with implementation requirements
    if generative_function:
        try:
            generative_function(nres=1, cres=0)
        except Exception as err:  # this is generic Exception on purpose
            errmsg = (
                'The `generative_function` provided is not compatible with '
                'the building process. Please read `build_conformers` docstring'
                ' for more details.'
                )
            raise IDPConfGenException(errmsg) from err

    # Start building process

    r_atom_labels, r_residue_numbers, r_residue_labels = \
        create_conformer_labels(r_input_seq, r_input_seq_3_letters)

    # first yields the labels for the caller to handle them
    # because later the conformers will be yielded and all share the
    # same labels
    yield r_atom_labels, r_residue_numbers, r_residue_labels

    r_num_atoms = len(r_atom_labels)
    #num_ij_pairs = num_atoms * (num_atoms - 1) // 2

    r_acoeff, r_bcoeff, r_charges_ij, r_bonds_ge_3_mask = \
        create_energy_func_params(r_atom_labels, r_residue_numbers, r_residue_labels)

    # /
    # Prepares alanine/proline template backbone
    ala_pro_seq = ''.join('A' if res not in ('P', 'G') else res for res in r_input_seq)
    ala_pro_seq_3l = translate_seq_to_3l(ala_pro_seq)
    atom_labels, residue_numbers, residue_labels = \
        create_conformer_labels(ala_pro_seq, ala_pro_seq_3l)

    num_atoms = len(atom_labels)
    num_ij_pairs = num_atoms * (num_atoms - 1) // 2

    ap_acoeff, ap_bcoeff, ap_charges_ij, ap_bonds_ge_3_mask = \
        create_energy_func_params(atom_labels, residue_numbers, residue_labels)

    # TODO: needed for other energy functions
    #atoms_VDW = 0.5 * sigmas_ii * 2**(1/6)

    # /
    # creates coordinate data-structures and,
    # creates index-translated boolean masks

    r_confmasks = init_confmasks(r_atom_labels)
    ap_confmasks = init_confmasks(atom_labels)

    # create coordinates and views
    coords = np.full((num_atoms, 3), NAN, dtype=np.float64)

    # +2 because of the dummy coordinates required to start building.
    # see later adding dummy coordinates to the structure seed
    bb = np.full((ap_confmasks.bb3.size + 2, 3), NAN, dtype=np.float64)
    bb_real = bb[2:, :]  # backbone coordinates without the dummies
    # coordinates for the carbonyl oxigen atoms
    bb_CO = np.full((ap_confmasks.COs.size, 3), NAN, dtype=np.float64)

    # notice that NHydrogen_mask does not see Prolines
    bb_NH = np.full((ap_confmasks.NHs.size, 3), NAN, dtype=np.float64)

    # the following array serves to avoid placing HN in Proline residues
    residue_labels_bb_simulating = residue_labels[ap_confmasks.bb3]

    ss_masks = create_sidechains_masks_per_residue(
        residue_numbers,
        atom_labels,
        backbone_atoms,
        )
    # ?

    all_atoms_coords = np.full((r_num_atoms, 3), NAN, dtype=np.float64)
    all_atoms_masks = init_confmasks(r_atom_labels)


    # /
    # creates seed coordinates:
    # 1st) a dummy atom at the y-axis to build the first atom
    # 2nd) N-terminal N-atom is at 0, 0, 0
    # 3rd) CA atom of the firs residue is at the x-axis
    # # coordinates are created always from the parameters in the core
    # # definitions of IDPConfGen
    dummy_CA_m1_coord = np.array((0.0, 1.0, 1.0))
    dummy_C_m1_coord = np.array((0.0, 1.0, 0.0))
    n_terminal_N_coord = np.array((0.0, 0.0, 0.0))
    #n_terminal_CA_coord = np.array((distances_N_CA[ala_pro_seq[0]], 0.0, 0.0))

    # seed coordinates array
    seed_coords = np.array((
        dummy_CA_m1_coord,
        dummy_C_m1_coord,
        n_terminal_N_coord,
        #n_terminal_CA_coord,
        ))
    # ?

    # /
    # prepares method binding
    bbi0_register = []
    bbi0_R_APPEND = bbi0_register.append
    bbi0_R_POP = bbi0_register.pop
    bbi0_R_CLEAR = bbi0_register.clear

    COi0_register = []
    COi0_R_APPEND = COi0_register.append
    COi0_R_POP = COi0_register.pop
    COi0_R_CLEAR = COi0_register.clear

    NHi0_register = []
    NHi0_R_APPEND = NHi0_register.append
    NHi0_R_POP = NHi0_register.pop
    NHi0_R_CLEAR = NHi0_register.clear

    res_R = []  # residue number register
    res_R_APPEND = res_R.append
    res_R_POP = res_R.pop
    res_R_CLEAR = res_R.clear
    # ?

    # /
    # required inits
    broke_on_start_attempt = False
    start_attempts = 0
    max_start_attempts = 500  # maximum attempts to start a conformer
    # because we are building from a experimental database there can be
    # some angle combinations that fail on our validation process from start
    # if this happens more than `max_start_attemps` the production is canceled.
    # ?

    # /
    # STARTS BUILDING
    conf_n = 1
    while 1:
        # prepares cycles for building process
        # cycles need to be regenerated every conformer because the first
        # atom build is a C and the last atom built is the CA, which breaks
        # periodicity
        # all these are dictionaries
        bond_lens = get_cycle_distances_backbone()
        bond_bend = get_cycle_bend_angles()

        # in the first run of the loop this is unnecessary, but is better to
        # just do it once than flag it the whole time
        coords[:, :] = NAN
        bb[:, :] = NAN
        bb_CO[:, :] = NAN
        bb_NH[:, :] = NAN
        for _mask, _coords in ss_masks:
            _coords[:, :] = NAN

        bb[:3, :] = seed_coords  # this contains a dummy coord at position 0

        # add N-terminal hydrogens to the origin
        coords[ap_confmasks.Hterm, :] = N_TERMINAL_H
        current_Hterm_coords = N_TERMINAL_H

        bbi = 1  # starts at 1 because there are two dummy atoms
        bbi0_R_CLEAR()
        bbi0_R_APPEND(bbi)

        COi = 0  # carbonyl atoms
        COi0_R_CLEAR()
        COi0_R_APPEND(COi)

        # NHi is 0 because
        # the first residue, remember the first residue as 'H1,2,3'
        # see bb_NH_builder definition
        NHi = 0
        NHi0_R_CLEAR()
        NHi0_R_APPEND(NHi)

        # residue integer number
        current_res_number = 0
        res_R_CLEAR()
        res_R_APPEND(current_res_number)

        backbone_done = False
        number_of_trials = 0
        # TODO: use or not to use number_of_trials2? To evaluate in future.
        number_of_trials2 = 0
        number_of_trials3 = 0
        # run this loop until a specific BREAK is triggered
        while 1:  # 1 is faster than True :-)

            # I decided to use an if-statement here instead of polymorph
            # the else clause to a `generative_function` variable because
            # the resulting overhead from the extra function call and
            # **kwargs handling was greater then the if-statement processing
            # https://pythonicthoughtssnippets.github.io/2020/10/21/PTS14-quick-in-if-vs-polymorphism.html
            if generative_function:
                agls = generative_function(
                    nres=RINT(1, 6),
                    #cres=abs(bbi - 2) // 3,  #TODO: remove
                    cres=calc_residue_num_from_index(bbi)
                    )

            else:
                # following `aligndb` function,
                # `angls` will always be cyclic with:
                # phi - psi - omega - phi - psi - omega - (...)
                agls = angles[RC(slices), :].ravel()

            # index at the start of the current cycle
            try:
                for torsion in agls:
                    # TODO: remove this comment
                    # #bbi -2 makes the match, because despite the first atom
                    # #being placed is a C, it is placed using the PHI angle of a
                    # #residue. And PHI is the first angle of that residue angle
                    # #set.
                    #current_res_number = abs(bbi - 2) // 3  #TODO: remove
                    current_res_number = calc_residue_num_from_index(bbi)
                    current_residue = ala_pro_seq[current_res_number]

                    # bbi is the same for bb_real and bb, but bb_real indexing
                    # is displaced by 1 unit from bb. So 2 in bb_read is the
                    # next atom of 2 in bb.
                    bb_real[bbi, :] = MAKE_COORD_Q_LOCAL(
                        bb[bbi - 1, :],
                        bb[bbi, :],
                        bb[bbi + 1, :],
                        next(bond_lens)[current_residue],
                        next(bond_bend)[current_residue],
                        torsion,
                        )
                    bbi += 1

            except IndexError:
                # IndexError happens when the backbone is complete
                # in this protocol the last atom build was a carbonyl C
                # bbi is the last index of bb + 1, and the last index of
                # bb_real + 2

                # activate flag to finish loop at the end
                backbone_done = True
                print('adding COO')

                # add the carboxyls
                coords[[ap_confmasks.OXT2, ap_confmasks.OXT1]] = \
                    MAKE_COORD_Q_COO_LOCAL(bb[-2, :], bb[-1, :])
                print(coords[[ap_confmasks.OXT2, ap_confmasks.OXT1]])

            
            finally:
                print('bbi register: ', bbi0_register[-1], bbi)

            # builds carbonyl atoms. Two situations can happen here:
            # 1) the backbone is not complete - the last atom is CA
            # 2) the backbone is complete - the last atom is C
            # whichever the case, with `bbi - 1` applying on the range bellow
            # will build carbonyls for the new chain expect for the last C
            # if the chain is completed.
            # this is so because range(X, Y, Z) equals range(X, Y-1, Z)
            # if Y-1 is not in range(X, Y, Z). And, this is the case for N-CA
            # pair.
            for k in range(bbi0_register[-1] + 1, bbi - 1, 3):
                bb_CO[COi, :] = MAKE_COORD_Q_CO_LOCAL(
                    bb_real[k - 1, :],
                    bb_real[k, :],
                    bb_real[k + 1, :],
                    )
                COi += 1

            # Adds N-H Hydrogens
            if len(bbi0_register) > 1:
                for k in range(bbi0_register[-2] + 2, bbi - 1, 3):

                    if residue_labels_bb_simulating[k] == 'PRO':
                        continue

                    # MAKE_COORD_Q_CO_LOCAL can be used for NH by giving
                    # disntace and bend parameters
                    print(bb_real[k - 1, :])
                    print(bb_real[k, :])
                    print(bb_real[k + 1, :])
                    bb_NH[NHi, :] = MAKE_COORD_Q_CO_LOCAL(
                        bb_real[k - 1, :],
                        bb_real[k, :],
                        bb_real[k + 1, :],
                        distance=DISTANCE_NH,
                        bend=BUILD_BEND_H_N_C,
                        )
                    assert not np.all(np.isnan(bb_NH[NHi, :]))
                    NHi += 1

            # Adds sidechain template structures
            # TODO: remove this if-statement
            if True:#with_sidechains:
                for res_i in range(res_R[-1], current_res_number):  # noqa: E501

                    print('res i for sidechain', res_i)

                    _sstemplate, _sidechain_idxs = \
                        SIDECHAIN_TEMPLATES[ala_pro_seq_3l[res_i]]

                    print('sidechain atoms', atom_labels[ap_confmasks.bb3][res_i * 3:res_i * 3 + 3])

                    sscoords = PLACE_SIDECHAIN_TEMPLATE(
                        bb_real[res_i * 3:res_i * 3 + 3, :],  # from N to C
                        _sstemplate,
                        )

                    ss_masks[res_i][1][:, :] = sscoords[_sidechain_idxs]

                # Transfers coords to the main coord array
                for _smask, _sidecoords in ss_masks[: current_res_number + backbone_done]:  # noqa: E501
                    coords[_smask] = _sidecoords

            # / Place coordinates for energy calculation
            #
            # use `bb_real` to do not consider the initial dummy atom
            coords[ap_confmasks.bb3] = bb_real
            coords[ap_confmasks.COs] = bb_CO
            coords[ap_confmasks.NHs] = bb_NH

            if False:#len(bbi0_register) == 1 and ala_pro_seq[0] != 'G':
                # rotates the N-terminal Hs only if it is the first
                # chunk being built

                # measure torsion angle reference H1 - HA
                _h1_ha_angle = \
                    CALC_TORSION_ANGLES(coords[ap_confmasks.H1_N_CA_HA, :])[0]

                # given any angle calculated along an axis, calculate how
                # much to rotate along that axis to place the
                # angle at 60 degrees
                _rot_angle = _h1_ha_angle % PI2 - RAD_60

                current_Hterm_coords = ROT_COORDINATES(
                    coords[ap_confmasks.Hterm, :],
                    VEC_1_X_AXIS,
                    _rot_angle,
                    )
                coords[ap_confmasks.Hterm, :] = current_Hterm_coords
            # ?

            # /
            # calc energy
            total_energy = calc_energy(
                coords,
                ap_acoeff,
                ap_bcoeff,
                ap_charges_ij,
                ap_bonds_ge_3_mask,
                )

            if False:#total_energy > 0:
                # reset coordinates to the original value
                # before the last chunk added

                # reset the same chunk maximum 5 times,
                # after that reset also the chunk before
                if number_of_trials > 50:
                    bbi0_R_POP()
                    COi0_R_POP()
                    NHi0_R_POP()
                    res_R_POP()
                    number_of_trials = 0
                    number_of_trials2 += 1

                if number_of_trials2 > 5:
                    bbi0_R_POP()
                    COi0_R_POP()
                    NHi0_R_POP()
                    res_R_POP()
                    number_of_trials2 = 0
                    number_of_trials3 += 1

                if number_of_trials3 > 5:
                    bbi0_R_POP()
                    COi0_R_POP()
                    NHi0_R_POP()
                    res_R_POP()
                    number_of_trials3 = 0

                try:
                    _bbi0 = bbi0_register[-1]
                    _COi0 = COi0_register[-1]
                    _NHi0 = NHi0_register[-1]
                    _resi0 = res_R[-1]
                except IndexError:
                    # if this point is reached,
                    # we erased until the beginning of the conformer
                    # discard conformer, something went really wrong
                    broke_on_start_attempt = True
                    break # conformer while loop, starts conformer from scratch

                # clean previously built protein chunk
                bb_real[_bbi0:bbi, :] = NAN
                bb_CO[_COi0:COi, :] = NAN
                bb_NH[_NHi0:NHi, :] = NAN

                # do the same for sidechains
                # ... it is not necessary because in the loop above
                # sidechain coordinates are replaces when repositioned

                # reset also indexes
                bbi = _bbi0
                COi = _COi0
                NHi = _NHi0
                current_res_number = _resi0

                # coords needs to be reset because size of protein next
                # chunks may not be equal
                coords[:, :] = NAN
                coords[ap_confmasks.Hterm, :] = current_Hterm_coords

                # prepares cycles for building process
                # this is required because the last chunk created may have been
                # the final part of the conformer
                if backbone_done:
                    bond_lens = get_cycle_distances_backbone()
                    bond_bend = get_cycle_bend_angles()

                # we do not know if the next chunk will finish the protein
                # or not
                backbone_done = False
                number_of_trials += 1
                continue  # send back to the CHUNK while loop

            # if the conformer is valid
            number_of_trials = 0
            bbi0_R_APPEND(bbi)
            COi0_R_APPEND(COi)
            NHi0_R_APPEND(NHi)
            # the residue where the build process stopped
            res_R_APPEND(current_res_number)

            if backbone_done:
                # this point guarantees all protein atoms are built
                break  # CHUNK while loop
        # END of CHUNK while loop, go up and build the next CHUNK

        if broke_on_start_attempt:
            start_attempts += 1
            if start_attempts > max_start_attempts:
                log.error(
                    'Reached maximum amount of re-starts. Canceling... '
                    f'Built a total of {conf_n} conformers.'
                    )
                return
            broke_on_start_attempt = False
            continue  # send back to the CHUNK while loop

        # we do not want sidechains at this point
        all_atoms_coords[all_atoms_masks.bb4] = coords[ap_confmasks.bb4]
        all_atoms_coords[all_atoms_masks.NHs] = coords[ap_confmasks.NHs]
        #TODO Hterm
        #all_atoms_coords[all_atoms_masks.Hterm] = coords[ap_confmasks.Hterm]
        all_atoms_coords[[all_atoms_masks.OXT2, all_atoms_masks.OXT1], :] = \
            coords[[ap_confmasks.OXT2, ap_confmasks.OXT1], :]


        print(all_atoms_coords[all_atoms_masks.NHs])


        if with_sidechains:

            all_atoms_coords[all_atoms_masks.non_Hs_non_OXT] = \
                calc_sidechains(
                    coords[ap_confmasks.bb4],
                    )

            total_energy = calc_energy(
                all_atoms_coords,
                r_acoeff,
                r_bcoeff,
                r_charges_ij,
                r_bonds_ge_3_mask,
                )

            if total_energy > 0:
                print('Conformer with WORST energy', total_energy)
                continue
            else:
                print(conf_n, total_energy)

        yield all_atoms_coords
        conf_n += 1


def calc_outer_multiplication_upper_diagonal_raw(data, result):
    """
    Calculate the upper diagonal multiplication with for loops.

    The use of for-loop based calculation avoids the creation of very
    large arrays using numpy outer derivates. This function is thought
    to be jut compiled.

    Does not create new data structure. It requires the output structure
    to be provided. Hence, modifies in place. This was decided so
    because this function is thought to be jit compiled and errors with
    the creation of very large arrays were rising. By passing the output
    array as a function argument, errors related to memory freeing are
    avoided.

    Parameters
    ----------
    data : an interable of Numbers, of length N

    result : a mutable sequence, either list of np.ndarray,
             of length N*(N-1)//2
    """
    c = 0
    len_ = len(data)
    for i in range(len_ - 1):
        for j in range(i + 1, len_):
            result[c] = data[i] * data[j]
            c += 1

    # assert result.size == (data.size * data.size - data.size) // 2
    # assert abs(result[0] - data[0] * data[1]) < 0.0000001
    # assert abs(result[-1] - data[-2] * data[-1]) < 0.0000001
    return


njit_calc_multiplication_upper_diagonal_raw = \
    njit(calc_outer_multiplication_upper_diagonal_raw)


def calc_sum_upper_diagonal_raw(data, result):
    """
    Calculate outer sum for upper diagonal with for loops.

    The use of for-loop based calculation avoids the creation of very
    large arrays using numpy outer derivates. This function is thought
    to be jut compiled.

    Does not create new data structure. It requires the output structure
    to be provided. Hence, modifies in place. This was decided so
    because this function is thought to be jit compiled and errors with
    the creation of very large arrays were rising. By passing the output
    array as a function argument, errors related to memory freeing are
    avoided.

    Parameters
    ----------
    data : an interable of Numbers, of length N

    result : a mutable sequence, either list of np.ndarray,
             of length N*(N-1)//2
    """
    c = 0
    len_ = len(data)
    for i in range(len_ - 1):
        for j in range(i + 1, len_):
            result[c] = data[i] + data[j]
            c += 1

    #assert result.size == (data.size * data.size - data.size) // 2
    #assert abs(result[0] - (data[0] + data[1])) < 0.0000001
    #assert abs(result[-1] - (data[-2] + data[-1])) < 0.0000001
    return


njit_calc_sum_upper_diagonal_raw = njit(calc_sum_upper_diagonal_raw)


def calc_LJ_energy(acoeff, bcoeff, dists_ij, to_eval_mask, NANSUM=np.nansum):
    """Calculates Lennard-Jones Energy."""
    # assert dists_ij.size == acoeff.size == bcoeff.size, (dists_ij.size, acoeff.size)
    # assert dists_ij.size == to_eval_mask.size

    ar = acoeff / (dists_ij ** 12)
    br = bcoeff / (dists_ij ** 6)
    energy_ij = ar - br

    # nansum is used because some values in dists_ij are expected to be nan
    # return energy_ij, NANSUM(energy_ij[to_eval_mask])
    return NANSUM(energy_ij[to_eval_mask])


njit_calc_LJ_energy = njit(calc_LJ_energy)


def create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        bonds_intra,
        bonds_inter,
        base_bool=False,
        ):
    """
    Create bool mask array identifying the pairs X bonds apart in ij pairs.

    Given `bonds_intra` and `bonds_inter` criteria, idenfities those ij
    atom pairs in N*(N-1)/2 condition (upper all vs all diagonal) that
    agree with the described bonds.

    Inter residue bonds are only considered for consecutive residues.

    Paramters
    ---------
    atom_labels : iterable, list or np.ndarray
        The protein atom labels. Ex: ['N', 'CA, 'C', 'O', 'CB', ...]

    residue_numbers : iterable, list or np.ndarray
        The protein residue numbers per atom in `atom_labels`.
        Ex: [1, 1, 1, 1, 1, 2, 2, 2, 2, ...]

    residue_labels : iterable, list or np.ndarray
        The protein residue labels per atom in `atom_labels`.
        Ex: ['Met', 'Met', 'Met', ...]

    Depends
    -------
    `gen_ij_pairs_upper_diagonal`
    `gen_atom_pair_connectivity_masks`
    """
    atom_labels_ij_gen = gen_ij_pairs_upper_diagonal(atom_labels)
    residue_numbers_ij_gen = gen_ij_pairs_upper_diagonal(residue_numbers)
    residue_labels_ij_gen = gen_ij_pairs_upper_diagonal(residue_labels)

    bonds_indexes_gen = gen_atom_pair_connectivity_masks(
        residue_labels_ij_gen,
        residue_numbers_ij_gen,
        atom_labels_ij_gen,
        bonds_intra,
        bonds_inter,
        )

    num_ij_pairs = len(atom_labels) * (len(atom_labels) - 1) // 2
    other_bool = not base_bool
    bonds_mask = np.full(num_ij_pairs, base_bool)

    for idx in bonds_indexes_gen:
        bonds_mask[idx] = other_bool

    return bonds_mask


def create_sidechains_masks_per_residue(
        residue_numbers,
        atom_labels,
        backbone_atoms,
        ):
    """
    Create a map of numeric indexing masks pointing to side chains atoms.

    Create separate masks per residue.

    Parameters
    ----------
    residue_numbers : np.ndarray, shape (N,)
        The atom residue numbers of the protein.

    atom_labels : np.ndarray, shape (N,)
        The atom labels of the protein.

    backbone_atoms : list or tuple
        The labels of all possible backbone atoms.

    Returns
    -------
    list of tuples of length 2
        List indexes refer to protein residues, index 0 is residue 1.
        Per residue, a tuple of length 2 is given. Tuple index 0 are
        the indexes of that residue sidechain atoms mapped to an array
        of the `atom_labels` and `residue_numbers` characteristics.
        The tuple index 1 is an array of length M, where M is the number
        of sidechain atoms for that residue, defaults to np.nan.
    """
    assert residue_numbers.size == atom_labels.size
    assert type(backbone_atoms) in (list, tuple, np.ndarray)

    ss = []
    ssa = ss.append

    is_backbone = np.isin(atom_labels, backbone_atoms)
    is_sidechain = np.logical_not(is_backbone)

    float64 = np.float64
    full = np.full
    logical_and = np.logical_and
    nan = np.nan
    where = np.where

    # sorted here is mandatory because ss indexes must follow the residue
    # numbering
    for resnum in sorted(list(set(residue_numbers))):
        is_residue = residue_numbers == resnum
        bool_result = logical_and(is_residue, is_sidechain)
        ssa((
            where(bool_result)[0],
            full((np.sum(bool_result), 3), nan, dtype=float64)
            ))

    return ss


def gen_3l_residue_labels_per_atom(
        input_seq_3letter,
        atom_labels,
        ):
    """
    Generate residue 3-letter labels per atom.

    Parameters
    ----------
    input_seq_3letter : list of 3letter residue codes
        Most not be a generator.

    atom_labels : list or tuple of atom labels
        Most not be a generator.

    Yields
    ------
    String of length 3
         The 3-letter residue code per atom.
    """
    _count_labels = Counter(atom_labels)['N']
    _len = len(input_seq_3letter)
    assert _count_labels == _len, (_count_labels, _len)

    counter = -1
    for atom in atom_labels:
        if atom == 'N':
            counter += 1
        yield input_seq_3letter[counter]


def gen_atom_pair_connectivity_masks(
        res_names_ij,
        res_num_ij,
        atom_names_ij,
        connectivity_intra,
        connectivity_inter,
        ):
    """
    Generate atom pair connectivity indexes.

    Given atom information for the ij pairs and connectivity criteria,
    yields the index of the ij pair if the pair is connected according
    to the connectivity criteria.

    For example, if the ij pair is covalently bonded, or 3 bonds apart,
    etc.

    Parameters
    ----------
    res_names_ij
    res_num_ij,
    atom_names_ij, iterables of the same length and synchronized information.

    connectivity_intra,
    connectivity_inter, dictionaries mapping atom labels connectivity

    Depends
    -------
    `are_connected`
    """
    zipit = zip(res_names_ij, res_num_ij, atom_names_ij)
    counter = 0
    for (rn1, _), (n1, n2), (a1, a2) in zipit:

        found_connectivity = are_connected(
            int(n1),
            int(n2),
            rn1,
            a1,
            a2,
            connectivity_intra,
            connectivity_inter,
            )

        if found_connectivity:
            yield counter

        counter += 1


def gen_ij_pairs_upper_diagonal(data):
    """
    Generate upper diagonal ij pairs in tuples.

    The diagonal is not considered.

    Yields
    ------
    tuple of length 2
        IJ pairs in the form of N*(N-1) / 2.
    """
    for i in range(len(data) - 1):
        for j in range(i + 1, len(data)):
            yield (data[i], data[j])


def gen_PDB_from_conformer(
        input_seq_3_letters,
        atom_labels,
        residues,
        coords,
        ALF=atom_line_formatter,
        AA1TO3=aa1to3,
        ROUND=np.round,
        ):
    """."""
    lines = []
    LINES_APPEND = lines.append
    ALF_FORMAT = ALF.format
    resi = -1

    # this is possible ONLY because there are no DOUBLE CHARS atoms
    # in the atoms that constitute a protein chain
    ATOM_LABEL_FMT = ' {: <3}'.format

    assert len(atom_labels) == coords.shape[0]

    for i in range(len(atom_labels)):

        if np.isnan(coords[i, 0]):
            continue

        if atom_labels[i] == 'N':
            resi += 1
            current_residue = input_seq_3_letters[resi]
            current_resnum = residues[i]

        atm = atom_labels[i].strip()
        ele = atm.lstrip('123')[0]

        if len(atm) < 4:
            atm = ATOM_LABEL_FMT(atm)

        LINES_APPEND(ALF_FORMAT(
            'ATOM',
            i,
            atm,
            '',
            current_residue,
            'A',
            current_resnum,
            '',
            coords[i, 0],
            coords[i, 1],
            coords[i, 2],
            0.0,
            0.0,
            '',
            ele,
            '',
            ))

    return '\n'.join(lines)


def gen_residue_number_per_atom(atom_labels, start=1):
    """
    Create a list of residue numbers based on atom labels.

    This is a contextualized function, not an abstracted one.
    Considers `N` to be the first atom of the residue.

    Yields
    ------
    ints
        The integer residue number per atom label.
    """
    assert atom_labels[0] == 'N', atom_labels[0]

    # creates a seamless interface between human and python 0-indexes
    start -= 1
    for al in atom_labels:
        if al == 'N':
            start += 1
        yield start


def make_list_atom_labels(input_seq, atom_labels_dictionary):
    """
    Make a list of the atom labels for an `input_seq`.

    Considers the N-terminal to be protonated H1 to H3.

    Parameters
    ----------
    input_seq : str
        1-letter amino-acid sequence.

    atom_labels_dictionary : dict
        The ORDERED atom labels per residue.

    Returns
    -------
    list
        List of consecutive atom labels for the protein.
    """
    labels = []
    LE = labels.extend

    first_residue_atoms = atom_labels_dictionary[input_seq[0]]

    # the first residue is a special case, we add here the three protons
    # for consistency with the forcefield
    # TODO: parametrize? multiple forcefields?
    for atom in first_residue_atoms:
        if atom == 'H':
            LE(('H1', 'H2', 'H3'))
        else:
            labels.append(atom)

    for residue in input_seq[1:]:
        LE(atom_labels_dictionary[residue])

    labels.append('OXT')

    assert Counter(labels)['N'] == len(input_seq)
    assert labels[-1] == 'OXT'
    assert 'H1' in labels
    assert 'H2' in labels
    assert 'H3' in labels
    return labels


def populate_ff_parameters_in_structure(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        ):
    """
    Creates a list of SIGMAS, EPSILONS, and CHARGES.

    Matches label information in the conformer to the respective values
    in the force field.
    """
    assert len(atom_labels) == len(residue_numbers) == len(residue_labels)

    sigmas_l = []
    epsilons_l = []
    charges_l = []

    sigmas_append = sigmas_l.append
    epsilons_append = epsilons_l.append
    charges_append = charges_l.append

    zipit = zip(atom_labels, residue_numbers, residue_labels)
    for atom_name, res_num, res_label in zipit:

        # adds C to the terminal residues
        if res_num == residue_numbers[-1]:
            res = 'C' + res_label
            was_in_C_terminal = True
            assert res.isupper() and len(res) == 4, res

        elif res_num == residue_numbers[0]:
            res = 'N' + res_label
            was_in_N_terminal = True
            assert res.isupper() and len(res) == 4, res

        else:
            res = res_label

        # TODO:
        # define protonation state in parameters
        if res_label.endswith('HIS'):
            res_label = res_label[:-3] + 'HIP'

        try:
            # force field atom type
            atype = force_field[res][atom_name]['type']

        # TODO:
        # try/catch is here to avoid problems with His...
        # for this purpose we are only using side-chains
        except KeyError:
            raise KeyError(tuple(force_field[res].keys()))

        ep = float(force_field[atype]['epsilon'])
        # epsilons by order of atom
        epsilons_append(ep)

        sig = float(force_field[atype]['sigma'])
        # sigmas by order of atom
        sigmas_append(sig)

        charge = float(force_field[res][atom_name]['charge'])
        charges_append(charge)

    assert len(epsilons_l) == len(sigmas_l) == len(charges_l), 'Lengths differ.'
    assert len(epsilons_l) == len(atom_labels),\
        'Lengths differ. There should be one epsilon per atom.'
    assert was_in_C_terminal, \
        'The C terminal residue was never computed. It should have.'
    assert was_in_N_terminal, \
        'The N terminal residue was never computed. It should have.'

    # brief theoretical review:
    # sigmas and epsilons are combined parameters self vs self (i vs i)
    # charges refer only to the self alone (i)
    sigmas_ii = np.array(sigmas_l)
    epsilons_ii = np.array(epsilons_l)
    charges_i = np.array(charges_l)

    assert epsilons_ii.shape == (len(epsilons_l),)
    assert sigmas_ii.shape == (len(sigmas_l),)
    assert charges_i.shape == (len(charges_l),)

    return sigmas_ii, epsilons_ii, charges_i


def are_connected(n1, n2, rn1, a1, a2, bonds_intra, bonds_inter):
    """
    Detect if a certain atom pair is bonded accordind to criteria.

    Considers only to the self residue and next residue
    """
    # requires
    assert isinstance(n1, int) and isinstance(n2, int), (type(n1), type(n2))
    assert all(isinstance(i, str) for i in (rn1, a1, a2)), \
        (type(i) for i in (rn1, a1, a2))
    assert all(isinstance(i, dict) for i in (bonds_intra, bonds_inter)), \
        (type(i) for i in (bonds_intra, bonds_inter))

    answer = (
        (n1 == n2 and a2 in bonds_intra[rn1][a1])
        or (
            n1 + 1 == n2
            and (
                a1 in bonds_inter  # this void KeyError
                and a2 in bonds_inter[a1]
                )
            )
        )

    assert isinstance(answer, bool)
    return answer


# NOT USED HERE
# kept to maintain compatibility with cli_validate.
# TODO: revisit cli_validate
def generate_vdW_data(
        atom_labels,
        residue_numbers,
        residue_labels,
        vdW_radii,
        bonds_apart=3,
        tolerance=0.4,
        ):
    # {{{
    """
    Generate van der Waals related data structures.

    Generated data structures are aligned to the output generated by
    scipy.spatial.distances.pdist used to compute all pairs distances.

    This function generates an (N,) array with the vdW thresholds aligned
    with the `pdist` output. And, a (N,) boolean array where `False`
    denotes pairs of covalently bonded atoms.

    Input should satisfy:

    >>> len(atom_labels) == len(residue_numbers) == len(residue_labels)

    Parameters
    ----------
    atom_labels : array or list
        The ordered list of atom labels of the target conformer upon
        which the returned values of this function will be used.

    residue_numbers : array or list
        The list of residue numbers corresponding to the atoms of the
        target conformer.

    residue_labels : array or list
        The list of residue 3-letter labels of each atom of the target
        conformer.

    vdW_radii : dict
        A dictionary containing the van der Waals radii for each atom
        type (element) present in `atom_labels`.

    bonds_apart : int
        The number of bonds apart to ignore vdW clash validation.
        For example, 3 means vdW validation will only be computed for
        atoms at least 4 bonds apart.

    tolerance : float
        The tolerance in Angstroms.

    Returns
    -------
    nd.array, dtype=np.float
        The vdW atom pairs thresholds. Equals to the sum of atom vdW
        radius.

    nd.array, dtype=boolean
        `True` where the distance between pairs must be considered.
        `False` where pairs are covalently bound.

    Note
    ----
    This function is slow, in the order of 1 or 2 seconds, but it is
    executed only once at the beginning of the building protocol.
    """
    # }}}
    assert len(atom_labels) == len(residue_numbers)
    assert len(atom_labels) == len(residue_labels)

    # we need to use re because hydrogen atoms start with integers some times
    atoms_char = re.compile(r'[CHNOPS]')
    findall = atoms_char.findall
    cov_topologies = generate_residue_template_topology(
        amber_pdbs,
        atom_labels_amber,
        add_OXT=True,
        add_Nterminal_H=True,
        )
    bond_structure_local = \
        expand_topology_bonds_apart(cov_topologies, bonds_apart)
    inter_connect_local = inter_residue_connectivities[bonds_apart]

    # adds OXT to the bonds connectivity, so it is included in the
    # final for loop creating masks
    #add_OXT_to_residue(bond_structure_local[residue_labels[-1]])

    # the following arrays/lists prepare the conformer label data in agreement
    # with the function scipy.spatial.distances.pdist, in order for masks to be
    # used properly
    #
    # creates atom label pairs
    atoms = np.array([
        (a, b)
        for i, a in enumerate(atom_labels, start=1)
        for b in atom_labels[i:]
        ])

    # creates vdW radii pairs according to atom labels
    vdw_pairs = np.array([
        (vdW_radii[findall(a)[0]], vdW_radii[findall(b)[0]])
        for a, b in atoms
        ])

    # creats the vdW threshold according to vdW pairs
    # the threshold is the sum of the vdW radii pair
    vdW_sums = np.power(np.sum(vdw_pairs, axis=1) - tolerance, 2)

    # creates pairs for residue numbers, so we know from which residue number
    # is each atom of the confronted pair
    res_nums_pairs = np.array([
        (a, b)
        for i, a in enumerate(residue_numbers, start=1)
        for b in residue_numbers[i:]
        ])

    # does the same as above but for residue 3-letter labels
    res_labels_pairs = (
        (_a, _b)
        for i, _a in enumerate(residue_labels, start=1)
        for _b in residue_labels[i:]
        )

    # we want to compute clashes to all atoms that are not covalentely bond
    # that that have not certain connectivity between them
    # first we assum True to all cases, and we replace to False where we find
    # a bond connection
    vdW_non_bond = np.full(len(atoms), True)

    counter = -1
    zipit = zip(atoms, res_nums_pairs, res_labels_pairs)
    for (a1, a2), (res1, res2), (l1, _) in zipit:
        counter += 1
        if res1 == res2 and a2 in bond_structure_local[l1][a1]:
            # here we found a bond connection
            vdW_non_bond[counter] = False

        # handling the special case of a covalent bond between two atoms
        # not part of the same residue
        elif a1 in inter_connect_local \
                and res2 == res1 + 1 \
                and a2 in inter_connect_local[a1]:
            vdW_non_bond[counter] = False

    return vdW_sums, vdW_non_bond


def calc_energy(coords, acoeff, bcoeff, charges_ij, bonds_ge_3_mask):
    CALC_DISTS = calc_all_vs_all_dists
    NANSUM = np.nansum
    distances_ij = CALC_DISTS(coords)

    energy_lj = njit_calc_LJ_energy(
        acoeff,
        bcoeff,
        distances_ij,
        bonds_ge_3_mask,
        )

    energy_elec_ij = charges_ij / distances_ij
    energy_elec = NANSUM(energy_elec_ij[bonds_ge_3_mask])

    # total energy
    return energy_lj + energy_elec


def create_energy_func_params(atom_labels, residue_numbers, residue_labels):
    # /
    # Prepares terms for the energy function
    # TODO: parametrize this.
    # the user should be able to chose different forcefields
    ff14SB = read_ff14SB_params()
    res_topology = generate_residue_template_topology(
        amber_pdbs,
        atom_labels_amber,
        add_OXT=True,
        add_Nterminal_H=True,
        )
    bonds_equal_3_intra = topology_3_bonds_apart(res_topology)
    bonds_le_2_intra = expand_topology_bonds_apart(res_topology, 2)

    # units in:
    # nm, kJ/mol, proton units
    sigmas_ii, epsilons_ii, charges_i = populate_ff_parameters_in_structure(
        atom_labels,
        residue_numbers,
        residue_labels,
        ff14SB,  # the forcefield dictionary
        )

    # this mask will deactivate calculations in covalently bond atoms and
    # atoms separated 2 bonds apart
    _bonds_ge_3_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        bonds_le_2_intra,
        bonds_le_2_inter,
        base_bool=True,
        )
    # on lun 21 dic 2020 17:33:31 EST, I tested for 1M sized array
    # the numeric indexing performed better than the boolean indexing
    # 25 ns versus 31 ns.
    bonds_ge_3_mask = np.where(_bonds_ge_3_mask)[0]

    # /
    # Prepares Coulomb and Lennard-Jones pre computed parameters:
    # calculates ij combinations using raw njitted functions because using
    # numpy outer variantes in very large systems overloads memory and
    # reduces performance.
    #
    num_ij_pairs = len(atom_labels) * (len(atom_labels) - 1) // 2
    # sigmas
    sigmas_ij_pre = np.empty(num_ij_pairs, dtype=np.float64)
    njit_calc_sum_upper_diagonal_raw(sigmas_ii, sigmas_ij_pre)
    #
    # epsilons
    epsilons_ij_pre = np.empty(num_ij_pairs, dtype=np.float64)
    njit_calc_multiplication_upper_diagonal_raw(epsilons_ii, epsilons_ij_pre)
    #
    # charges
    charges_ij = np.empty(num_ij_pairs, dtype=np.float64)
    njit_calc_multiplication_upper_diagonal_raw(charges_i, charges_ij)

    # mixing rules
    epsilons_ij = epsilons_ij_pre ** 0.5
    sigmas_ij = sigmas_ij_pre * 5  # mixing + nm to Angstrom converstion

    acoeff = 4 * epsilons_ij * (sigmas_ij ** 12)
    bcoeff = 4 * epsilons_ij * (sigmas_ij ** 6)

    charges_ij *= 0.25  # dielectic constant

    # The mask to identify ij pairs exactly 3 bonds apart is needed for the
    # special scaling factor of Coulomb and LJ equations
    # This mask will be used only aftert the calculation of the CLJ params
    bonds_exact_3_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        bonds_equal_3_intra,
        bonds_equal_3_inter,
        )

    # this is the Lennard-Jones special case, where scaling factors are applied
    # to atoms bonded 3 bonds apart, known as the '14' cases.
    # 0.4 was calibrated manually, until I could find a conformer
    # within 50 trials dom 20 dic 2020 13:16:50 EST
    # I believe, because we are not doing local minimization here, we
    # cannot be that strick with the 14 scaling factor, and a reduction
    # factor of 2 is not enough
    acoeff[bonds_exact_3_mask] *= float(ff14SB['lj14scale']) * 0.2  # was 0.4
    bcoeff[bonds_exact_3_mask] *= float(ff14SB['lj14scale']) * 0.2
    charges_ij[bonds_exact_3_mask] *= float(ff14SB['coulomb14scale'])

    #del sigmas_ij_pre, epsilons_ij_pre, epsilons_ij, sigmas_ij
    #del bonds_exact_3_mask, _bonds_ge_3_mask, ff14SB
    assert len(acoeff) == len(bcoeff) == len(charges_ij)

    return acoeff, bcoeff, charges_ij, bonds_ge_3_mask
    # ?


if __name__ == "__main__":
    libcli.maincli(ap, main)
