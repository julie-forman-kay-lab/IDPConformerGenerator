"""
Builds IDP conformers.

Build from a database of torsion angles and secondary structure
information. Database is as created by `idpconfgen torsions` CLI.

USAGE:
    $ idpconfgen build -db torsions.json -seq MMMMMMM...

"""
import argparse
import sys
from functools import partial
from itertools import cycle
from multiprocessing import Pool, Queue
# from numbers import Number
#from random import choice as randchoice
from random import randint
from time import time

import numpy as np
from numba import njit

from idpconfgen import Path, log
from idpconfgen.core.build_definitions import (
    backbone_atoms,
    build_bend_H_N_C,
    distance_C_O,
    distance_H_N,
    forcefields,
    n_terminal_h_coords_at_origin,
    sidechain_templates,
    )
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs import libcli
from idpconfgen.libs.libbuild import (
    prepare_slice_dict,
    read_db_to_slices_single_secondary_structure,
    compute_sidechains,
    create_sidechains_masks_per_residue,
    get_cycle_bond_type,
    get_cycle_distances_backbone,
    init_conflabels,
    init_confmasks,
    prepare_energy_function,
    read_db_to_slices,
    )
from idpconfgen.libs.libcalc import (
    calc_residue_num_from_index,
    calc_torsion_angles,
    make_coord_Q,
    make_coord_Q_COO,
    make_coord_Q_planar,
    make_seq_probabilities,
    place_sidechain_template,
    rotate_coordinates_Q_njit,
    rrd10_njit,
    )
from idpconfgen.libs.libhigherlevel import bgeo_reduce
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libparse import (
    get_seq_chunk_njit,
    get_seq_next_residue_njit,
    get_trimer_seq_njit,
    remap_sequence,
    translate_seq_to_3l,
    )
from idpconfgen.libs.libpdb import atom_line_formatter
from idpconfgen.logger import S, T, init_clean_files, pre_msg, report_on_crash


_file = Path(__file__).myparents()
LOGFILESNAME = 'idpconfgen_build'

# Global variables needed to build conformers.
# Why are global variables needed?
# I use global variables to facilitate distributing conformer creation
# processes across multiple cores. In this way cores can read global variables
# fast and with non-significant overhead.

# Bond Geometry library variables
# if __name__ == '__main__', these will be populated in main()
# else will be populated in conformer_generator
# populate_globals() populates these variables once called.
BGEO_path = Path(_file, 'core', 'data', 'bgeo.tar')
BGEO_full = {}
BGEO_trimer = {}
BGEO_res = {}

# SLICES and ANGLES will be populated in main() with the torsion angles.
# it is not expected SLICES or ANGLES to be populated anywhere else.
# The slice objects from where the builder will feed to extract torsion
# chunks from ANGLES.
ANGLES = None
SLICES = []
SLICEDICT_MONOMERS = None
SLICEDICT_XMERS = None
XMERPROBS = None
GET_ADJ = None

# keeps a record of the conformer numbers written to disk across the different
# cores
CONF_NUMBER = Queue()
RANDOMSEEDS = Queue()

# The conformer building process needs data structures for two different
# identities: the all-atom representation of the input sequence, and the
# corresponding Ala/Gly/Pro template uppon which the coordinates will be built.
# These variables are defined at the module level so they serve as global
# variables to be read by the different process during multiprocessing. Reading
# from global variables is performant in Python multiprocessing. This is the
# same strategy as applied for SLICES and ANGLES.
ALL_ATOM_LABELS = None
ALL_ATOM_MASKS = None
ALL_ATOM_EFUNC = None
TEMPLATE_LABELS = None
TEMPLATE_MASKS = None
TEMPLATE_EFUNC = None


def are_globals():
    """Assess if global variables needed for building are populated."""
    return all((
        ALL_ATOM_LABELS,
        ALL_ATOM_MASKS,
        ALL_ATOM_EFUNC,
        TEMPLATE_LABELS,
        TEMPLATE_MASKS,
        TEMPLATE_EFUNC,
        BGEO_full,
        BGEO_trimer,
        BGEO_res,
        ))


# CLI argument parser parameters
_name = 'build'
_help = 'Builds conformers from database.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)


ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
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
    help='The Conformer residue sequence. String or FASTA file.',
    required=True,
    nargs='?',
    action=libcli.SeqOrFasta,
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

_ffchoice = list(forcefields.keys())
ap.add_argument(
    '-ff',
    '--forcefield',
    help=(
        'Forcefield parameters and atom labels. '
        f'Defaults to {_ffchoice[0]}.'
        ),
    choices=_ffchoice,
    default=_ffchoice[0],
    )

ap.add_argument(
    '-bgeo_path',
    '--bgeo_path',
    help=(
        'Path to the bond geometry database as generated by `bgeo` CLI .'
        'Defaults to `None`, uses the internal library.'
        ),
    default=None,
    )

ap.add_argument(
    '-et',
    '--energy-threshold',
    help='The energy threshold after which conformers/chunks will be rejected',
    default=10.0,
    type=float,
    )

ap.add_argument(
    '-xp',
    '--xmers-probs',
    help='The Xmers relative probabilities for selection',
    default=False,
    nargs='+',
    action=libcli.ListOfInts,
    )


libcli.add_argument_random_seed(ap)
libcli.add_argument_ncores(ap)


def main(
        input_seq,
        database,
        dssp_regexes=r'L+',#r'(?=(L{2,6}))',
        func=None,
        forcefield=None,
        bgeo_path=None,
        nconfs=1,
        ncores=1,
        random_seed=0,
        xmers_probs,
        **kwargs,  # other kwargs target energy function, for example.
        ):
    """
    Execute main client logic.

    Distributes over processors.
    """
    init_clean_files(log, LOGFILESNAME)
    log.info(f'input sequence: {input_seq}')

    # Calculates how many conformers are built per core
    if nconfs < ncores:
        ncores = 1
        conformers_per_core = nconfs
        remaining_confs = 0
    else:
        conformers_per_core = nconfs // ncores
        # in case nconfs is not multiple of ncores, builds the remaining confs
        # at the end
        remaining_confs = nconfs % ncores

    log.info(
        f'running in {ncores} cores with '
        f'{remaining_confs} remaining confs'
        )

    # populates globals
    #if False:
    #    global ANGLES, SLICES
    #    _slices, ANGLES = read_db_to_slices(database, dssp_regexes, ncores=ncores)
    #    SLICES.extend(_slices)
    # #

    # we use a dictionary because chunks will be evaluated to exact match
    global ANGLES, SLICEDICT_MONOMERS, SLICEDICT_XMERS, XMERPROBS, GET_ADJ
    primary, ANGLES = read_db_to_slices_single_secondary_structure(database, dssp_regexes)
    SLICEDICT_XMERS = prepare_slice_dict(primary, input_seq)
    SLICEDICT_MONOMERS = SLICEDICT_XMERS.pop(1)

    XMERPROBS = make_seq_probabilities(
        xmers_probs if xmers_probs else list(SLICEDICT_XMERS.keys()),
        reverse=True,
        )

    GET_ADJ = get_adjacent_angles(
        list(SLICEDICT_XMERS.keys()),
        XMERPROBS,
        input_seq,
        ANGLES,
        SLICEDICT_XMERS,
        SLICEDICT_MONOMERS,
        )

    populate_globals(
        input_seq=input_seq,
        bgeo_path=bgeo_path or BGEO_path,
        forcefield=forcefields[forcefield],
        **kwargs)

    # create different random seeds for the different cores
    # seeds created to the cores based on main seed are predictable
    for i in range(ncores + bool(remaining_confs)):
        RANDOMSEEDS.put(random_seed + i)

    # creates a queue of numbers that will serve all subprocesses.
    # Used to name the output files, conformer_1, conformer_2, ...
    for i in range(1, nconfs + 1):
        CONF_NUMBER.put(i)

    # prepars execution function
    consume = partial(
        _build_conformers,
        input_seq=input_seq,  # string
        nconfs=conformers_per_core,  # int
        **kwargs,
        )

    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_prefix=_name,
        )

    start = time()
    with Pool(ncores) as pool:
        imap = pool.imap(execute, range(ncores))
        for _ in imap:
            pass

    if remaining_confs:
        execute(conformers_per_core * ncores, nconfs=remaining_confs)

    log.info(f'{nconfs} conformers built in {time() - start:.3f} seconds')


def populate_globals(
        *,
        input_seq=None,
        bgeo_path=BGEO_path,
        forcefield=None,
        **efunc_kwargs):
    """
    Populate global variables needed for building.

    Currently, global variables include:

    BGEO_full
    BGEO_trimer
    BGEO_res
    ALL_ATOM_LABELS, ALL_ATOM_MASKS, ALL_ATOM_EFUNC
    TEMPLATE_LABELS, TEMPLATE_MASKS, TEMPLATE_EFUNC

    Parameters
    ----------
    bgeo_path : str or Path
        The path pointing to a bond geometry library as created by the
        `bgeo` CLI.

    forcefield : str
        A key in the `core.build_definitions.forcefields` dictionary.
    """
    if not isinstance(input_seq, str):
        raise ValueError(
            '`input_seq` not valid. '
            f'Expected string found {type(input_seq)}'
            )

    global BGEO_full, BGEO_trimer, BGEO_res

    BGEO_full.update(read_dictionary_from_disk(bgeo_path))
    _1, _2 = bgeo_reduce(BGEO_full)
    BGEO_trimer.update(_1)
    BGEO_res.update(_2)
    del _1, _2
    assert BGEO_full
    assert BGEO_trimer
    assert BGEO_res
    # this asserts only the first layer of keys
    assert list(BGEO_full.keys()) == list(BGEO_trimer.keys()) == list(BGEO_res.keys())  # noqa: E501

    # populates the labels
    global ALL_ATOM_LABELS, ALL_ATOM_MASKS, ALL_ATOM_EFUNC
    global TEMPLATE_LABELS, TEMPLATE_MASKS, TEMPLATE_EFUNC

    topobj = forcefield(add_OXT=True, add_Nterminal_H=True)

    ALL_ATOM_LABELS = init_conflabels(input_seq, topobj.atom_names)
    TEMPLATE_LABELS = init_conflabels(remap_sequence(input_seq), topobj.atom_names)  # noqa: E501

    ALL_ATOM_MASKS = init_confmasks(ALL_ATOM_LABELS.atom_labels)
    TEMPLATE_MASKS = init_confmasks(TEMPLATE_LABELS.atom_labels)

    ALL_ATOM_EFUNC = prepare_energy_function(
        ALL_ATOM_LABELS.atom_labels,
        ALL_ATOM_LABELS.res_nums,
        ALL_ATOM_LABELS.res_labels,
        topobj,
        **efunc_kwargs)

    TEMPLATE_EFUNC = prepare_energy_function(
        TEMPLATE_LABELS.atom_labels,
        TEMPLATE_LABELS.res_nums,
        TEMPLATE_LABELS.res_labels,
        topobj,
        **efunc_kwargs)

    del topobj
    return


# private function because it depends on the global `CONF_NUMBER`
# which is assembled in `main()`
def _build_conformers(
        *args,
        input_seq=None,
        conformer_name='conformer',
        nconfs=1,
        **kwargs,
        ):
    """Arrange building of conformers and saves them to PDB files."""
    ROUND = np.round

    # TODO: this has to be parametrized for the different HIS types
    input_seq_3_letters = translate_seq_to_3l(input_seq)

    builder = conformer_generator(
        input_seq=input_seq,
        random_seed=RANDOMSEEDS.get(),
        **kwargs)

    atom_labels, residue_numbers, residue_labels = next(builder)

    for _ in range(nconfs):

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
        energy_threshold=10,
        bgeo_path=None,
        forcefield=None,
        lj_term=True,
        coulomb_term=False,
        random_seed=0,
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

    sidechain_method : str
        The method used to build/pack sidechains over the backbone
        structure. Defaults to `faspr`.
        Expects a key in `libs.libbuild.compute_sidechains`.

    bgeo_path : str of Path
        Path to a bond geometry library as created by `bgeo` CLI.

    lj_term : bool
        Whether to compute the Lennard-Jones term during building and
        validation. If false, expect a physically meaningless result.
    """
    if not isinstance(input_seq, str):
        raise ValueError(f'`input_seq` must be given! {input_seq}')
    if sidechain_method not in compute_sidechains:
        raise ValueError(
            f'{sidechain_method} not in `compute_sidechains`. '
            f'Expected {list(compute_sidechains.keys())}.'
            )

    log.info(f'random seed: {random_seed}')
    np.random.seed(random_seed)
    seed_report = pre_msg(f'seed {random_seed}')

    # prepares protein sequences
    all_atom_input_seq = input_seq
    template_input_seq = remap_sequence(all_atom_input_seq)
    template_seq_3l = translate_seq_to_3l(template_input_seq)
    del input_seq

    BUILD_BEND_H_N_C = build_bend_H_N_C
    CALC_TORSION_ANGLES = calc_torsion_angles
    DISTANCE_NH = distance_H_N
    DISTANCE_C_O = distance_C_O
    ISNAN = np.isnan
    GET_TRIMER_SEQ = get_trimer_seq_njit
    MAKE_COORD_Q_COO_LOCAL = make_coord_Q_COO
    MAKE_COORD_Q_PLANAR = make_coord_Q_planar
    MAKE_COORD_Q_LOCAL = make_coord_Q
    NAN = np.nan
    NORM = np.linalg.norm
    N_TERMINAL_H = n_terminal_h_coords_at_origin
    PI2 = np.pi * 2
    PLACE_SIDECHAIN_TEMPLATE = place_sidechain_template
    RAD_60 = np.radians(60)
    RC = np.random.choice
    RINT = randint
    ROT_COORDINATES = rotate_coordinates_Q_njit
    RRD10 = rrd10_njit
    SIDECHAIN_TEMPLATES = sidechain_templates
    angles = ANGLES
    slices = SLICES
    global BGEO_full
    global BGEO_trimer
    global BGEO_res
    global ALL_ATOM_LABELS
    global ALL_ATOM_MASKS
    global ALL_ATOM_EFUNC
    global TEMPLATE_LABELS
    global TEMPLATE_MASKS
    global TEMPLATE_EFUNC
    global XMERPROBS
    global SLICEDICT_MONOMERS
    global SLICEDICT_XMERS
    global GET_ADJ

    # these flags exist to populate the global variables in case they were not
    # populated yet. Global variables are populated through the main() function
    # if the script runs as CLI. Otherwise, if conformer_generator() is imported
    # and used directly, the global variables need to be configured here.
    if not are_globals():
        if forcefield not in forcefields:
            raise ValueError(
                f'{forcefield} not in `forcefields`. '
                f'Expected {list(forcefields.keys())}.'
                )
        populate_globals(
            input_seq=all_atom_input_seq,
            bgeo_path=bgeo_path or BGEO_path,
            forcefield=forcefields[forcefield],
            lj_term=lj_term,
            coulomb_term=coulomb_term,
            )

    # semantic exchange for speed al readibility
    with_sidechains = not(disable_sidechains)

    if with_sidechains:
        build_sidechains = compute_sidechains[sidechain_method](all_atom_input_seq)  # noqa: E501

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

    # yields atom labels
    # all conformers generated will share these labels
    yield (
        ALL_ATOM_LABELS.atom_labels,
        ALL_ATOM_LABELS.res_nums,
        ALL_ATOM_LABELS.res_labels,
        )

    all_atom_num_atoms = len(ALL_ATOM_LABELS.atom_labels)
    template_num_atoms = len(TEMPLATE_LABELS.atom_labels)

    all_atom_coords = np.full((all_atom_num_atoms, 3), NAN, dtype=np.float64)
    template_coords = np.full((template_num_atoms, 3), NAN, dtype=np.float64)

    # +2 because of the dummy coordinates required to start building.
    # see later adding dummy coordinates to the structure seed
    bb = np.full((TEMPLATE_MASKS.bb3.size + 2, 3), NAN, dtype=np.float64)
    bb_real = bb[2:, :]  # backbone coordinates without the dummies

    # coordinates for the carbonyl oxigen atoms
    bb_CO = np.full((TEMPLATE_MASKS.COs.size, 3), NAN, dtype=np.float64)

    # notice that NHydrogen_mask does not see Prolines
    bb_NH = np.full((TEMPLATE_MASKS.NHs.size, 3), NAN, dtype=np.float64)
    bb_NH_idx = np.arange(len(bb_NH))
    # Creates masks and indexes for the `for` loop used to place NHs.
    # The first residue has no NH, prolines have no NH.
    non_pro = np.array(list(template_input_seq)[1:]) != 'P'
    # NHs index numbers in bb_real
    bb_NH_nums = np.arange(3, (len(template_input_seq) - 1) * 3 + 1, 3)[non_pro]
    bb_NH_nums_p1 = bb_NH_nums + 1
    assert bb_NH.shape[0] == bb_NH_nums.size == bb_NH_idx.size

    # sidechain masks
    # this is sidechain agnostic, works for every sidechain, yet here we
    # use only ALA, PRO, GLY - Mon Feb 15 17:29:20 2021
    ss_masks = create_sidechains_masks_per_residue(
        TEMPLATE_LABELS.res_nums,
        TEMPLATE_LABELS.atom_labels,
        backbone_atoms,
        )
    # ?

    # /
    # creates seed coordinates:
    # because the first torsion angle of a residue is the omega, we need
    # to prepare 2 dummy atoms to simulate the residue -1, so that the
    # first omega can be placed. There is no need to setup specific
    # positions, just to create a place upon which the build atom
    # routine can create a new atom from a torsion.
    dummy_CA_m1_coord = np.array((0.0, 1.0, 1.0))
    dummy_C_m1_coord = np.array((0.0, 1.0, 0.0))
    n_terminal_N_coord = np.array((0.0, 0.0, 0.0))

    # seed coordinates array
    seed_coords = np.array((
        dummy_CA_m1_coord,
        dummy_C_m1_coord,
        n_terminal_N_coord,
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
        bond_lens = get_cycle_distances_backbone()
        bond_type = get_cycle_bond_type()

        # in the first run of the loop this is unnecessary, but is better to
        # just do it once than flag it the whole time
        template_coords[:, :] = NAN
        bb[:, :] = NAN
        bb_CO[:, :] = NAN
        bb_NH[:, :] = NAN
        for _mask, _coords in ss_masks:
            _coords[:, :] = NAN

        bb[:3, :] = seed_coords  # this contains a dummy coord at position 0

        # add N-terminal hydrogens to the origin

        bbi = 1  # starts at 1 because there are two dummy atoms
        bbi0_R_CLEAR()
        bbi0_R_APPEND(bbi)

        COi = 0  # carbonyl atoms
        COi0_R_CLEAR()
        COi0_R_APPEND(COi)

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
            #print(bbi)

            # I decided to use an if-statement here instead of polymorph
            # the else clause to a `generative_function` variable because
            # the resulting overhead from the extra function call and
            # **kwargs handling was greater then the if-statement processing
            # https://pythonicthoughtssnippets.github.io/2020/10/21/PTS14-quick-in-if-vs-polymorphism.html
            if generative_function:
                agls = generative_function(
                    nres=RINT(1, 6),
                    cres=calc_residue_num_from_index(bbi)
                    )

            else:
                # following `aligndb` function,
                # `angls` will always be cyclic with:
                # omega - phi - psi - omega - phi - psi - (...)
                #agls = angles[RC(slices), :].ravel()
                # agls = angles[:, :].ravel()

                #primer_template = get_idx_primer_njit(
                #    all_atom_input_seq,
                #    RC(slices_dict_keys, p=XMERPROBS),
                #    calc_residue_num_from_index(bbi - 1),
                #    )

                # algorithm for adjacent building
                # TODO
                # primer_template here is used temporarily, and needs to be
                # removed when get_adj becomes an option
                primer_template, agls = GET_ADJ(bbi - 1)

            # index at the start of the current cycle
            #print(primer_template)
            PRIMER = cycle(primer_template)
            try:
                for (omg, phi, psi) in zip(agls[0::3], agls[1::3], agls[2::3]):

                    current_res_number = calc_residue_num_from_index(bbi - 1)

                    # assert the residue being built is of the same nature as the one in the angles
                    # TODO: remove this assert
                    assert all_atom_input_seq[current_res_number] == next(PRIMER)

                    curr_res, tpair = GET_TRIMER_SEQ(
                        all_atom_input_seq,
                        current_res_number,
                        )
                    torpair = f'{RRD10(phi)},{RRD10(psi)}'

                    for torsion_angle in (omg, phi, psi):

                        _bt = next(bond_type)

                        try:
                            _bend_angle = RC(BGEO_full[_bt][curr_res][tpair][torpair])  # noqa: E501
                        except KeyError:
                            try:
                                _bend_angle = RC(BGEO_trimer[_bt][curr_res][tpair])  # noqa: E501
                            except KeyError:
                                _bend_angle = RC(BGEO_res[_bt][curr_res])

                        _bond_lens = next(bond_lens)[curr_res]

                        bb_real[bbi, :] = MAKE_COORD_Q_LOCAL(
                            bb[bbi - 1, :],
                            bb[bbi, :],
                            bb[bbi + 1, :],
                            _bond_lens,
                            _bend_angle,
                            torsion_angle,
                            )
                        bbi += 1

                    try:
                        co_bend = RC(BGEO_full['Ca_C_O'][curr_res][tpair][torpair])  # noqa: E501
                    except KeyError:
                        try:
                            co_bend = RC(BGEO_trimer['Ca_C_O'][curr_res][tpair])
                        except KeyError:
                            co_bend = RC(BGEO_res['Ca_C_O'][curr_res])

                    bb_CO[COi, :] = MAKE_COORD_Q_PLANAR(
                        bb_real[bbi - 3, :],
                        bb_real[bbi - 2, :],
                        bb_real[bbi - 1, :],
                        distance=DISTANCE_C_O,
                        bend=co_bend
                        )
                    COi += 1

            except IndexError:
                # IndexError happens when the backbone is complete
                # in this protocol the last atom build was a carbonyl C
                # bbi is the last index of bb + 1, and the last index of
                # bb_real + 2

                # activate flag to finish loop at the end
                backbone_done = True

                # add the carboxyls
                template_coords[TEMPLATE_MASKS.cterm] = \
                    MAKE_COORD_Q_COO_LOCAL(bb[-2, :], bb[-1, :])

            # Adds N-H Hydrogens
            # Not a perfect loop. It repeats for Hs already placed.
            # However, was a simpler solution than matching the indexes
            # and the time cost is not a bottle neck.
            _ = ~ISNAN(bb_real[bb_NH_nums_p1, 0])
            for k, j in zip(bb_NH_nums[_], bb_NH_idx[_]):

                bb_NH[j, :] = MAKE_COORD_Q_PLANAR(
                    bb_real[k + 1, :],
                    bb_real[k, :],
                    bb_real[k - 1, :],
                    distance=DISTANCE_NH,
                    bend=BUILD_BEND_H_N_C,
                    )

            # Adds sidechain template structures
            for res_i in range(res_R[-1], current_res_number + 1):  # noqa: E501

                _sstemplate, _sidechain_idxs = \
                    SIDECHAIN_TEMPLATES[template_seq_3l[res_i]]

                sscoords = PLACE_SIDECHAIN_TEMPLATE(
                    bb_real[res_i * 3:res_i * 3 + 3, :],  # from N to C
                    _sstemplate,
                    )

                ss_masks[res_i][1][:, :] = sscoords[_sidechain_idxs]

            # Transfers coords to the main coord array
            for _smask, _sidecoords in ss_masks[:current_res_number + 1]:
                template_coords[_smask] = _sidecoords

            # / Place coordinates for energy calculation
            #
            # use `bb_real` to do not consider the initial dummy atom
            template_coords[TEMPLATE_MASKS.bb3] = bb_real
            template_coords[TEMPLATE_MASKS.COs] = bb_CO
            template_coords[TEMPLATE_MASKS.NHs] = bb_NH

            if len(bbi0_register) == 1:
                # places the N-terminal Hs only if it is the first
                # chunk being built
                _ = PLACE_SIDECHAIN_TEMPLATE(bb_real[0:3, :], N_TERMINAL_H)
                template_coords[TEMPLATE_MASKS.Hterm, :] = _[3:, :]
                current_Hterm_coords = _[3:, :]
                del _

                if template_input_seq[0] != 'G':
                    # rotates only if the first residue is not an
                    # alanie

                    # measure torsion angle reference H1 - HA
                    _h1_ha_angle = CALC_TORSION_ANGLES(
                        template_coords[TEMPLATE_MASKS.H1_N_CA_CB, :]
                        )[0]

                    # given any angle calculated along an axis, calculate how
                    # much to rotate along that axis to place the
                    # angle at 60 degrees
                    _rot_angle = _h1_ha_angle % PI2 - RAD_60

                    current_Hterm_coords = ROT_COORDINATES(
                        template_coords[TEMPLATE_MASKS.Hterm, :],
                        template_coords[1] / NORM(template_coords[1]),
                        _rot_angle,
                        )

                    template_coords[TEMPLATE_MASKS.Hterm, :] = current_Hterm_coords  # noqa: E501
            # ?

            total_energy = TEMPLATE_EFUNC(template_coords)
            #print(bbi, total_energy)

            if total_energy > energy_threshold:
                #print('---------- energy positive')
                # reset coordinates to the original value
                # before the last chunk added

                # reset the same chunk maximum 5 times,
                # after that reset also the chunk before
                try:
                    if number_of_trials > 30:
                        bbi0_R_POP()
                        COi0_R_POP()
                        res_R_POP()
                        number_of_trials = 0
                        number_of_trials2 += 1

                    if number_of_trials2 > 5:
                        bbi0_R_POP()
                        COi0_R_POP()
                        res_R_POP()
                        number_of_trials2 = 0
                        number_of_trials3 += 1

                    if number_of_trials3 > 5:
                        bbi0_R_POP()
                        COi0_R_POP()
                        res_R_POP()
                        number_of_trials3 = 0

                    _bbi0 = bbi0_register[-1]
                    _COi0 = COi0_register[-1]
                    _resi0 = res_R[-1]
                except IndexError:
                    # if this point is reached,
                    # we erased until the beginning of the conformer
                    # discard conformer, something went really wrong
                    broke_on_start_attempt = True
                    break  # conformer while loop, starts conformer from scratch

                # clean previously built protein chunk
                bb_real[_bbi0:bbi, :] = NAN
                bb_CO[_COi0:COi, :] = NAN

                # reset also indexes
                bbi = _bbi0
                COi = _COi0
                current_res_number = _resi0

                # coords needs to be reset because size of protein next
                # chunks may not be equal
                template_coords[:, :] = NAN
                template_coords[TEMPLATE_MASKS.Hterm, :] = current_Hterm_coords

                # prepares cycles for building process
                # this is required because the last chunk created may have been
                # the final part of the conformer
                if backbone_done:
                    bond_lens = get_cycle_distances_backbone()
                    bond_type = get_cycle_bond_type()

                # we do not know if the next chunk will finish the protein
                # or not
                backbone_done = False
                number_of_trials += 1
                continue  # send back to the CHUNK while loop

            # if the conformer is valid
            number_of_trials = 0
            bbi0_R_APPEND(bbi)
            COi0_R_APPEND(COi)
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
        all_atom_coords[ALL_ATOM_MASKS.bb4] = template_coords[TEMPLATE_MASKS.bb4]  # noqa: E501
        all_atom_coords[ALL_ATOM_MASKS.NHs] = template_coords[TEMPLATE_MASKS.NHs]  # noqa: E501
        all_atom_coords[ALL_ATOM_MASKS.Hterm] = template_coords[TEMPLATE_MASKS.Hterm]  # noqa: E501
        all_atom_coords[ALL_ATOM_MASKS.cterm, :] = template_coords[TEMPLATE_MASKS.cterm, :]  # noqa: E501

        if with_sidechains:

            all_atom_coords[ALL_ATOM_MASKS.non_Hs_non_OXT] = build_sidechains(
                template_coords[TEMPLATE_MASKS.bb4],
                )

            total_energy = ALL_ATOM_EFUNC(all_atom_coords)

            if total_energy > energy_threshold:
                log.info(seed_report(
                    f'Conformer with WORST energy {total_energy}'))
                continue
            else:
                log.info(seed_report(
                    f'finished conf: {conf_n} with energy {total_energy}'))

        #print(all_atom_coords)
        yield all_atom_coords
        conf_n += 1


def gen_PDB_from_conformer(
        input_seq_3_letters,
        atom_labels,
        residues,
        coords,
        ALF=atom_line_formatter,
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

    atom_i = 1
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
            atom_i,
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

        atom_i += 1

    return '\n'.join(lines)



def get_adjacent_angles(
        options,
        probs,
        seq,
        db,
        slice_dict,
        slicemonomers,
        RC=np.random.choice,
        ):
    """
    Get angles to build the next adjacent protein chunk.

    Parameters
    ----------
    options : list
        The length of the possible chunk sizes.

    probs : list
        A list with the relative probabilites to select from `options`.

    seq : str
        The conformer sequence.

    db : dict-like
        The angle omega/phi/psi database.

    slice_dict : dict-like
        A dictionary containing the chunks strings as keys and as values
        lists with slice objects.

    slicemonomers : dict-like
        The same as `slice_dict` but containing information only for
        monomers.
    """

    max_opt = max(options)

    def func(aidx):
        #print('Entering new angle collection ----------------------------', aidx)
        # calculates the current residue number from the atom index
        cr = calc_residue_num_from_index(aidx)

        # chooses the size of the chunk from pre-configured range of sizes
        plen = RC(options, p=probs)

        # defines the chunk identity accordingly
        primer_template = get_seq_chunk_njit(seq, cr, plen)
        next_residue = get_seq_next_residue_njit(seq, cr, plen)

        # recalculates the plen to avoid plen/template inconsistencies that
        # occur if the plen is higher then the number of
        # residues until the end of the protein.
        plen = len(primer_template)

        # Now the tricky loop to accommodate the very special case of prolines.
        template_not_found_once = False
        while 1 < plen <= max_opt:
            #print('L/PT/NR/Bool: ', plen, primer_template, next_residue, template_not_found_once)
            # Confirms the next residue to the template is NOT a proline
            if next_residue != 'P':
                try:
                    # attempts to retrieve angles to the template, exact sequence
                    # match
                    agls = db[RC(slice_dict[plen][primer_template]), :].ravel()
                    # if the angles are found, there is nothing more to do
                    break

                # in case the template is not present in the database. It may
                # occur for pentamers, perhaps some tetramers.
                except KeyError as err:
                    # reduces the size of the template by 1
                    # if primer_template[-1] == 'P':
                    #     plen -= 2
                    #     primer_template = primer_template[:-2]
                    # else:
                    #print('KeyError found, decreasing plen')
                    plen -= 1
                    next_residue = primer_template[-1]
                    primer_template = primer_template[:-1]
                    # if a certain primer template was not found, that means
                    # no larger primer will be found. Therefore, any change in
                    # plen or primer_template has to go towards reducing the
                    # length.
                    template_not_found_once = True
                    continue
            else:
                # if the template was not found, makes no sense to increase its
                # size, hence, we reduce it.
                #print('it has a proline')
                if template_not_found_once:
                    plen -= 1
                    next_residue = primer_template[-1]
                    primer_template = primer_template[:-1]
                # if we never actually attempted retrieving angles, so we
                # increase template size.
                else:
                    plen += 1
                    primer_template = get_seq_chunk_njit(seq, cr, plen)
                    next_residue = get_seq_next_residue_njit(seq, cr, plen)
        # beautiful while/else.
        # only executes if the template grew to large or too short.
        else:
            #print('entered the while else')
            # template grew too large, most likely because of the presence of
            # many consecutive prolines. Reduce the template to two residues
            # search for the pair considering the consecutive residue, but build
            # only the first residue in the pair. That is, the current residue
            if plen > max_opt:
                plen = 2
                primer_template = get_seq_chunk_njit(seq, cr, plen)
                agls = db[RC(slice_dict[plen][primer_template]), :].ravel()[:3]
                #print('sending just one', plen, primer_template)

            elif plen == 1:
                primer_template = seq[cr]
                # this should only happen at the end of the chain
                agls = db[RC(slicemonomers[primer_template]), :].ravel()
                #print('the conformer ended', cr)
            else:
                log.debug(plen, primer_template, seq[cr:cr + plen])
                # raise AssertionError to avoid `python -o` silencing
                raise AssertionError('The code should not arrive here')

        assert primer_template
        assert not np.any(np.isnan(agls))
        return primer_template, agls

    return func




if __name__ == "__main__":
    libcli.maincli(ap, main)
