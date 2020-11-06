"""
Builds IDP conformers.

Build from a database of torsion angles and secondary structure
information. Database is as created by `idpconfgen torsions` CLI.

USAGE:
    $ idpconfgen build -db torsions.json -seq MMMMMMM...

"""
import argparse
import re
from collections import Counter
from functools import partial
from itertools import cycle
from multiprocessing import Pool
from random import choice as randchoice
from random import randint
from time import time

import numpy as np

from idpconfgen import log
from idpconfgen.core.build_definitions import (
    add_OXT_to_residue,
    atom_labels,
    build_bend_angles_CA_C_Np1,
    build_bend_angles_Cm1_N_CA,
    build_bend_angles_N_CA_C,
    build_bend_H_N_C,
    distance_H_N,
    distances_C_Np1,
    distances_CA_C,
    distances_N_CA,
    expand_topology_bonds_apart,
    generate_residue_template_topology,
    inter_residue_connectivities,
    sidechain_templates,
    )
from idpconfgen.core.definitions import aa1to3, vdW_radii_dict
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs import libcli
from idpconfgen.libs.libcalc import (
    calc_all_vs_all_dists_square,
    make_coord_Q,
    make_coord_Q_CO,
    make_coord_Q_COO,
    place_sidechain_template,
    )
from idpconfgen.libs.libfilter import aligndb, regex_search
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libpdb import atom_line_formatter
from idpconfgen.libs.libtimer import timeme


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

#libcli.add_argument_vdWb(ap)
libcli.add_argument_vdWr(ap)
libcli.add_argument_vdWt(ap)
libcli.add_argument_ncores(ap)


SLICES = []
ANGLES = None


def main(
        input_seq,
        database,
        conformer_name='conformer',
        dssp_regexes=r'(?=(L{2,6}))',
        func=None,
        nconfs=1,
        ncores=1,
        vdW_bonds_apart=3,
        vdW_tolerance=0.4,
        vdW_radii='tsai1999',
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

    # prepars execution function
    execute = partial(
        main_exec,
        input_seq=input_seq,  # string
        nconfs=core_chunks,  # int
        conformer_name=conformer_name,  # string
        vdW_bonds_apart=vdW_bonds_apart,
        vdW_tolerance=vdW_tolerance,
        vdW_radii=vdW_radii,
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


def main_exec(
        execution_run,
        input_seq,
        conformer_name='conformer',
        generative_function=None,
        vdW_bonds_apart=3,
        vdW_tolerance=0.4,
        vdW_radii='tsai1999',
        nconfs=1,
        ):
    """
    Build conformers.

    *Note*: `execution_run` is a positional parameter in order to maintain
    operability with the multiprocessing operations in `main`, sorry for
    that :-)

    Parameters
    ----------
    execution_run : int
        A zero or positive integer. `execution_run` will number the
        conformers accordingly to the CPU job number which generates
        them.

        `execution_run` dictate the numbering suffix of the PDBs saved
        to disk.

        If a single CPU is being used, just use 0 for convenience.

    input_seq : str
        The FASTA sequence of the protein being built. This parameter
        needs to be accurate and is extremely important to create an
        appropriate data structure able to contain all conformer's
        atoms.
        Example: "MAGERDDAPL".

    conformer_name : str, optional
        The prefix name of the PDB file saved to disk.
        Defaults to 'conformer'.

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

    nconfs : int
        The number of conformers to build.
    """
    BUILD_BEND_H_N_C = build_bend_H_N_C
    CALC_DISTS = calc_all_vs_all_dists_square
    COUNT_NONZERO = np.count_nonzero
    DISTANCE_NH = distance_H_N
    LOGICAL_AND = np.logical_and
    MAKE_COORD_Q_COO_LOCAL = make_coord_Q_COO
    MAKE_COORD_Q_CO_LOCAL = make_coord_Q_CO
    MAKE_COORD_Q_LOCAL = make_coord_Q
    NAN = np.nan
    RC = randchoice
    RINT = randint
    ROUND = np.round
    angles = ANGLES
    slices = SLICES

    # tests generative function complies with implementation requirements
    if generative_function:
        try:
            generative_function(nres=1, cres=0)
        except Exception as err:  # this is generic Exception on purpose
            errmsg = (
                'The `generative_function` provided is not compatible with '
                'the building process. Please read `main_exec` docstring '
                'for more details.'
                )
            raise IDPConfGenException(errmsg) from err

    # Start building process

    # prepares data based on the input sequence
    # considers sidechain all-atoms
    atom_labels = np.array(generate_atom_labels(input_seq))
    num_atoms = len(atom_labels)
    residue_numbers = np.array(generate_residue_numbers(atom_labels))
    residue_labels = np.array(generate_residue_labels(input_seq))

    assert len(residue_numbers) == num_atoms
    assert len(residue_labels) == num_atoms, (len(residue_labels), num_atoms)

    vdW_sums, vdW_non_bond = generate_vdW_data(
        atom_labels,
        residue_numbers,
        residue_labels,
        vdW_radii_dict[vdW_radii],
        vdW_bonds_apart,
        vdW_tolerance,
        )

    # creates index-translated boolean masks
    bb_mask = np.where(np.isin(atom_labels, ('N', 'CA', 'C')))[0]
    carbonyl_mask = np.where(atom_labels == 'O')[0]
    NHydrogen_mask = np.where(atom_labels == 'H')[0]
    OXT_index = np.where(atom_labels == 'OXT')[0][0]

    # the last O value is removed because it is built in the same step
    # OXT is built
    OXT1_index = carbonyl_mask[-1]
    carbonyl_mask = carbonyl_mask[:-1]

    # create coordinates and views
    coords = np.full((num_atoms, 3), NAN, dtype=np.float64)

    # +1 because of the dummy coordinate required to start building.
    # see later
    bb = np.full((bb_mask.size + 1, 3), NAN, dtype=np.float64)
    bb_real = bb[1:, :]  # backbone coordinates without the dummy
    # coordinates for the carbonyl oxigen atoms
    bb_CO = np.full((carbonyl_mask.size, 3), NAN, dtype=np.float64)

    # notice that NHydrogen_mask does not see Prolines
    bb_NH = np.full((NHydrogen_mask.size, 3), NAN, dtype=np.float64)

    # the following array serves to avoid placing HN in Proline residues
    residue_labels_bb_simulating = residue_labels[bb_mask]

    # creates seed coordinates:
    # 1st) a dummy atom at the y-axis to build the first atom
    # 2nd) N-terminal N-atom is at 0, 0, 0
    # 3rd) CA atom of the firs residue is at the x-axis
    # # coordinates are created always from the parameters in the core
    # # definitions of IDPConfGen
    dummy_CA_m1_coord = np.array((0.0, 1.0, 0.0))
    n_terminal_N_coord = np.array((0.0, 0.0, 0.0))
    n_terminal_CA_coord = np.array((distances_N_CA[input_seq[0]], 0.0, 0.0))

    # seed coordinates array
    seed_coords = np.array((
        dummy_CA_m1_coord,
        n_terminal_N_coord,
        n_terminal_CA_coord,
        ))

    # side chain templates and masks
    # create a lits of tuples, where index 0 is a boolean mask placing
    # the side chain in coords, and index 1 is an array of the sidechain atoms
    # excluded the 4 atoms of the backbone

    # this is a list but could be a dictionary because key are indexes
    # the sublists must be lists and not tuples because they are modified
    # during the building process
    ss = [
        [
            np.where(residue_numbers == _resnum)[0][4:],
            np.full((_natoms - 4, 3), NAN, dtype=np.float64),
            ]
        for _resnum, _natoms in sorted(Counter(residue_numbers).items())
        ]

    # the last residue will contain an extra atom OXT, which needs to be
    # removed
    ss[-1][0] = ss[-1][0][:-1]
    ss[-1][1] = ss[-1][1][:-1]
    assert ss[-1][0].size == ss[-1][1].shape[0]

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

    broke_on_start_attempt = False
    start_attempts = 0
    max_start_attempts = 5  # maximum attempts to start a conformer
    # because we are building from a experimental database there can be
    # some angle combinations that fail on our validation process from start
    # if this happens more than `max_start_attemps` the production is canceled.

    # STARTS BUILDING
    # aligns the indexes so that conformers can be named properly in
    # multicore operations
    start_conf = nconfs * execution_run
    end_conf = start_conf + nconfs
    conf_n = start_conf
    # for conf_n in range(start_conf, end_conf):
    while conf_n < end_conf:
        print('building conformer: ', conf_n)

        # prepares cycles for building process
        # cycles need to be regenerated every conformer because the first
        # atom build is a C and the last atom built is the CA, which breaks
        # periodicity
        # all these are dictionaries
        bond_lens = cycle((
            distances_CA_C,  # used for PHI
            distances_C_Np1,  # used for PSI
            distances_N_CA,  # used for OMEGA
            ))
        bond_bend = cycle((
            build_bend_angles_N_CA_C,  # used for PHI
            build_bend_angles_CA_C_Np1,  # used for PSI
            build_bend_angles_Cm1_N_CA,  # used for OMEGA
            ))

        # in the first run of the loop this is unnecessary, but is better to
        # just do it once than flag it the whole time
        coords[:, :] = NAN
        bb[:, :] = NAN
        bb_CO[:, :] = NAN
        bb_NH[:, :] = NAN
        for _mask, _coords in ss:
            _coords[:, :] = NAN

        bb[:3, :] = seed_coords  # this contains a dummy coord at position 0

        # IMPLEMENT SIDECHAINS HERE

        bbi = 2  # starts at 2 because the first 3 atoms are already placed
        bbi0_R_CLEAR()
        bbi0_R_APPEND(bbi)

        COi = 0  # carbonyl atoms
        COi0_R_CLEAR()
        COi0_R_APPEND(COi)

        # NHi is 0 because it will be applied in bb_NH_builder which ignores
        # the first residue, see bb_NH_builder definition
        NHi = 1
        NHi0_R_CLEAR()
        NHi0_R_APPEND(NHi)

        # residue integer number
        current_res_number = 0
        res_R_CLEAR()
        res_R_APPEND(current_res_number)

        backbone_done = False
        number_of_trials = 0
        number_of_trials2 = 0
        # run this loop until a specific BREAK is triggered
        while True:

            # I decided to use an if-statement here instead of polymorph
            # the else clause to a `generative_function` variable because
            # the resulting overhead from the extra function call and
            # **kwargs handling was greater then the if-statement processing
            # https://pythonicthoughtssnippets.github.io/2020/10/21/PTS14-quick-in-if-vs-polymorphism.html
            if generative_function:
                agls = generative_function(nres=RINT(1, 6), cres=bbi - 2 // 3)

            else:
                # following `aligndb` function,
                # `angls` will always be cyclic with:
                # phi - psi - omega - phi - psi - omega - (...)
                agls = angles[RC(slices), :].ravel()

            # index at the start of the current cycle
            try:
                for torsion in agls:
                    # bbi -2 makes the match, because despite the first atom
                    # being placed is a C, it is placed using the PHI angle of a
                    # residue. And PHI is the first angle of that residue angle
                    # set.
                    current_res_number = (bbi - 2) // 3
                    current_residue = input_seq[current_res_number]

                    # bbi is the same for bb_real and bb, but bb_real indexing
                    # is displaced by 1 unit from bb. So 2 in bb_read is the
                    # next atom of 2 in bb.
                    bb_real[bbi, :] = MAKE_COORD_Q_LOCAL(
                        bb[bbi - 2, :],
                        bb[bbi - 1, :],
                        bb[bbi, :],
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

                # add the carboxyls
                coords[[OXT1_index, OXT_index]] = MAKE_COORD_Q_COO_LOCAL(
                    bb[-2, :],
                    bb[-1, :],
                    )

            # builds carbonyl atoms. Two situations can happen here:
            # 1) the backbone is not complete - the last atom is CA
            # 2) the backbone is complete - the last atom is C
            # whichever the case, with `bbi - 1` applying on the range bellow
            # will build carbonyls for the new chain expect for the last C
            # if the chain is completed.
            # this is so because range(X, Y, Z) equals range(X, Y-1, Z)
            # if Y-1 is not in range(X, Y, Z). And, this is the case for N-CA
            # pair.
            for k in range(bbi0_register[-1], bbi - 1, 3):
                bb_CO[COi, :] = MAKE_COORD_Q_CO_LOCAL(
                    bb_real[k - 1, :],
                    bb_real[k, :],
                    bb_real[k + 1, :],
                    )
                COi += 1

            # Adds N-H Hydrogens
            for k in range(bbi0_register[-1] + 1, bbi, 3):

                if residue_labels_bb_simulating[k] == 'PRO':
                    continue

                # MAKE_COORD_Q_CO_LOCAL can be used for NH by giving
                # disntace and bend parameters
                bb_NH[NHi, :] = MAKE_COORD_Q_CO_LOCAL(
                    bb_real[k - 1, :],
                    bb_real[k, :],
                    bb_real[k + 1, :],
                    distance=DISTANCE_NH,
                    bend=BUILD_BEND_H_N_C,
                    )
                NHi += 1

            # Adds sidechain template structures
            for res_i in range(res_R[-1], current_res_number + backbone_done):
                sscoords = place_sidechain_template(
                    bb_real[res_i * 3:res_i * 3 + 3, :],  # from N to C
                    sidechain_templates[aa1to3[input_seq[res_i]]],
                    )
                ss[res_i][1][:, :] = sscoords

            # Transfers coords to the main coord array
            for _smask, _sidecoords in ss[: current_res_number + backbone_done]:
                coords[_smask] = _sidecoords
            coords[bb_mask] = bb_real  # do not consider the initial dummy atom
            coords[carbonyl_mask] = bb_CO
            coords[NHydrogen_mask] = bb_NH

            # note that CALC_DISTS computes the square (without sqrt)
            distances = CALC_DISTS(coords)

            # compatibly, vdW_sums, considers squared distances
            clash = distances < vdW_sums

            # it is actually faster to compute everything and select only
            # those relevant
            valid_clash = LOGICAL_AND(clash, vdW_non_bond)

            # count number of True occurrences
            has_clash = COUNT_NONZERO(valid_clash)

            if has_clash:
                # reset coordinates to the original value
                # before the last chunk added

                # reset the same chunk maximum 5 times,
                # after that reset also the chunk before
                if number_of_trials > 5:
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
                    break  # while loop

                # clean previously built protein chunk
                bb_real[_bbi0:bbi, :] = NAN
                bb_CO[_COi0:COi, :] = NAN
                bb_NH[_NHi0:NHi, :] = NAN

                # do the same for sidechains
                # ...

                # reset also indexes
                bbi = _bbi0
                COi = _COi0
                NHi = _NHi0
                current_res_number = _resi0

                # coords needs to be reset because size of protein next
                # chunks may not be equal
                coords[:, :] = NAN

                # prepares cycles for building process
                # this is required because the last chunk created may have been
                # the final part of the conformer
                if backbone_done:
                    bond_lens = cycle((
                        distances_CA_C,
                        distances_C_Np1,
                        distances_N_CA,
                        ))
                    bond_bend = cycle((
                        build_bend_angles_N_CA_C,
                        build_bend_angles_CA_C_Np1,
                        build_bend_angles_Cm1_N_CA,
                        ))

                # we do not know if the next chunk will finish the protein
                # or not
                backbone_done = False
                number_of_trials += 1
                continue  # send back to the while loop

            # if the conformer is valid
            number_of_trials = 0
            bbi0_R_APPEND(bbi)
            COi0_R_APPEND(COi)
            NHi0_R_APPEND(NHi)
            # the residue where the build process stopped
            res_R_APPEND(current_res_number)

            if backbone_done:
                # this point guarantees all protein atoms are built
                break  # while loop
        # END of while loop

        if broke_on_start_attempt:
            start_attempts += 1
            if start_attempts > max_start_attempts:
                log.error(
                    'Reached maximum amount of starts. Canceling... '
                    f'{end_conf - conf_n} of {end_conf - start_conf} '
                    'conformers were not built.'
                    )
                return
            broke_on_start_attempt = False
            continue

        # until sidechains are implemented this is needed
        relevant = np.logical_not(np.isnan(coords[:, 0]))

        pdb_string = gen_PDB_from_conformer(
            input_seq,
            atom_labels[relevant],
            residue_numbers[relevant],
            ROUND(coords[relevant], decimals=3),
            )

        fname = f'{conformer_name}_{conf_n}.pdb'
        with open(fname, 'w') as fout:
            fout.write(pdb_string)

        conf_n += 1

    return


def gen_PDB_from_conformer(
        input_seq,
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
    ATOM_LABEL_FMT = ' {: <3}'.format

    for i in range(len(atom_labels)):

        if atom_labels[i] == 'N':
            resi += 1
            current_residue = input_seq[resi]
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
            AA1TO3[current_residue],
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


def generate_atom_labels(input_seq, AL=atom_labels):
    """."""
    labels = []
    LE = labels.extend

    for residue in input_seq:
        LE(AL[residue])

    labels.append('OXT')

    return labels


def generate_residue_numbers(atom_labels, start=1):
    """
    Create a list of residue numbers based on atom labels.

    Considers `N` to be the first atom of the residue.
    If this is not the case, the output can be meaningless.

    Returns
    -------
    list of ints
    """
    if atom_labels[0] == 'N':
        start -= 1  # to compensate the +=1 implementation in the for loop

    residues = []
    RA = residues.append

    for al in atom_labels:
        if al == 'N':
            start += 1
        RA(start)

    return residues


def generate_residue_labels(
        input_seq,
        atom_labels=atom_labels,
        aa1to3=aa1to3,
        ):
    """Generate residue 3-letter labels."""
    residue_labels = []
    RL_APPEND = residue_labels.append

    for res_1_letter in input_seq:
        res3 = aa1to3[res_1_letter]
        for _ in range(len(atom_labels[res_1_letter])):
            RL_APPEND(res3)

    # To compensate for the OXT atom
    residue_labels.append(res3)

    return residue_labels


def generate_vdW_data(
        atom_labels,
        residue_numbers,
        residue_labels,
        vdW_radii,
        bonds_apart=3,
        tolerance=0.4,
        ):
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
    assert len(atom_labels) == len(residue_numbers)
    assert len(atom_labels) == len(residue_labels)

    # we need to use re because hydrogen atoms start with integers some times
    atoms_char = re.compile(r'[CHNOPS]')
    findall = atoms_char.findall
    cov_topologies = generate_residue_template_topology()
    bond_structure_local = \
        expand_topology_bonds_apart(cov_topologies, bonds_apart)
    inter_connect_local = inter_residue_connectivities[bonds_apart]

    # adds OXT to the bonds connectivity, so it is included in the
    # final for loop creating masks
    add_OXT_to_residue(bond_structure_local[residue_labels[-1]])

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


if __name__ == "__main__":
    libcli.maincli(ap, main)
