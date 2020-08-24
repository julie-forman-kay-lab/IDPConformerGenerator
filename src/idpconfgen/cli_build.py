"""
Builds.

USAGE:
    $ idpconfgen build DB

"""
import argparse
import sys
from functools import partial
from itertools import cycle
from multiprocessing import Pool, Manager
from random import choice
from time import time

import numpy as np

from idpconfgen import log
from idpconfgen.core.build_definitions import (
    atom_labels,
    build_bend_angles_CA_C_Np1,
    build_bend_angles_Cm1_N_CA,
    build_bend_angles_N_CA_C,
    distances_C_Np1,
    distances_CA_C,
    distances_N_CA,
    )
from idpconfgen.core.definitions import aa1to3
from idpconfgen.libs import libcli
from idpconfgen.libs.libcalc import (
    make_coord_Q,
    make_coord_Q_CO,
    make_coord_Q_COO,
    )
from idpconfgen.libs.libfilter import aligndb, regex_search
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libpdb import atom_line_formatter
from idpconfgen.libs.libtimer import timeme
from idpconfgen.libs.libvalidate import validate_conformer_for_builder


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
    help='Number of conformers to build',
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

    _slices, ANGLES = read_db_to_slices(database, dssp_regexes)
    SLICES.extend(_slices)

    # prepars execution function
    execute = partial(
        main_exec,
        input_seq=input_seq,  # string
        #slices=slices_multi,
        #database=database,  # path string
        #dssp_regexes=dssp_regexes,  # list of strings
        nconfs=core_chunks,  # int
        conformer_name=conformer_name,  # string
        )

    start = time()
    with Pool(ncores) as pool:
        imap = pool.imap(execute, range(ncores))
        for _ in imap:
            pass

    if remaining_chunks:
        execute(core_chunks * ncores, nconfs=remaining_chunks)

    log.info(f'{nconfs} conformers built in {time() - start:.3f} seconds')


def read_db_to_slices(database, dssp_regexes):
    # reads db dictionary from disk
    db = read_dictionary_from_disk(database)
    log.info(f'Read DB with {len(db)} entries')

    # reads and prepares IDPConfGen data base
    timed = partial(timeme, aligndb)
    pdbs, angles, dssp, resseq = timed(db)

    # searchs for slices in secondary structure, according to user requests
    timed = partial(timeme, regex_search)
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
        #slices,
        #database,
        conformer_name='conformer',
        #dssp_regexes=r'(?=(L{2,6}))',
        nconfs=1,
        ):
    """Prepare data base and builds conformers."""
    # bring global to local scope
    MAKE_COORD_Q_COO_LOCAL = make_coord_Q_COO
    MAKE_COORD_Q_CO_LOCAL = make_coord_Q_CO
    MAKE_COORD_Q_LOCAL = make_coord_Q
    NAN = np.nan
    RC = choice
    ROUND = np.round
    VALIDATE_CONF_LOCAL = validate_conformer_for_builder
    angles = ANGLES
    slices = SLICES

    # Start building process

    # prepares data based on the input sequence
    # considers sidechain all-atoms
    atom_labels = np.array(generate_atom_labels(input_seq))
    num_atoms = len(atom_labels)
    residue_numbers = np.array(generate_residue_numbers(atom_labels))

    # creates masks
    bb_mask = np.isin(atom_labels, ('N', 'CA', 'C'))
    carbonyl_mask = atom_labels == 'O'
    OXT_index = np.argwhere(atom_labels == 'OXT')[0][0]
    # replces the last TRUE value by false because that is a carboxyl
    # and not a carbonyl
    OXT1_index = np.argwhere(carbonyl_mask)[-1][0]
    carbonyl_mask[OXT1_index] = False

    # create coordinates and views
    coords = np.full((num_atoms, 3), NAN, dtype=np.float64)
    # +1 because of the dummy coordinate required to start building.
    # see later
    bb = np.full((np.sum(bb_mask) + 1, 3), NAN, dtype=np.float64)
    bb_real = bb[1:, :]  # backbone coordinates without the dummy
    # coordinates for the carbonyl oxigen atoms
    bb_CO = np.full((np.sum(carbonyl_mask), 3), NAN, dtype=np.float64)

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

    bbi0_register = []
    bbi0_R_APPEND = bbi0_register.append
    bbi0_R_POP = bbi0_register.pop
    bbi0_R_CLEAR = bbi0_register.clear

    COi0_register = []
    COi0_R_APPEND = COi0_register.append
    COi0_R_POP = COi0_register.pop
    COi0_R_CLEAR = COi0_register.clear

    # STARTS BUILDING
    # aligns the indexes so that conformers can be named properly in
    # multicore operations
    start_conf = nconfs * execution_run
    end_conf = start_conf + nconfs
    for conf_n in range(start_conf, end_conf):

        # prepares cycles for building process
        # cycles need to be regenerated every conformer because the first
        # atom build is a CA and the last atom built is the C, which breaks
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

        bb[:3, :] = seed_coords  # this contains a dummy coord at position 0

        # IMPLEMENT SIDECHAINS HERE

        bbi = 2  # starts at 2 because the first 3 atoms are already placed
        bbi0_R_CLEAR()
        bbi0_R_APPEND(bbi)

        COi = 0  # carbonyl atoms
        COi0_R_CLEAR()
        COi0_R_APPEND(COi)

        backbone_done = False
        number_of_trials = 0
        # run this loop until a specific BREAK is triggered
        while True:

            # following aligndb function, `angls` will always be cyclic with:
            # phi - psi - omega - phi - psi - omega - (...)
            agls = angles[RC(slices), :].ravel()

            # index at the start of the current cycle
            try:
                for torsion in agls:
                    # bbi -2 makes the match, because despite the first atom
                    # being placed is a C, it is placed using the PHI angle of a
                    # residue. And PHI is the first angle of that residue angle
                    # set.
                    current_residue = input_seq[(bbi - 2) // 3]

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
                coords[[OXT_index, OXT1_index]] = MAKE_COORD_Q_COO_LOCAL(
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

            # IMPLEMENT SIDE CHAIN CONSTRUCTION HERE

            # validate conformer current state
            coords[bb_mask] = bb_real  # do not consider the initial dummy atom
            coords[carbonyl_mask] = bb_CO

            energy = VALIDATE_CONF_LOCAL(
                coords,
                atom_labels,
                residue_numbers,
                bb_mask,
                carbonyl_mask,
                )

            if energy > 0:  # not valid
                # reset coordinates to the original value
                # before the last chunk added

                # reset the same chunk maximum 5 times,
                # after that reset also the chunk before
                if number_of_trials > 5:
                    bbi0_R_POP()
                    COi0_R_POP()
                    number_of_trials = 0

                try:
                    _bbi0 = bbi0_register[-1]
                    _COi0 = COi0_register[-1]
                except IndexError:
                    # if this point is reached,
                    # we erased until the beginning of the conformer
                    # discard conformer, something went really wrong
                    sys.exit('total error')  # change this to a functional error

                # clean previously built protein chunk
                bb_real[_bbi0:bbi, :] = NAN
                bb_CO[_COi0:COi, :] = NAN

                # do the same for sidechains
                # ...

                # reset also indexes
                bbi = _bbi0
                COi = _COi0

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

            if backbone_done:
                # this point guarantees all protein atoms are built
                break  # while loop
        # END of while loop

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


if __name__ == "__main__":
    libcli.maincli(ap, main)
