"""
Builds.

USAGE:
    $ idpconfgen build DB

"""
import argparse
import re
import sys
from functools import partial
from itertools import cycle
from multiprocessing import Pool
from random import choice as RC

import numpy as np

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libcalc import make_coord_Q, make_coord_Q_CO, make_coord_Q_COO
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libfilter import (
    aligndb,
    regex_search,
    )
from idpconfgen.libs.libtimer import timeme
from idpconfgen.core.build_definitions import (
    build_bend_angles_CA_C_Np1,
    build_bend_angles_Cm1_N_CA,
    build_bend_angles_N_CA_C,
    distances_N_CA,
    distances_CA_C,
    distances_C_Np1,
    atom_labels,
    )

from idpconfgen.core.definitions import aa1to3
from idpconfgen.libs.libpdb import atom_line_formatter
from idpconfgen.libs.libvalidate import validate_conformer_for_builder
from idpconfgen.libs.libmulticore import pool_function


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


def main(
        input_seq,
        database,
        func=None,
        dssp_regexes=r'(?=(L{2,6}))',
        nconfs=1,
        conformer_name='conformer',
        ncores=1,
        ):
    """."""
    core_chunks = nconfs // ncores
    remaining_chunks = nconfs % ncores

    #execute = partial(
    main_exec(
        0,
        input_seq=input_seq,  # string
        database=database,  # path string
        dssp_regexes=dssp_regexes,  # list of strings
        nconfs=core_chunks,  # int
        conformer_name=conformer_name, # string
        )

    #from time import time
    #start = time()
    #with Pool(ncores) as pool:
    #    imap = pool.imap(execute, range(ncores))
    #    for _ in imap:
    #        pass

    #if remaining_chunks:
    #    execute(core_chunks * ncores, nconfs=remaining_chunks)
    #print(time() - start)


def main_exec(
        execution_run,
        input_seq,
        database,
        dssp_regexes=r'(?=(L{2,6}))',
        nconfs=1,
        conformer_name='conformer',
        ROUND=np.round,
        ):
    """."""
    # bring global to local scope
    MAKE_COORD_Q_LOCAL = make_coord_Q
    MAKE_COORD_Q_CO_LOCAL = make_coord_Q_CO
    MAKE_COORD_Q_COO_LOCAL = make_coord_Q_COO
    VALIDATE_CONF_LOCAL = validate_conformer_for_builder
    NAN = np.nan


    db = read_dictionary_from_disk(database)
    ldb = len(db)
    log.info(f'Read DB with {ldb} entries')

    # reads and aligns IDPConfGen data base
    timed = partial(timeme, aligndb)
    pdbs, angles, dssp, resseq = timed(db)

    # seachs for slices in secondary structure
    #timed = partial(timeme, regex_search)
    #slices = []
    #if isinstance(dssp_regexes, str):
    #    dssp_regexes = [dssp_regexes]

    #for dssp_regex_string in dssp_regexes:
    #    slices.extend(timed(dssp, dssp_regex_string))
    slices = [pdbs[list(pdbs.keys())[0]]]

    log.info(f'Found {len(slices)} indexes for {dssp_regexes}')


    # building

    # prepares data based on the input sequence
    #len_conf = len(input_seq)  # number of residues
    atom_labels = np.array(generate_atom_labels(input_seq))  # considers sidechain all-atoms
    num_atoms = len(atom_labels)
    residue_numbers = np.array(generate_residue_numbers(atom_labels))

    # creates masks
    bb_mask = np.isin(atom_labels, ('N', 'CA', 'C'))
    carbonyl_mask = np.isin(atom_labels, ('O',))
    OXT_index = np.argwhere(atom_labels == 'OXT')[0][0]
    # replces the last TRUE value by false because that is a carboxyl
    # and nota carbonyl
    OXT1_index = np.argwhere(carbonyl_mask)[-1][0]
    carbonyl_mask[OXT1_index] = False

    # create coordinates and views
    coords = np.full((num_atoms, 3), NAN, dtype=np.float32)
    bb = np.full((np.sum(bb_mask) + 1, 3), NAN, dtype=np.float64)
    bb_real = bb[1:, :]  # without the dummy
    bb_CO = np.full((np.sum(carbonyl_mask), 3), NAN, dtype=np.float64)

    # places seed coordinates
    # coordinates are created always from the parameters in the core
    # definitions of IDPConfGen
    # first atom (N-terminal) is at 0, 0, 0
    # second atom (CA of the firs residue) is at the x-axis
    dummy_CA_m1_coord = np.array((0.0, 1.0, 0.0))
    n_terminal_N_coord = np.array((0.0, 0.0, 0.0))
    n_terminal_CA_coord = np.array((distances_N_CA[input_seq[0]], 0.0, 0.0))

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
    start_conf = nconfs * execution_run
    end_conf = start_conf + nconfs
    for conf_n in range(start_conf, end_conf):

        # prepares cycles for building process
        bond_lens = cycle((distances_CA_C, distances_C_Np1, distances_N_CA))
        bond_bend = cycle((
            build_bend_angles_N_CA_C,
            build_bend_angles_CA_C_Np1,
            build_bend_angles_Cm1_N_CA,
            ))

        # in the first run of the loop this is unnecessary, but is better of
        # just do it once than flag it the whole time
        coords[:, :] = NAN
        bb[:, :] = NAN
        bb_CO[:, :] = NAN

        bb[:3, :] = seed_coords  # this contains a dummy coord at position 0
        # SIDECHAINS HERE

        bbi = 2  # starts at 2 because the first 3 atoms are already placed
        bbi0_R_CLEAR()
        bbi0_R_APPEND(bbi)

        # and needs to adjust with the += assignment inside the loop
        COi = 0  # carbonyl atoms
        COi0_R_CLEAR()
        COi0_R_APPEND(COi)

        backbone_done = False
        number_of_trials = 0
        # run this loop until a specific BREAK is triggered
        while True:

            # the slice [1:-2] removes the first phi and the last psi and omega
            # from the group of angles. These angles are not needed because the
            # implementation always follows the order: psi-omega-phi(...)
            agls = angles[RC(slices), :].ravel()#[1:-2]

            # index at the start of the current cycle
            try:
                for torsion in agls:
                    # bbi -2 makes the match
                    current_residue = input_seq[(bbi - 2) // 3]
                    bbend_a = next(bond_lens)[current_residue]
                    print(current_residue, bbend_a)

                    bb_real[bbi, :] = MAKE_COORD_Q_LOCAL(
                        bb[bbi - 2, :],
                        bb[bbi - 1, :],
                        bb[bbi, :],
                        bbend_a,
                        next(bond_bend)[current_residue],
                        torsion,
                        )
                    bbi += 1
            except IndexError:
                # here bbi is the last index + 1
                backbone_done = True

                for k in range(bbi0_register[-1], bbi - 1, 3):
                    bb_CO[COi, :] = MAKE_COORD_Q_CO_LOCAL(
                        bb_real[k - 1, :],
                        bb_real[k, :],
                        bb_real[k + 1, :],
                        )
                    COi += 1

                coords[[OXT_index, OXT1_index]] = MAKE_COORD_Q_COO_LOCAL(
                    bb[-2, :],
                    bb[-1, :],
                    )

            # this else performs if the for loop concludes properly
            #else:
            # when backbone completes,
            # adds carbonyl oxygen atoms
            # after this loop all COs are added for the portion of BB
            # added previously
            else:
                for k in range(bbi0_register[-1], bbi, 3):
                    bb_CO[COi, :] = MAKE_COORD_Q_CO_LOCAL(
                        bb_real[k - 1, :],
                        bb_real[k, :],
                        bb_real[k + 1, :],
                        )
                    COi += 1

            #if backbone_done:  # make terminal carboxyl coordinates
            #    coords[[OXT_index, OXT1_index]] = MAKE_COORD_Q_COO_LOCAL(
            #        bb[-2, :],
            #        bb[-1, :],
            #        )

            # add sidechains here.....

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

            if False:#energy > 0:  # not valid
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
                bb_real[_bbi0:bbi, :] = np.nan
                bb_CO[_COi0:COi, :] = np.nan

                # do the same for sidechains
                # ...

                # reset also indexes
                bbi = _bbi0
                COi = _COi0

                # coords needs to be reset because size of protein next
                # chunks may not be equal
                coords[:, :] = np.nan

                # review!!!!!!!!
                # prepares cycles for building process
                bond_lens = cycle((distances_CA_C, distances_C_Np1, distances_N_CA))
                bond_bend = cycle((
                    build_bend_angles_N_CA_C,
                    build_bend_angles_CA_C_Np1,
                    build_bend_angles_Cm1_N_CA,
                    ))

                backbone_done = False
                number_of_trials += 1
                continue

            number_of_trials = 0
            bbi0_R_APPEND(bbi)
            COi0_R_APPEND(COi)

            if backbone_done:
                # this point guarantees all protein atoms are built
                break
        # END of while loop

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
            #print(f'saved: {fname}')

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
    pdb_text = ''

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
