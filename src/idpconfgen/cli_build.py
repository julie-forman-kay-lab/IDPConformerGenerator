"""
Builds.

USAGE:
    $ idpconfgen build DB

"""
import argparse
import re
import random
from functools import partial
from itertools import cycle

import numpy as np

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libcalc import make_coord_Q
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libfilter import (
    aligndb,
    index_single_regex_overlap,
    )
from idpconfgen.libs.libtimer import timeme
from idpconfgen.core.definitions import (
    distance_N_CA,
    distance_CA_C,
    distance_C_Np1,
    average_N_CA_C,
    average_CA_C_Np1,
    average_Cm1_N_CA,
    average_CA_C_O,
    )
from idpconfgen.core.exceptions import IDPConfGenException


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
    'database',
    help='The IDPConfGen database.',
    )

ap.add_argument(
    'input_seq',
    )


def main(input_seq, database, func=None):
    """."""

    db = read_dictionary_from_disk(database)
    ldb = len(db)
    _log = f'Read DB with {ldb} entries'
    log.info(_log)

    timed = partial(timeme, aligndb)
    pdbs, angles, dssp, resseq = timed(db)

    # esta parte tengo que ponerla de parametro externo
    regex = re.compile(r'(?=(L{6}))')
    timed = partial(timeme, index_single_regex_overlap)
    loops_6 = timed(dssp, regex)
    _log = f'Found {len(loops_6)} indexes for loops of 6'
    log.info(_log)

    # building

    # prepares data based on the input sequence
    len_conf = len(input_seq)  # number of residues
    atom_labels = generate_atom_labels(input_seq)  # considers sidechain all-atoms
    num_atoms = len(atom_labels)
    residue_numbers = generate_residue_numbers(atom_labels)
    coords = np.zeros((num_atoms, 3), dtype=np.float32)

    # creates masks
    bb_mask = np.isin(atom_labels, ('N', 'CA', 'C'))
    carbonyl_mask = np.isin(atom_labels, ('O',))

    # creates views
    bb = np.zeros((np.sum(bb_mask), 3), dtype=np.float32)
    bb_o = np.zeros((np.sum(carbonyl_mask), 3), dtype=np.float32)

    # places seed coordinates
    # coordinates are created always from the parameters in the core
    # definitions of IDPConfGen

    # first atom (N-terminal) is at 0, 0, 0 - no need to add it because
    # coords was created with np.zeros :-)

    # second atom (CA of the firs residue) is at the x-axis
    bb[1, :] = (distance_N_CA, 0.0, 0.0)

    # third atom (C of the first residue) needs to be computed according
    # to the bond length parameters and bend angle.
    bb[2, :] = make_coord_Q(
        np.array((0.0, distance_N_CA, 0.0)),  # dummy coordinate used only here
        bb[0, :],
        bb[1, :],
        distance_CA_C,
        average_N_CA_C,
        0,  # null torsion angle
        )

    # prepares cycles for building process
    bond_lens = cycle([distance_C_Np1, distance_N_CA, distance_CA_C])
    bond_bend = cycle([average_CA_C_Np1, average_Cm1_N_CA, average_N_CA_C])

    RC = random.choice
    ALL = np.all
    i = 3  # starts at 2 because the first 3 atoms are already placed
    # and needs to adjust with the += assignment inside the loop
    last_O = 0  # carbonyl atoms
    while True:
        agls = angles[RC(loops_6), :].ravel()[1:-2]
        i0 = i
        try:
            for torsion in agls:
                bb[i, :] = make_coord_Q(
                    bb[i - 3, :],
                    bb[i - 2, :],
                    bb[i - 1, :],
                    next(bond_lens),
                    next(bond_bend),
                    torsion,
                    )
                i += 1

            else:
                # when backbone completes,
                # adds carbonyl oxygen atoms
                for k in range(i0, i, 3):
                    bb_o[last_O, :] = make_coord_Q(
                        bb[k - 2, :],
                        bb[k - 1, :],
                        bb[k, :],
                        distance_N_CA,  # to change
                        average_CA_C_O,
                        0,
                        )
                    last_O += 1
        except IndexError:
            # i-1, to discard the last carboxyl
            for k in range(i0, i-1, 3):
                bb_o[last_O, :] = make_coord_Q(
                    bb[k - 2, :],
                    bb[k - 1, :],
                    bb[k, :],
                    distance_N_CA,  # to change
                    average_CA_C_O,
                    0,
                    )
                last_O += 1
            break


# TODO:
# clashes
# sidechains
# carboxylo final
# correct distance C_O
# correct bend angle CA_C_O

    coords[bb_mask] = bb
    coords[carbonyl_mask] = bb_o
    print(coords)

def generate_atom_labels(input_seq):
    return ('N', 'CA', 'C', 'O') * len(input_seq)


def generate_residue_numbers(atom_labels):
    return [
        number
        for number, atom_label in enumerate(atom_labels, start=1)
        if atom_label == 'N'
        ]


#class Conformer:
#    def __init__(self, input_seq):
#        len_conf = len(input_seq)
#        num_atoms = get_num_atoms(input_seq)
#
#        atom_labels = get_atom_labels(input_seq)
#        residue_numbers = get_residue_numbers(input_seq)
#
#        bb_mask = np.isin(atom_labels, ('N', 'CA', 'C'))
#        bb_o_mask = np.isin(atom_labels, ('O',))
#
#        coords = np.zeros((num_atoms, 3), dtype=np.float32)
#        bb = coords[bb_mask]
#        bb_o = coords[bb_o_mask]
#        bb[1, :] = (distance_N_CA, 0.0, 0.0)
#        bb[2, :] = make_coord_Q(
#            np.array((0.0, distance_N_CA, 0.0)),
#            bb[0, :],
#            bb[1, :],
#            distance_CA_C,
#            average_N_CA_C,
#            0,
#            )
#
#    def is_complete(self):
#        return np.all(bb)





    #from time import time
    #start = time()
    #for i in range(500_000):
    #    agl = random.choice(loops_6)
    #    agls = angles[agl, :]
    #print(time() - start)

    #return


def build(input_seq):
    return



if __name__ == "__main__":
    libcli.maincli(ap, main)
