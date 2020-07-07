"""Contain functions to filter information from the DB."""
import re

import numpy as np


def aligndb(db):
    """Aligns IDPConfGen DB."""
    phi = []
    psi = []
    omg = []
    dssp = []
    resseq = []  # FASTA sequence
    pdbs = {}
    PSIE = psi.extend
    PHIE = phi.extend
    OMGE = omg.extend
    DA = dssp.append
    RS = resseq.append

    current = 0
    for pdb, data in db.items():

        len_segment = len(data['fasta'])
        _PHITMP = [np.nan] + data['phi'] + [np.nan]
        _PSITMP = data['psi'] + [np.nan, np.nan]
        _OMGTMP = data['omega'] + [np.nan, np.nan]

        l_PHITMP = len(_PHITMP)
        if sum([l_PHITMP, len(_PSITMP), len(_OMGTMP)]) / 3 != l_PHITMP:
            print(f'lengths of angles are not the same, ignoring... {pdb}')
            continue

        if len_segment + 1 != l_PHITMP:
            print(f"Seq length and angles don't match, ignoring... {pdb}")
            print(f'>>> {len_segment} {l_PHITMP} {pdb}')
            continue

        pdbs[pdb] = current
        # +1 because resseq will be concatenated with '|'
        # can't avoid +1 becase of the need to set the next starting integer
        current += len_segment + 1
        PHIE(_PHITMP)
        PSIE(_PSITMP)
        OMGE(_OMGTMP)

        DA(data['dssp'])
        RS(data['fasta'])

    _resseq = '|'.join(resseq)
    _dssp = '|'.join(dssp)
    _angles = np.array([phi, psi, omg], dtype=np.float32).T

    # this trick is used so that the last value of indexes is the
    # excluded index need to slice the aligned datastructure
    indexes = list(pdbs.values())
    values = indexes + [len(indexes)]

    # j - 1 is needed because j alone points to the '|'
    pdbs = {
        k: slice(i, j - 1)
        for k, i, j in zip(pdbs.keys(), values[:-1], values[1:])
        }

    return pdbs, _angles, _dssp, _resseq



def index_single_regex_overlap(sequence, regex):
    """
    Returns the indexes of re.finditer.

    """
    return [slice(m.start(1), m.end(1)) for m in regex.finditer(sequence)]
