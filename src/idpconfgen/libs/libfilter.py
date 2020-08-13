"""Contain functions to filter information from the DB."""
import re

import numpy as np


def aligndb(db, NAN=np.nan):
    """Aligns IDPConfGen DB."""
    phi, psi, omg, dssp, resseq = [], [], [], [], []
    pdbs = {}
    PSIE = psi.extend
    PHIE = phi.extend
    OMGE = omg.extend
    DA = dssp.append
    RA = resseq.append

    current = 0
    for pdb, data in db.items():

        len_segment = len(data['fasta'])
        _PHITMP = [NAN] + data['phi'] + [NAN]
        _PSITMP = data['psi'] + [NAN, NAN]
        _OMGTMP = data['omega'] + [NAN, NAN]

        l_PHITMP = len(_PHITMP)
        if sum([l_PHITMP, len(_PSITMP), len(_OMGTMP)]) / 3 != l_PHITMP:
            print(f'lengths of angles are not the same, ignoring... {pdb}')
            continue

        # +1 because NAN are added as spacers
        len_segment_w_spacer = len_segment + 1
        if len_segment_w_spacer != l_PHITMP:
            print(f"Seq length and angles don't match, ignoring... {pdb}")
            print(f'>>> {len_segment} {l_PHITMP} {pdb}')
            continue

        pdbs[pdb] = slice(current, current + len_segment)
        # +1 because resseq will be concatenated with '|'
        # can't avoid +1 because of the need to set the next starting integer
        current += len_segment_w_spacer

        PHIE(_PHITMP)
        PSIE(_PSITMP)
        OMGE(_OMGTMP)

        DA(data['dssp'])
        RA(data['fasta'])

    _resseq = '|'.join(resseq)
    _dssp = '|'.join(dssp)
    _angles = np.array((phi, psi, omg), dtype=np.float32).T

    return pdbs, _angles, _dssp, _resseq



def index_single_regex_overlap(sequence, regex, start_offset=1):
    """
    Returns the indexes of re.finditer.

    Offset accounts for the difference betwee the explored sequence
    and the target data where to use the produced slice.

    For example, to use on torsion angles table,
    the stop slice found in the 
    """
    return [
        slice(m.start(1) + start_offset, m.end(1))
        for m in regex.finditer(sequence)
        ]
