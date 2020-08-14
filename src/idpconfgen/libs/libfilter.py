"""Contain functions to filter information from the DB."""
import re

import numpy as np

from idpconfgen.core.exceptions import IDPConfGenException


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


# # regex to compute
# forward with overlap
# forward no overlap


def regex_search(sequence, regex):
    """Search for regex in sequence)."""
    overlap_r = re.compile('\(\?\=\(.+\)')
    rex = re.compile(regex)

    if overlap_r.findall(regex):
        slices = regex_forward_with_overlap(sequence, rex)

    else:
        # do without overlap
        slices = regex_forward_no_overlap(sequence, rex)

    return slices



def regex_forward_no_overlap(sequence, regex):
    """
    Searches for regex forward without overlap.

    Examples:

        r'L'
        r'L{2}'
        r'L{1,3}'

    In the first example, returns all indexes of single char L without
        overlap.

    On the second example, returns all indexes of entire 'LL' sequences
        without overlap. So, 'LLL' returns only slice(0, 2).

    3) returns all sequences with 'LLL' without overlap, if a terminal 'LL'
        is found, returns that. Same for 'L' if found at the end.

    Using expressions such as r'(?=(L))' give not correct results.
    """
    # m.span() is used for regexes without overlap
    # using m.start(1) would not work here.
    # See regex_forward_with_overlap
    result = [slice(*m.span()) for m in regex.finditer(sequence)]
    if not result:
        raise IDPConfGenException(f'No matches found for: {regex.pattern}')
    else:
        return result


def regex_forward_with_overlap(sequence, regex):
    """
    Returns the indexes of re.finditer.

    Accepts only regex expressions with overlap, for example:

        r'(?=(L{3}))'
    """
    result = [slice(m.start(1), m.end(1)) for m in regex.finditer(sequence)]
    if not result:  # no match found
        # I decided to raise an error here because of the functionality of
        # IDPConfGen. It can lead to hidden bugs and errors if a regex without
        # matches is given, at least for now
        raise IDPConfGenException()#errmsg=f'No matches found for: {regex.pattern}')
    else:
        return result
