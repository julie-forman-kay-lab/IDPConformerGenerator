"""Contain functions to filter information from the DB."""
import itertools as it
import re
from functools import partial

import numpy as np

from idpconfgen import log
from idpconfgen.libs.libmulticore import pool_function


REGEX_OVERLAP = re.compile(r'\(\?\=\(.+\)')
REGEX_RANGE = re.compile(r'(\{\d\,\d\}|\{\d\}|\{\d\,\}|\{\,\d\})')
REGEX_RANGE_CHAR = re.compile(r'\w\{')


def aligndb(db, NAN=np.nan):
    """Aligns IDPConfGen DB."""
    phi, psi, omg, dssp, resseq = [], [], [], [], []
    pdbs = {}
    PSIE = psi.extend
    PHIE = phi.extend
    OMGE = omg.extend
    DA = dssp.append
    RA = resseq.append
    spacer = 1  # algorithm definition

    current = 0
    for pdb, data in db.items():

        len_segment = len(data['fasta'])
        _PHITMP = [NAN] + data['phi'] + [NAN]
        _PSITMP = data['psi'] + [NAN, NAN]
        _OMGTMP = data['omega'] + [NAN, NAN]

        # +1 because NAN are added as spacers
        len_segment_w_spacer = len_segment + spacer

        correct_lengths = [
            len_segment_w_spacer,
            len(data['dssp']) + spacer,
            len(_PHITMP),
            len(_PSITMP),
            len(_OMGTMP),
            ]
        if sum(correct_lengths) / len(correct_lengths) != len_segment_w_spacer:
            log.error(
                'number of residues, SS chars and angles does not match, '
                f'ignoring... {pdb}'
                )
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


def regex_search(sequence, regex_string, rex_range=REGEX_RANGE, **kwargs):
    """
    Search for regex in sequence.

    Return
    ------
    list of slices
    """
    # if a range exists in regex_string
    # range is defined by default by: L{1}, L{1,5} situations
    # the following functions ensures searchs goes both directions
    if rex_range.findall(regex_string):
        result = regex_range(sequence, regex_string, **kwargs)

    else:
        overlap = regex_has_overlap(regex_string)
        func = regex_forward_with_overlap \
            if overlap \
            else regex_forward_no_overlap
        result = func(sequence, regex_string)

    assert isinstance(result, list)
    assert all(isinstance(S, slice) for S in result)  # heavy and slow!!
    return result


def regex_range(sequence, regex_string, ncores=1):
    """
    Find slices of sequence where regex_string applies.

    Searches forward and backwards by means of exploring ranges
    as individual regexes.
    """
    # prefix and suffix regex for overlap
    if regex_has_overlap(regex_string):
        pre, suf, func = r'(?=(', r'))', regex_forward_with_overlap
    else:
        pre, suf, func = '', '', regex_forward_no_overlap

    # creates all possible regex combinations without ranges
    # from the rangex identified in the original regex_string
    regex_combinations = make_regex_combinations_from_ranges(
        regex_string,
        pre=pre,
        suf=suf,
        )

    exec_pool = pool_function(
        partial(func, sequence),
        regex_combinations,
        ncores=ncores,
        )

    slices = []
    for result in exec_pool:
        slices.extend(result)

    return slices


def regex_has_overlap(regex_string, overlap_rex=REGEX_OVERLAP):
    """
    Find if a `regex_string` defines overlap.

    Parameters
    ----------
    regex_string : str
        The regex string.

    overlap_fmt : str
        The regex to find overlap in regex_string.

    Returns
    -------
    bool
    """
    return bool(overlap_rex.findall(regex_string))


def make_regex_combinations_from_ranges(regex_string, **kwargs):
    """."""
    ranges, chars = make_ranges(regex_string)
    return make_regex_combinations(ranges, chars, **kwargs)


def make_ranges(
        regex_string,
        rang_rex=REGEX_RANGE,
        char_rex=REGEX_RANGE_CHAR,
        max_range=30,
        ):
    """
    Define a set of ranges and characters from `regex_string`.

    Examples
    --------
        >>> make_range('L{1}H{1,2}')
        [range(1, 2), range(1, 3)], ['L', 'H']

        >>> make_range('L{1,}H{,2}')
        [range(1, None), range(1, 3)], ['L', 'H']
    """
    # requires
    assert isinstance(regex_string, str)

    # examples of i according to `prev`
    # 'L{', 'H{'
    chars = [i[:-1] for i in char_rex.findall(regex_string)]

    # yes, I tried to write this in a list comprehension, but these are
    # 3 for loops. I thought is just more readable to make it explicit
    rangs = (i.strip('{}') for i in rang_rex.findall(regex_string))
    ranges = []
    for trange in rangs:
        # examples of trange:
        # '1', '1,', '1,2', ',2'
        # to make ranges compatible with regex, start should be equal to 1
        # and because regex ranges are inclusive, 1 must be added to stop
        ts = trange.split(',')

        # 'or 1' covers {,3} cases, that yield ['', '3']
        start = int(ts[0] or 1)

        try:
            # 'or max_range' covers {3,} cases, that yield ['3', '']
            end = int(ts[1] or max_range)
        except IndexError:
            # covers {1} no-range situations
            end = start

        ranges.append(range(start, end + 1))

    # ensures
    assert isinstance(ranges, list)
    assert isinstance(chars, list)
    assert len(ranges) == len(chars)
    return ranges, chars


def make_regex_combinations(ranges, chars, pre=None, suf=None):
    """
    Make combinations of regexes from `ranges` and `chars`.

    This function is not a general abstraction. Is instead a partial
    abstraction within the problem of IDPConfGen.

    Parameters
    ----------
    ranges : list of range objects
        Ranges where to start and stop searching in the regex.
        Ranges should follow regex conventions, usually 1-indexed
        and stop inclusive.

    chars : list of 1 letter chars
        The chars that will be searched in ranges.

    pre, suf : str
        Strings to add as prefix and suffix of the generates ranges.
    """
    pre = pre or ''
    suf = suf or ''
    regex_combinations = []
    for range_tuple in it.product(*ranges):
        c_regex = (c + '{' + str(ii) + '}' for ii, c in zip(range_tuple, chars))
        regex_combinations.append(f"{pre}{''.join(c_regex)}{suf}")
    return regex_combinations


def regex_forward_no_overlap(sequence, regex):
    r"""
    Search for regex forward without overlap.

    Examples
    --------

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
    regex_c = re.compile(regex)
    return [slice(*m.span()) for m in regex_c.finditer(sequence)]


def regex_forward_with_overlap(sequence, regex):
    """
    Find matches for regex in sequence of chars.

    Considers regex defines overlap.

    Accepts only regex expressions with overlap, for example:

        r'(?=(L{3}))'

    Returns
    -------
    list of slices
        Where slice.start and slice.stop follow Python conventions.
    """
    regex_c = re.compile(regex)
    return [slice(m.start(1), m.end(1)) for m in regex_c.finditer(sequence)]
