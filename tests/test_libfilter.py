"""Test libfilter functions."""
import re

import numpy as np
import pytest

from idpconfgen.libs.libfilter import (
    aligndb,
    make_ranges,
    make_regex_combinations,
    regex_forward_no_overlap,
    regex_forward_with_overlap,
    regex_has_overlap,
    regex_range,
    regex_search,
    )


@pytest.mark.parametrize(
    'in1,out',
    [
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5",
                    },
                },
            (
                {"pdb1": slice(0, 6, None)},
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL",
                "MMMMMM",
                ),
            ),
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "HHHHHH",
                    "fasta": "AAAAAA",
                    "omega": [5, 5, 5, 5, 5],
                    "psi": [4, 4, 4, 4, 4],
                    "phi": [6, 6, 6, 6, 6],
                    "resids": "7,8,9,10,11,12",
                    },
                },
            (
                {
                    "pdb1": slice(0, 6, None),
                    "pdb2": slice(7, 13, None),
                    },
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    [np.nan, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL|HHHHHH",
                "MMMMMM|AAAAAA",
                ),
            ),
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "HHHHHH",
                    "fasta": "AAAAAA",
                    "omega": [5, 5, 5, 5, 5],
                    "psi": [4, 4, 4, 4, 4],
                    "phi": [6, 6, 6, 6, 6],
                    "resids": "7,8,9,10,11,12",
                    },
                "pdb3": {
                    "dssp": "EEEEEE",
                    "fasta": "WWWWWW",
                    "omega": [8, 8, 8, 8, 8],
                    "psi": [7, 7, 7, 7, 7],
                    "phi": [9, 9, 9, 9, 9],
                    "resids": "13,14,15,16,17,18",
                    },
                },
            (
                {
                    "pdb1": slice(0, 6, None),
                    "pdb2": slice(7, 13, None),
                    "pdb3": slice(14, 20, None),
                    },
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    [np.nan, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    [np.nan, 7, 8],
                    [9, 7, 8],
                    [9, 7, 8],
                    [9, 7, 8],
                    [9, 7, 8],
                    [9, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL|HHHHHH|EEEEEE",
                "MMMMMM|AAAAAA|WWWWWW",
                ),
            ),
        ],
    )
def test_aligndb(in1, out):
    """."""
    assert isinstance(in1, dict)
    pdbs, angles, dssp, fasta = aligndb(in1)
    assert pdbs == out[0]
    assert angles.shape == out[1].shape
    assert np.array_equal(angles, out[1], equal_nan=True)
    assert dssp == out[2]
    assert fasta == out[3]


@pytest.mark.parametrize(
    'in1, out',
    [
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "HHHHHH",
                    "fasta": "AAAAA",
                    "omega": [5, 5, 5, 5, 5],
                    "psi": [4, 4, 4, 4, 4],
                    "phi": [6, 6, 6, 6, 6],
                    "resids": "7,8,9,10,11,12",
                    },
                },
            (
                {
                    "pdb1": slice(0, 6, None),
                    },
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL",
                "MMMMMM",
                ),
            ),
        ]
    )
def test_aligndb_continue_1(in1, out):
    """Test ignores when number of angles doesn't match seq length."""
    assert isinstance(in1, dict)
    pdbs, angles, dssp, fasta = aligndb(in1)
    assert pdbs == out[0]
    assert angles.shape == out[1].shape
    assert np.array_equal(angles, out[1], equal_nan=True)
    assert dssp == out[2]
    assert fasta == out[3]


@pytest.mark.parametrize(
    'in1, out',
    [
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "HHHHH",
                    "fasta": "AAAAAA",
                    "omega": [5, 5, 5, 5, 5],
                    "psi": [4, 4, 4, 4, 4],
                    "phi": [6, 6, 6, 6, 6],
                    "resids": "7,8,9,10,11,12",
                    },
                },
            (
                {
                    "pdb1": slice(0, 6, None),
                    },
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL",
                "MMMMMM",
                ),
            ),
        ]
    )
def test_aligndb_continue_2(in1, out):
    """
    Test ignores when sizes differ.

    In other words, number of angles, residues or DSSP are different.
    """
    assert isinstance(in1, dict)
    pdbs, angles, dssp, fasta = aligndb(in1)
    assert pdbs == out[0]
    assert angles.shape == out[1].shape
    assert np.array_equal(angles, out[1], equal_nan=True)
    assert dssp == out[2]
    assert fasta == out[3]


@pytest.mark.parametrize(
    'seq, regex, expected_slices',
    [
        # ################## FIST Example
        # L between 1 and 3 length.
        (
            "LLLLLLLLHHLL",
            r"(?=(L{1,3}))",
            [
                slice(0, 3, None),
                slice(1, 4, None),
                slice(2, 5, None),
                slice(3, 6, None),
                slice(4, 7, None),
                slice(5, 8, None),
                slice(6, 8, None),
                slice(7, 8, None),
                slice(10, 12, None),
                slice(11, 12, None),
                ],
            ),
        # ################## SECOND Example
        # L of length 1
        (
            "HHHLLHHH",
            r"(?=(L))",
            [
                slice(3, 4),
                slice(4, 5),
                ],
            ),
        ],
    )
def test_regex_forward_with_overlap(seq, regex, expected_slices):
    """."""
    result = regex_forward_with_overlap(seq, re.compile(regex))

    assert len(result) == len(expected_slices)
    # the line bellow does not report on which slice is wrong
    # assert all(r == e for r, e in zip(result, expected_slices))
    # verbose for loop does report properly :-)
    for r, e in zip(result, expected_slices):
        assert r == e, (r, e)


@pytest.mark.parametrize(
    'seq, regex, expected_slices',
    [
        # ################## FIST Example
        # L between 1 and 3 length.
        (
            "LLLLLLLLHHLL",
            r"L{3}",
            [
                slice(0, 3, None),
                slice(3, 6, None),
                ],
            ),
        # ################## SECOND Example
        # L of length 1
        (
            "HHHLLHHH",
            r"H",
            [
                slice(0, 1),
                slice(1, 2),
                slice(2, 3),
                slice(5, 6),
                slice(6, 7),
                slice(7, 8),
                ],
            ),
        # ################## THIRD Example
        # L of length 1
        (
            "HHHLLHHH",
            r"H{2}",
            [
                slice(0, 2),
                slice(5, 7),
                ],
            ),
        # ################## FORTH Example
        # L of length 1
        (
            "HHHLLHHH",
            r"H{2}L{2}H{1,2}",
            [
                slice(1, 7),
                ],
            ),
        ],
    )
def test_regex_forward_no_overlap(seq, regex, expected_slices):
    """."""
    result = regex_forward_no_overlap(seq, re.compile(regex))

    assert len(result) == len(expected_slices)
    # the line bellow does not report on which slice is wrong
    # assert all(r == e for r, e in zip(result, expected_slices))
    # verbose for loop does report properly :-)
    for r, e in zip(result, expected_slices):
        assert r == e, (r, e)


@pytest.mark.parametrize(
    'in1,max_range,expected_ranges,expected_chars',
    [
        (
            'L{1}H{2}',
            20,
            [range(1, 2), range(2, 3)],
            ['L', 'H'],
            ),
        (
            'L{1,}L{8,9}H{5,}',
            20,
            [range(1, 21), range(8, 10), range(5, 21)],
            ['L', 'L', 'H'],
            ),
        (
            'L{,3}',
            10,
            [range(1, 4)],
            ['L'],
            ),
        ]
    )
def test_make_range(in1, max_range, expected_ranges, expected_chars):
    """."""
    ranges, chars = make_ranges(in1, max_range=max_range)

    assert len(ranges) == len(expected_ranges)

    for i, j in zip(ranges, expected_ranges):
        assert i == j

    assert len(chars) == len(expected_chars)

    for i, j in zip(chars, expected_chars):
        assert i == j


@pytest.mark.parametrize(
    'in_ranges, in_chars, in_pre, in_suf, expected',
    [
        (
            [range(4, 8)],
            'X',
            None,
            None,
            ['X{4}', 'X{5}', 'X{6}', 'X{7}'],
            ),
        (
            [range(1, 2), range(5, 9)],
            ['L', 'E'],
            r'(?=(',
            r'))',
            [
                '(?=(L{1}E{5}))',
                '(?=(L{1}E{6}))',
                '(?=(L{1}E{7}))',
                '(?=(L{1}E{8}))',
                ],
            ),
        ]
    )
def test_make_regex_combinations(
        in_ranges,
        in_chars,
        in_pre,
        in_suf,
        expected,
        ):
    """."""
    result = make_regex_combinations(in_ranges, in_chars, in_pre, in_suf)
    assert len(result) == len(expected)
    for i, j in zip(result, expected):
        assert i == j, (i, j)


@pytest.mark.parametrize(
    'in1, expected',
    [
        ('L{4,5}', False),
        ('(?=(L{4,5}))', True),
        ],
    )
def test_regex_has_overlap(in1, expected):
    """Test has overlap."""
    result = regex_has_overlap(in1)
    assert result == expected


# this list will be used in test_regex_search and test_regex_range
# ORDER MATHERS!
REGEX_SEARCH_w_RANGE = [
    (
        'LLLHHLLLE|HLLLEELLL',
        'L{3}',
        [
            slice(0, 3),
            slice(5, 8),
            slice(11, 14),
            slice(16, 19),
            ],
        ),
    (
        'LLLHHLLLE|HLLLEELLL',
        'L{1}H{2}',
        [slice(2, 5)],
        ),
    (
        'LLLHHLLLE|HLLLEELLL',
        'L{1,3}H{1,2}L{1,3}',
        [
            slice(2, 6),
            slice(2, 7),
            slice(2, 8),
            slice(1, 6),
            slice(1, 7),
            slice(1, 8),
            slice(0, 6),
            slice(0, 7),
            slice(0, 8),
            ],
        ),
    (
        'LLLHHLLLE|HLLLEELLL',
        'L{1}E{2}L{2,3}',
        [
            slice(13, 18),
            slice(13, 19),
            ],
        ),
    (
        'LLLHHLLLE|HLLLEELLL',
        '(?=(L{1}E{2}L{2,3}))',
        [
            slice(13, 18),
            slice(13, 19),
            ],
        ),
    (
        'LLLHHLLLE|HLLLEELLL',
        '(?=(L{2}))',
        [
            slice(0, 2),
            slice(1, 3),
            slice(5, 7),
            slice(6, 8),
            slice(11, 13),
            slice(12, 14),
            slice(16, 18),
            slice(17, 19),
            ],
        ),
    (
        'LLLHHLLLE|HLLLEELLL',
        'L{2}',
        [
            slice(0, 2),
            slice(5, 7),
            slice(11, 13),
            slice(16, 18),
            ],
        ),
    (
        'LLLHHLLLE|HLLLEELLL',
        'H{1,2}',
        [
            slice(3, 4),
            slice(4, 5),
            slice(10, 11),
            slice(3, 5),
            ],
        ),
    ]


@pytest.mark.parametrize(
    'sequence, regex_string, expected',
    # ORDER MATHERS!
    REGEX_SEARCH_w_RANGE + [
        (
            'LLLHHLLLE|HLLLEELLL',
            'LHHL',
            [slice(2, 6)],
            ),
        (
            'LLLHHLLLE|HLLLEELLL',
            'LL',
            [
                slice(0, 2),
                slice(5, 7),
                slice(11, 13),
                slice(16, 18),
                ],
            ),
        (
            'LLLHHLLLE|HLLLEELLL',
            '(?=(LL))',
            [
                slice(0, 2),
                slice(1, 3),
                slice(5, 7),
                slice(6, 8),
                slice(11, 13),
                slice(12, 14),
                slice(16, 18),
                slice(17, 19),
                ],
            ),
        ]
    )
def test_regex_search(sequence, regex_string, expected):
    """
    Test regex_search in single core.

    Single core is important because results are ordered.
    """
    result = regex_search(sequence, regex_string, ncores=1)
    assert len(result) == len(expected)
    for i, j in zip(result, expected):
        assert i == j, (i, j)


@pytest.mark.parametrize(
    'sequence, regex_string, expected',
    REGEX_SEARCH_w_RANGE,
    )
def test_regex_ranges(sequence, regex_string, expected):
    """Test regex_ranges function in multicore."""
    result = regex_range(sequence, regex_string, ncores=None)
    assert len(result) == len(expected)
    # here because is multicore the results are note sorted in order
    for i in result:
        assert i in expected, i
