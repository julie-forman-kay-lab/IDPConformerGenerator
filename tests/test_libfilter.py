"""Test libfilter functions."""
import re

import numpy as np
import pytest

from idpconfgen.libs.libfilter import (
    REGEX_RANGE,
    aligndb,
    make_any_overlap_regex,
    make_helix_overlap_regex,
    make_loop_overlap_regex,
    make_overlap_regex,
    make_ranges,
    make_regex_combinations,
    make_strand_overlap_regex,
    regex_forward_no_overlap,
    regex_forward_with_overlap,
    regex_has_overlap,
    regex_range,
    regex_search,
    )


@pytest.mark.parametrize(
    'regex,expected',
    [
        ('(?=(L{2,6}))', True),
        ('(?=(L{,6}))', True),
        ('(?=(L{6,}))', True),
        ('(?=(L{6}))', True),
        ('(?=([L]{6}))', True),
        ('(?=([L]{1,6}))', True),
        ('(?=([LH]{2,6}))', True),
        ('(?=([LHE]{,6}))', True),
        ('(?=([LHE]{3,6}))', True),
        ('(?=(LLLH))', False),
        ('L+', False),
        ('L{1,5}', True),
        ]
    )
def test_REX_RANGE_regex(regex, expected):
    """Test global variable regex range."""
    result = REGEX_RANGE.findall(regex)
    assert bool(result) == expected


@pytest.mark.parametrize(
    'in1,out',
    [
        (
            {
                "pdb1": {
                    "dssp": "XLLLLX",
                    "fasta": "XMMMMX",
                    "omega": [2, 2, 2, 2, 0],
                    "psi": [0, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 0],
                    "resids": "1,2,3,4,5",
                    },
                },
            (
                {"pdb1": slice(0, 4, None)},
                np.array([
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLL",
                "MMMM",
                ),
            ),
        (
            {
                "pdb1": {
                    "dssp": "XLLLLX",
                    "fasta": "XMMMMX",
                    "omega": [2, 2, 2, 2, 0],
                    "phi": [3, 3, 3, 3, 0],
                    "psi": [0, 1, 1, 1, 1],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "XHHHHX",
                    "fasta": "XAAAAX",
                    "omega": [5, 5, 5, 5, 0],
                    "phi": [6, 6, 6, 6, 0],
                    "psi": [0, 4, 4, 4, 4],
                    "resids": "7,8,9,10,11,12",
                    },
                },
            (
                {
                    "pdb1": slice(0, 4, None),
                    "pdb2": slice(5, 9, None),
                    },
                np.array([
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [np.nan, np.nan, np.nan],
                    [5, 6, 4],
                    [5, 6, 4],
                    [5, 6, 4],
                    [5, 6, 4],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLL|HHHH",
                "MMMM|AAAA",
                ),
            ),
        (
            {
                "pdb1": {
                    "dssp": "XLLLLX",
                    "fasta": "XMMMMX",
                    "omega": [2, 2, 2, 2, 0],
                    "phi": [3, 3, 3, 3, 0],
                    "psi": [0, 1, 1, 1, 1],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "XHHHHX",
                    "fasta": "XAAAAX",
                    "omega": [5, 5, 5, 5, 0],
                    "phi": [6, 6, 6, 6, 0],
                    "psi": [0, 4, 4, 4, 4],
                    "resids": "7,8,9,10,11,12",
                    },
                "pdb3": {
                    "dssp": "XEEEEX",
                    "fasta": "XWWWWX",
                    "omega": [8, 8, 8, 8, 0],
                    "phi": [9, 9, 9, 9, 0],
                    "psi": [0, 7, 7, 7, 7],
                    "resids": "13,14,15,16,17,18",
                    },
                },
            (
                {
                    "pdb1": slice(0, 4, None),
                    "pdb2": slice(5, 9, None),
                    "pdb3": slice(10, 14, None),
                    },
                np.array([
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [np.nan, np.nan, np.nan],
                    [5, 6, 4],
                    [5, 6, 4],
                    [5, 6, 4],
                    [5, 6, 4],
                    [np.nan, np.nan, np.nan],
                    [8, 9, 7],
                    [8, 9, 7],
                    [8, 9, 7],
                    [8, 9, 7],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLL|HHHH|EEEE",
                "MMMM|AAAA|WWWW",
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
                    "dssp": "xLLLLX",
                    "fasta": "xMMMMx",
                    "omega": [2, 2, 2, 2, 0],
                    "psi": [0, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 0],
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
                    "pdb1": slice(0, 4, None),
                    },
                np.array([
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLL",
                "MMMM",
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
                    "dssp": "xLLLLx",
                    "fasta": "xMMMMx",
                    "omega": [2, 2, 2, 2, 0],
                    "psi": [0, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 0],
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
                    "pdb1": slice(0, 4, None),
                    },
                np.array([
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [2, 3, 1],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLL",
                "MMMM",
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
    result = regex_forward_no_overlap(seq, regex)

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
        (
            r'(?=(L{1,6}E{1,10}L{1,6}))',
            20,
            [range(1, 7), range(1, 11), range(1, 7)],
            ['L', 'E', 'L'],
            ),
        (
            r'(?=([LHE]{1,5}))',
            20,
            [range(1, 6)],
            ['LHE'],  # add this to the other function, is passing alist
            ),
        (
            r'(?=([LHE]{1,5}K{3}))',
            20,
            [range(1, 6), range(3, 4)],
            ['LHE', 'K'],
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
            [range(1, 6)],
            ['LHE'],
            r'(?=(',
            r'))',
            [
                r'(?=([LHE]{1}))',
                r'(?=([LHE]{2}))',
                r'(?=([LHE]{3}))',
                r'(?=([LHE]{4}))',
                r'(?=([LHE]{5}))',
                ]
            ),
        (
            [range(1, 6), range(3, 4)],
            ['LHE', 'K'],
            r'(?=(',
            r'))',
            [
                r'(?=([LHE]{1}K{3}))',
                r'(?=([LHE]{2}K{3}))',
                r'(?=([LHE]{3}K{3}))',
                r'(?=([LHE]{4}K{3}))',
                r'(?=([LHE]{5}K{3}))',
                ]
            ),
        (
            [range(1, 6)],
            'LHE',
            r'(?=(',
            r'))',
            [
                r'(?=([LHE]{1}))',
                r'(?=([LHE]{2}))',
                r'(?=([LHE]{3}))',
                r'(?=([LHE]{4}))',
                r'(?=([LHE]{5}))',
                ]
            ),
        (
            [range(1, 2), range(5, 9)],
            ['L', 'E'],
            r'(?=(',
            r'))',
            [
                r'(?=(L{1}E{5}))',
                r'(?=(L{1}E{6}))',
                r'(?=(L{1}E{7}))',
                r'(?=(L{1}E{8}))',
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
        ('(?=(L))', True),
        ('(?=([LHE]{2,5}))', True),
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
    (
        'LLLHHLLLE|HLLLEELLL',
        'K{1,2}',
        [],
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


@pytest.mark.skip(reason="something halts this test")
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


@pytest.mark.parametrize(
    's,pair,expected',
    [
        ("L", (1, 5), "(?=([L]{1,5}))"),
        ("H", (2, 8), "(?=([H]{2,8}))"),
        ("E", (3, 10), "(?=([E]{3,10}))"),
        ("X", (4, 5), "(?=([X]{4,5}))"),
        ("CDF", (4, 5), "(?=([CDF]{4,5}))"),
        ]
    )
def test_makeoverlap_regex(s, pair, expected):
    """Test make overlap regex."""
    result = make_overlap_regex(s, pair)
    assert result == expected


@pytest.mark.parametrize(
    'pair,expected',
    [
        ((1, 5), "(?=([L]{1,5}))"),
        ((2, 8), "(?=([L]{2,8}))"),
        ]
    )
def test_makeoverlap_regex_loop(pair, expected):
    """Test make overlap regex."""
    result = make_loop_overlap_regex(pair)
    assert result == expected


@pytest.mark.parametrize(
    'pair,expected',
    [
        ((1, 5), "(?=([H]{1,5}))"),
        ((2, 8), "(?=([H]{2,8}))"),
        ]
    )
def test_makeoverlap_regex_helix(pair, expected):
    """Test make overlap regex."""
    result = make_helix_overlap_regex(pair)
    assert result == expected


@pytest.mark.parametrize(
    'pair,expected',
    [
        ((1, 5), "(?=([E]{1,5}))"),
        ((4, 8), "(?=([E]{4,8}))"),
        ]
    )
def test_makeoverlap_regex_stand(pair, expected):
    """Test make overlap regex."""
    result = make_strand_overlap_regex(pair)
    assert result == expected


@pytest.mark.parametrize(
    'pair,expected',
    [
        ((1, 5), "(?=([LHE]{1,5}))"),
        ((4, 8), "(?=([LHE]{4,8}))"),
        ]
    )
def test_makeoverlap_regex_any(pair, expected):
    """Test make overlap regex."""
    result = make_any_overlap_regex(pair)
    assert result == expected


@pytest.mark.parametrize(
    's,pair',
    [
        ("L", (5, 1)),
        ("H", (-2, 8)),
        ]
    )
def test_makeoverlap_regex_with_error(s, pair):
    """Test make overlap regex."""
    with pytest.raises(ValueError):
        make_overlap_regex(s, pair)
