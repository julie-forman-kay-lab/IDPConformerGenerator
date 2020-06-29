"""Test libparse module."""
from pathlib import Path as Path_

import pytest

from idpconfgen import Path, log
from idpconfgen.core.exceptions import DSSPParserError
from idpconfgen.libs import libparse

from . import tcommons


def test_t2s_dispacher():
    """
    Test type2string dispacher.

    All values should return a string.
    """
    types = {
        str: 'hello',
        bytes: b'world',
        type(Path()): Path(__file__),
        type(Path_()): Path(__file__),
        }

    for k, v in types.items():
        assert isinstance(libparse.type2string[k](v), str)


@pytest.mark.parametrize(
    'in1, expected',
    [
        (
            'LLLHEEE',
            [['L', slice(0, 3)], ['H', slice(3, 4)], ['E', slice(4, 7)]]
            ),
        ]
    )
def test_group_by(in1, expected):
    """Group sequence by equality."""
    result = libparse.group_by(in1)
    assert result == expected


def test_group_runs_1():
    """Test group runs."""
    num = [-2, -1, 0, 1, 2, 3, 5, 6, 7]
    expected = [
        [-2, -1, 0, 1, 2, 3],
        [5, 6, 7],
        ]

    for i, result in enumerate(libparse.group_runs(num)):
        assert result == expected[i]
    else:
        assert i == 1


def test_group_runs_2():
    """Test group runs."""
    num = [-2, -1, 0, 1, 2, 3, 5, 6, 7]
    expected = [num]

    for i, result in enumerate(libparse.group_runs(num, tolerance=2)):
        assert result == expected[i]
    else:
        assert i == 0


def test_parse_dssp():
    """Parse DSSP from example."""
    expected = [
        (
            b'  S      T ',
            b'SRMPSPPMPSS',
            [str(i).encode() for i in (list(range(47, 56)) + [348, 349])],
            ),
        (
            b'   EEEE ',
            b'AAPRLSFL',
            [str(i).encode() for i in range(395, 403)],
            ),
        (
            b'    ',
            b'XXXX',
            [str(i).encode() for i in range(501, 505)],
            ),
        ]

    with tcommons.example_dssp.open('rb') as dssp_data:
        for i, result in enumerate(libparse.parse_dssp(dssp_data.read())):
            assert result == expected[i]
        else:
            assert i == 2


def test_parse_dssp_IndexError():
    """
    Test IndexError.

    Obtained in case `dssp` file as not the `#` sperator.
    """
    with open(Path(tcommons.data_folder, 'wrong.dssp'), 'rb') as fin:
        data = fin.read()
        with pytest.raises(DSSPParserError):
            for _ in libparse.parse_dssp(data):
                pass


@pytest.mark.parametrize(
    'in1, expected',
    [
        ('A', {'A', 'a'}),
        ('aa', {'AA', 'Aa', 'aA', 'aa'}),
        ('a_', {'A_', 'a_'}),
        ]
    )
def test_case_combinations(in1, expected):
    """Test cases combinations."""
    result = libparse.sample_case(in1)
    assert result == expected


@pytest.fixture
def BTE_A_results():
    """
    Return the expected dictionary for split_pdb_by_dssp.

    Considers using the `reduced` kwarg.
    """
    results = {
        "1BTE_A_seg0": {
            "dssp": "LLLEEEEEELLHHHHLLLLEEEEELLL",
            "fasta": "ETQECLFFNANWERDRTNQTGVEPCYG",
            "resids": "7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33",  # noqa: E501
            },
        "1BTE_A_seg1": {
            "dssp": "LLEEEEEEEEELLEEEEEEEEEELLLLLLLLLLLEEELLLLLLLEEEEELLLLHHHLEEELLLL",  # noqa: E501
            "fasta": "RRHCFATWKNISGSIEIVKQGCWLDDINCYDRTDCIEKKDSPEVYFCCCEGNMCNEKFSYFPEM",  # noqa: E501
            "resids": "38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101"  # noqa: E501
            },
        }
    return results


def test_split_pdb_by_dssp(BTE_A_results):
    """Test mkdssp with split."""
    sgen = libparse.split_pdb_by_dssp(
        Path(tcommons.data_folder, '1BTE_A.pdb'),
        Path(tcommons.data_folder, '1BTE_A.dssp').read_bytes(),
        reduced=True,
        )

    _perform_dssp_split(sgen, BTE_A_results)


def _perform_dssp_split(sgen, expected):
    """
    Separate logic to test PDB split by DSSP.

    Logic is separated because there are two entry points:
        * the direct function: split_pdb_by_dssp
        * by mkdssp: mkdssp_w_split
    """
    counter = 0
    for fout, d, pdbpart in sgen:
        with open(tcommons.data_folder / f'{fout}.pdb', 'rb') as PART:
            assert pdbpart == PART.read()

        assert expected[fout] == d
        counter += 1

    assert counter == 2

    with pytest.raises(StopIteration):
        next(sgen)


def test_split_pdb_by_dssp_minimum():
    """Test mkdssp with split."""
    results = {
        # this seg0 actually corresponds to seg1 if `minimum` is not used
        "1BTE_A_seg0": {
            "dssp": "LLEEEEEEEEELLEEEEEEEEEELLLLLLLLLLLEEELLLLLLLEEEEELLLLHHHLEEELLLL",  # noqa: E501
            "fasta": "RRHCFATWKNISGSIEIVKQGCWLDDINCYDRTDCIEKKDSPEVYFCCCEGNMCNEKFSYFPEM",  # noqa: E501
            "resids": "38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101"  # noqa: E501
            },
        }

    sgen = libparse.split_pdb_by_dssp(
        Path(tcommons.data_folder, '1BTE_A.pdb'),
        Path(tcommons.data_folder, '1BTE_A.dssp').read_bytes(),
        reduced=True,
        minimum=40,
        )

    counter = 0
    for fout, d, pdbpart in sgen:
        with open(tcommons.data_folder / '1BTE_A_seg1.pdb', 'rb') as PART:
            assert pdbpart == PART.read()
        counter += 1
        assert results[fout] == d

    assert counter == 1  # make sure generator had only 1 element

    with pytest.raises(StopIteration):
        next(sgen)


def test_mkdssp_w_split(BTE_A_results):
    """."""
    try:
        sgen = libparse.mkdssp_w_split(
            Path(tcommons.data_folder, '1BTE_A.pdb'),
            'mkdssp',
            reduced=True,
            )
    except FileNotFoundError:
        # in case mkdssp is not installed in the running machine
        return

    _perform_dssp_split(sgen, BTE_A_results)


def test_pop_dict_differece():
    """Pops dict difference keys."""
    d1 = dict.fromkeys(range(10))
    d2 = dict.fromkeys(range(5))
    libparse.pop_difference_with_log(d1, d2)
    assert d1.keys() == d2.keys()


def test_pop_dict_differece_2():
    """Pops dict difference keys."""
    d1 = dict.fromkeys(range(10))
    d2 = dict.fromkeys(range(10))
    libparse.pop_difference_with_log(d1, d2)
    assert d1.keys() == d2.keys()
