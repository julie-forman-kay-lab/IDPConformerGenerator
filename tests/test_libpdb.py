"""Test libpdb."""
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs import libpdb
from idpconfgen.libs.libpdb import (
    PDBID,
    PDBIDFactory,
    PDBList,
    is_pdb,
    )

from .tcommons import pdb_example


@pytest.mark.parametrize(
    'is_func,structure,expected',
    [
        (is_pdb, pdb_example, True),
        # we need a False example here, review this function.
        ]
    )
def test_is_function(is_func, structure, expected):
    """Test wether file is pdb."""
    assert is_func(structure.read_text()) == expected


@pytest.mark.parametrize(
    'atom,element,expected',
    [
        ('N', 'N', ' N  '),
        ('CO', 'C', ' CO '),
        ('CA', 'CA', 'CA  '),
        ('CHA', 'C', ' CHA'),
        ('HD11', 'H', 'HD11'),
        ('FE', 'FE', 'FE  '),
        ('OD1', 'O', ' OD1'),
        ]
    )
def test_format_atoms(atom, element, expected):
    """Test format atom."""
    assert libpdb.format_atom_name(atom, element) == expected


def test_format_atoms_error():
    """Test error with wrong atom."""
    with pytest.raises(EXCPTS.PDBFormatError):
        libpdb.format_atom_name('FES', 'FE')


def test_format_chain():
    """Test format chain."""
    assert libpdb.format_chainid('AAA') == 'A'


@given(st.text())
def test_nothing(txt):
    """Test nothing function."""
    assert txt == libpdb._nothing(txt)


class TestPDBID:
    """Test PDBID class."""

    def test_PDBID1(self):
        """Simplest test."""
        pdbid = PDBID('1ABC', chain='D')
        assert pdbid.name == '1ABC'

    def test_PDBID_null_chain(self):
        """Test without chain id."""
        pdbid = PDBID('1ABC')
        assert pdbid.name == '1ABC'
        assert pdbid.chain is None

    def test_PDBID_str(self):
        """Test parsing underscore."""
        pdbid = PDBID('1ABC', chain='Z')
        assert str(pdbid) == '1ABC_Z'

    def test_PDBID_str2(self):
        """Test string conversion without chain."""
        pdbid = PDBID('1ABC')
        assert str(pdbid) == '1ABC'

    def test_repr(self):
        """Test repr with chain."""
        pdbid = PDBID('1ABC', chain='Z')
        assert repr(pdbid) == "PDBID(name='1ABC', chain='Z')"

    def test_repr2(self):
        """Test repr without chain."""
        pdbid = PDBID('1ABC')
        assert repr(pdbid) == "PDBID(name='1ABC')"

    def test_equality(self):
        """Test equality."""
        pdbid1 = PDBID('1ABC', chain='Z')
        pdbid2 = PDBID('1ABC', chain='Z')
        assert pdbid1 == pdbid2

    def test_lower(self):
        """Test comparison: lower."""
        pdbid1 = PDBID('1ABC', chain='X')
        pdbid2 = PDBID('1ABC', chain='Z')
        assert pdbid1 < pdbid2

    def test_higher(self):
        """Test comparison: higher."""
        pdbid1 = PDBID('1ABC', chain='Z')
        pdbid2 = PDBID('1ABC', chain='X')
        assert pdbid1 > pdbid2


class TestPDBIDFactory:
    """Test PDBIDFactory."""

    def test_parse_XXXX(self):
        """Test XXXX parsing."""
        pdbid, chain = PDBIDFactory._parse_XXXX('1ABC')
        assert (pdbid, chain) == ('1ABC', None)

    def test_parse_XXXXC(self):
        """Test XXXXC parsing."""
        pdbid, chain = PDBIDFactory._parse_XXXXC('1ABCYY2 ')
        assert (pdbid, chain) == ('1ABC', 'YY2')

    def test_parse_XXXX_C(self):
        """Test XXXX_C parsing."""
        pdbid, chain = PDBIDFactory._parse_XXXX_C('1ABC_UUI    something ')
        assert (pdbid, chain) == ('1ABC', 'UUI')

    def test_parse1(self):
        """Test string PDBID input."""
        pdbid0 = PDBIDFactory('1ABC')
        pdbid1 = PDBID('1ABC')
        assert pdbid0 == pdbid1

    def test_parse2(self):
        """Test trailing whitespace."""
        pdbid0 = PDBIDFactory('1ABC ')
        pdbid1 = PDBID('1ABC')
        assert pdbid0 == pdbid1

    def test_parse3(self):
        """Test chain ID one character."""
        pdbid0 = PDBIDFactory('1ABCX')
        pdbid1 = PDBID('1ABC', chain='X')
        assert pdbid0 == pdbid1

    def test_parse4(self):
        """Test long string with more information."""
        pdbid0 = PDBIDFactory('1ABCX            something')
        pdbid1 = PDBID('1ABC', chain='X')
        assert pdbid0 == pdbid1

    def test_parse5(self):
        """Test 3 chars chain ID."""
        pdbid0 = PDBIDFactory('1ABCXYZ')
        pdbid1 = PDBID('1ABC', chain='XYZ')
        assert pdbid0 == pdbid1

    def test_parse6(self):
        """Test newline char."""
        pdbid0 = PDBIDFactory('1ABCXYZ\n')
        pdbid1 = PDBID('1ABC', chain='XYZ')
        assert pdbid0 == pdbid1

    def test_parse7(self):
        """Test path string."""
        pdbid0 = PDBIDFactory('some/folder/PDBs/1ABC_D.pdb')
        pdbid1 = PDBID('1ABC', chain='D')
        assert pdbid0 == pdbid1

    def test_parse8(self):
        """Test path Path."""
        pdbid0 = PDBIDFactory(Path('some', 'fldr', 'PDBs', '1ABC_DEF.pdb'))
        pdbid1 = PDBID('1ABC', chain='DEF')
        assert pdbid0 == pdbid1

    def test_parse_error_1(self):
        """Test parser not found."""
        with pytest.raises(EXCPTS.PDBIDFactoryError):
            PDBIDFactory('#*(!@#!@')


class TestPDBList:
    """Test PDBList class."""

    def test_1(self):
        """Test single element."""
        pdblist = PDBList(('1ABC',))
        assert pdblist[0] == PDBID('1ABC')

    def test_repr(self):
        """Test repr."""
        pdblist = PDBList(('1ABC',))
        r = repr(pdblist)
        expected = "PDBList(\n    PDBID(name='1ABC'))\n"
        assert r == expected

    def test_str(self):
        """Test string representation."""
        pdblist = PDBList(('1ABC',))
        s = str(pdblist)
        assert s == 'PDBList with 1 element(s).'

    def test_1_2(self):
        """Test intantiating from cls."""
        pdblist = PDBList(('1ABC',))
        pdb2 = PDBList(pdblist)
        assert pdb2 is pdblist

    def test_2(self):
        """Test multiple elements of different nature."""
        names = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            )

        pdblist = PDBList(names)
        pdbids = [PDBIDFactory(name) for name in names]
        assert len(pdblist) == len(names)
        assert pdblist == set(pdbids)

    def test_3(self):
        """Test with comments."""
        names = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            '# some comment',
            )

        pdblist = PDBList(names)
        pdbids = [PDBIDFactory(name) for name in names[:-1]]
        assert len(pdblist) == len(names) - 1
        assert pdblist == set(pdbids)

    def test_difference(self):
        """Test difference between two PDBLists."""
        names1 = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            )
        names2 = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            )

        pdblist1 = PDBList(names1)
        pdblist2 = PDBList(names2)

        pdblist3 = PDBList(names1[-1:])

        assert pdblist1.difference(pdblist2) == pdblist3

    def test_to_tuple(self):
        """Test conversion of PDBList to tuple."""
        names_sorted = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            )

        names_unsorted = (
            '2ABCZ',
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            '1ABC',
            )

        pdblist1 = PDBList(names_sorted)
        pdblist2 = PDBList(names_unsorted)

        assert pdblist1 == pdblist2
        assert isinstance(pdblist2.to_tuple(), tuple)
        assert tuple(pdblist1) == pdblist2.to_tuple()

    def test_write(self):
        """Test writing PDBList."""
        names_sorted = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            )
        pdblist1 = PDBList(names_sorted)
        fout = Path('pdblistfout.list')
        pdblist1.write(fout)
        result = fout.read_text()
        fout.unlink()
        expected = (
            '1ABC\n'
            '2ABC_Z\n'
            '3RFC_YYY\n'
            '3RFC_ZZZ'
            )

        assert result == expected
