"""Test libpdb."""
from pathlib import Path

from idpconfgen.libs.libpdb import (
    PDBID,
    PDBList,
    PDBIDFactory,
    )

class TestPDBID:
    
    def test_PDBID1(self):
        pdbid = PDBID('1ABC', chain='D')
        assert pdbid.name == '1ABC'
    
    def test_PDBID_null_chain(self):
        pdbid = PDBID('1ABC')
        assert pdbid.name == '1ABC'
        assert pdbid.chain is None
    
    def test_PDBID_str(self):
        pdbid = PDBID('1ABC', chain='Z')
        assert str(pdbid) == '1ABC_Z'
        
    def test_PDBID_str2(self):
        pdbid = PDBID('1ABC')
        assert str(pdbid) == '1ABC'
    
    def test_repr(self):
        pdbid = PDBID('1ABC', chain='Z')
        assert repr(pdbid) == "PDBID(name='1ABC', chain='Z')"
    
    def test_repr2(self):
        pdbid = PDBID('1ABC')
        assert repr(pdbid) == "PDBID(name='1ABC')"
    
    def test_equality(self):
        pdbid1 = PDBID('1ABC', chain='Z')
        pdbid2 = PDBID('1ABC', chain='Z')
        assert pdbid1 == pdbid2
    
    def test_lower(self):
        pdbid1 = PDBID('1ABC', chain='X')
        pdbid2 = PDBID('1ABC', chain='Z')
        assert pdbid1 < pdbid2
    
    def test_higher(self):
        pdbid1 = PDBID('1ABC', chain='Z')
        pdbid2 = PDBID('1ABC', chain='X')
        assert pdbid1 > pdbid2



class TestPDBIDFactory:
    
    def test_parse1(self):
        pdbid0 = PDBIDFactory('1ABC')
        pdbid1 = PDBID('1ABC')
        assert pdbid0 == pdbid1
    
    def test_parse2(self):
        pdbid0 = PDBIDFactory('1ABC ')
        pdbid1 = PDBID('1ABC')
        assert pdbid0 == pdbid1
    
    def test_parse3(self):
        pdbid0 = PDBIDFactory('1ABCX')
        pdbid1 = PDBID('1ABC', chain='X')
        assert pdbid0 == pdbid1
    
    def test_parse4(self):
        pdbid0 = PDBIDFactory('1ABCX            something')
        pdbid1 = PDBID('1ABC', chain='X')
        assert pdbid0 == pdbid1
    
    def test_parse5(self):
        pdbid0 = PDBIDFactory('1ABCXYZ')
        pdbid1 = PDBID('1ABC', chain='XYZ')
        assert pdbid0 == pdbid1
    
    def test_parse6(self):
        pdbid0 = PDBIDFactory('1ABCXYZ\n')
        pdbid1 = PDBID('1ABC', chain='XYZ')
        assert pdbid0 == pdbid1
    
    def test_parse7(self):
        pdbid0 = PDBIDFactory('some/folder/PDBs/1ABC_D.pdb')
        pdbid1 = PDBID('1ABC', chain='D')
        assert pdbid0 == pdbid1

    def test_parse8(self):
        pdbid0 = PDBIDFactory(Path('some', 'fldr', 'PDBs', '1ABC_DEF.pdb'))
        pdbid1 = PDBID('1ABC', chain='DEF')
        assert pdbid0 == pdbid1


class TestPDBList:
    """Test PDBList class."""

    def test_1(self):
        """Test single element."""
        pdblist = PDBList(('1ABC',))
        assert pdblist[0] == PDBID('1ABC')
    
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