"""
Test data file parsers.

Parsers are coded in libparse.

Currently tests:

    - DSSPParser
"""
import pytest

from idpconfgen import Path
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.core import definitions as DEFS
from idpconfgen.libs import libparse

from . import tcommons


class TestLITA():
    
    data = [
        'i like to program in python',
        'u like to program in java',
        'w like to program in csharp',
        ]

    def test_list_index_to_array_1(self):
        """Test index."""
        parsed = libparse.list_index_to_array(self.data, index=0)
        assert list(parsed) == ['i', 'u', 'w']

    def test_list_index_to_array_2(self):
        """Test slice."""
        parsed = libparse.list_index_to_array(self.data, start=21, stop=25)
        assert list(parsed) == ['pyth', 'java', 'csha']

    def test_list_index_to_array_3(self):
        """Test ValueError."""
        with pytest.raises(ValueError):
            libparse.list_index_to_array(self.data)
   
    def test_list_index_to_array_4(self):
        """Test IndexError."""
        with pytest.raises(IndexError):
            libparse.list_index_to_array(self.data, index=400)

    def test_list_index_to_array5(self):
        """Test slice object."""
        parsed = libparse.list_index_to_array(self.data, sObj=slice(21, 25, 1))
        assert list(parsed) == ['pyth', 'java', 'csha']
    
    def test_list_index_to_array6(self):
        """Test slice object."""
        with pytest.raises(ValueError):
            libparse.list_index_to_array(self.data, sObj=5)


class TestDSSPParser:
    """Test DSSPParser."""
    dssp_file = tcommons.data_folder / '1ABC_D.dssp'
    dssp_wrong = tcommons.data_folder / 'wrong.dssp'
    dssp_wrong2 = tcommons.data_folder / 'wrong2.dssp'
    dsspdata = dssp_file.read_text()

    def test_dssp1(self):
        """Test initiation empty args.""" 
        libparse.DSSPParser()
    
    def test_dssp1_1(self):
        """Test pdbid attribute.""" 
        dobj = libparse.DSSPParser(pdbid='1ABC')
        assert dobj.pdbid == '1ABC'

    def test_dssp2(self):
        """Test initiation from data."""
        libparse.DSSPParser(data=self.dsspdata)

    def test_dssp3(self):
        """Test initiation from file."""
        libparse.DSSPParser(fin=self.dssp_file)

    def test_dssp4(self):
        """Test secondary structure parsing."""
        dobj = libparse.DSSPParser(fin=self.dssp_file)

        expected = [' ', 'H', 'G', 'I', 'B', 'E', ' ', 'T', 'S', ' ']
        assert list(dobj.ss) == expected
   
    def test_dssp5(self):
        """Test DSSPParserError."""
        with pytest.raises(EXCPTS.DSSPParserError):
            dobj = libparse.DSSPParser(fin=self.dssp_wrong)

    def test_dssp6(self):
        """Test find index."""
        data = self.dssp_file.read_text().split('\n')
        index = libparse.DSSPParser._finds_data_index(data)
        assert index == 28

    def test_dssp7(self):
        """
        Test wrong ss capturing.

        The input here has the '#' wrongly placed.
        Captures if any of the registered chars does not belong
        to the possible DSSP secondary structure chars.
        """
        with pytest.raises(EXCPTS.DSSPSecStructError):
            libparse.DSSPParser(fin=self.dssp_wrong2)

    def test_dssp8(self):
        """Test file does not exist exception."""
        with pytest.raises(FileNotFoundError):
            libparse.DSSPParser(fin='doesnotexist.dssp')
