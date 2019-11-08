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
from idpconfgen.libs import libparse, libpdb

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
    dssp_file2 = tcommons.data_folder / '1ABC_E.dssp'
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

    def test_dssp9(self):
        """Test FASTA extraction."""
        dobj = libparse.DSSPParser(fin=self.dssp_file)
        expected = ['K', 'K', 'V', 'K', 'V', 'S', 'H', 'R', 'S', 'H']
        assert list(dobj.fasta) == expected

    def test_dssp10(self):
        """Test equality with pdbid None."""

        dobj1 = libparse.DSSPParser(fin=self.dssp_file)
        dobj2 = libparse.DSSPParser(fin=self.dssp_file)

        assert dobj1 == dobj2

    def test_dssp11(self):
        """Test equality with pdbid."""

        dobj1 = libparse.DSSPParser(fin=self.dssp_file, pdbid='1ABC')
        dobj2 = libparse.DSSPParser(fin=self.dssp_file, pdbid='1ABC')

        assert dobj1 == dobj2
    
    def test_dssp12(self):
        """Test differ with pdbid different."""

        dobj1 = libparse.DSSPParser(fin=self.dssp_file, pdbid='1ABC')
        dobj2 = libparse.DSSPParser(fin=self.dssp_file, pdbid='1XXX')

        assert dobj1 != dobj2

    def test_dssp13(self):
        """Test differ with data different."""

        dobj1 = libparse.DSSPParser(fin=self.dssp_file, pdbid='1ABC')
        dobj2 = libparse.DSSPParser(fin=self.dssp_file2, pdbid='1ABC')

        assert dobj1 != dobj2


class TestDSSPMediator:
    
    dssp_file = tcommons.data_folder / '1ABC_D.dssp'
    pdb_file = tcommons.data_folder / '1A12.pdb'
    dsspdata = dssp_file.read_text()

    def test_dsspmediator_0(self):
        
        dssp_mediated = libparse.DSSPMulticoreMediator((
            self.pdb_file,
            self.dsspdata,
            ))
        assert isinstance(dssp_mediated, libparse.DSSPParser)
    
    def test_dsspmediator_1(self):
        
        dssp_parser = libparse.DSSPParser(
            data=self.dsspdata,
            pdbid=libpdb.PDBIDFactory(self.pdb_file),
            )

        dssp_mediated = libparse.DSSPMulticoreMediator((
            self.pdb_file,
            self.dsspdata,
            ))

        assert dssp_parser == dssp_mediated


class Test_dssp_ss_saver:
    
    dssp_fileD = tcommons.data_folder / '1ABC_D.dssp'
    dssp_fileE = tcommons.data_folder / '1ABC_E.dssp'
    
    obj1 = libparse.DSSPParser(
        fin=dssp_fileD,
        pdbid='1ABC_D',
        )
    obj2 = libparse.DSSPParser(
        fin=dssp_fileE,
        pdbid='1ABC_E',
        )
    
    def test_concat_dsspparsers_1(self):
        """Test concatenation of DSSPParsers."""

        output = libparse._concatenate_ss_from_dsspparsers([
            self.obj1,
            self.obj2,
            ])

        expected = ['1ABC_D| HGIBE TS ', '1ABC_E| HGIBE TS']
        assert expected == output

    def test_export_ss(self):
        libparse.export_ss_from_DSSP(self.obj1, self.obj2,
            output=tcommons.folder_output / 'dssp.database'
            )
    
        ofile = Path(tcommons.folder_output / 'dssp.database')
        results = ofile.read_text()

        expected = '1ABC_D| HGIBE TS \n1ABC_E| HGIBE TS\n'
        assert expected == results

