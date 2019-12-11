"""Test lib fragment."""
import copy

import pytest

from idpconfgen import Path
from idpconfgen.libs import libfragment as LF

from . import tcommons


class TestResidueAngle:
    """Test ResidueAngle."""
    
    @pytest.fixture
    def residue_angle_a(self):
        """Example fixture of a ResidueAngle object."""
        return LF.ResidueAngle(
            pdbid='XXZ',
            residue='R',
            dssp='L',
            phi=10.0,
            psi=11.0,
            omega=12.0,
            )
    
    @pytest.fixture
    def residue_angle_b(self):
        """Example fixture of a ResidueAngle object."""
        return LF.ResidueAngle(
            pdbid='XXZ',
            residue='R',
            dssp='L',
            phi=10.0,
            psi=11.0,
            omega=12.0,
            )

    @pytest.fixture
    def residue_angle_c(self):
        """Example fixture of a ResidueAngle object."""
        return LF.ResidueAngle(
            pdbid='AAZ',
            residue='R',
            dssp='L',
            phi=10.0,
            psi=11.0,
            omega=12.0,
            )

    def test_init_TypeError(self):
        """Test init class with empty args raises TypeError."""
        with pytest.raises(TypeError):
            LF.ResidueAngle()
    
    def test_equality(self, residue_angle_a, residue_angle_b):
        """Test equality between two objects."""
        assert residue_angle_a == residue_angle_b

    def test_inequality(self, residue_angle_a, residue_angle_c):
        assert residue_angle_a != residue_angle_c

    def test_string(self, residue_angle_a):
        s = '  XXZ  R L None None None  572.958  630.254  687.549'
        # this is definitively dangerous!!
        assert s == str(residue_angle_a)


class TestFragmentAngleDB:
    """Test FragmentAngleDB."""
    
    @pytest.fixture
    def file_LVALL_sample(self):
        """Path to a angle db TXT file."""
        return Path(tcommons.data_folder, 'LVALL_sample')

    def test_init(self):
        """Test class can be initiated empty."""
        LF.FragmentAngleDB()
    
    @pytest.fixture
    def data_text_file(self, file_LVALL_sample):
        """Read the example file into a list of lists of strings."""
        return LF.FragmentAngleDB.read_text_file(file_LVALL_sample)
    
    def test_static_read_text_len_1(self, data_text_file):
        """Test number of fragments of strings."""
        assert len(data_text_file) == 3

    def test_static_read_text_len_2(self, data_text_file):
        """Test length of the first fragment of strings."""
        assert len(data_text_file[0]) == 12

    def test_static_read_text_len_3(self, data_text_file):
        """Test length of the second fragment of strings."""
        assert len(data_text_file[1]) == 10

    def test_static_read_text_len_4(self, data_text_file):
        """Test length of the third fragment of strings."""
        assert len(data_text_file[2]) == 9

    def test_static_read_text_5(self, data_text_file):
        """Test all items in fragment of strings are strings."""
        assert all(isinstance(i, str) for b in data_text_file for i in b)
    
    def test_static_read_text_6(self, data_text_file):
        """Test length of residue string after split."""
        assert len(data_text_file[0][0].split()) == 9
    
    @pytest.fixture
    def fragmentdb(self, file_LVALL_sample):
        """Read example file to a FragmentAngleDB object."""
        return LF.FragmentAngleDB.from_file(file_LVALL_sample)
    
    def test_fragnentdb_db_0(self, fragmentdb):
        """Test class type."""
        assert isinstance(fragmentdb, LF.FragmentAngleDB)

    def test_fragnentdb_db_1(self, fragmentdb):
        """Test size of the first fragment."""
        assert len(fragmentdb[0]) == 12
    
    def test_fragnentdb_db_2(self, fragmentdb):
        """Test all elements of a fragment are ResidueAngle instances."""
        assert all(isinstance(i, LF.ResidueAngle) for i in fragmentdb[0])

    def test_fragnentdb_db_3(self, fragmentdb):
        """Test phi attribute is read as float."""
        assert isinstance(fragmentdb[0][0].phi, float)
    
    def test_fragnentdb_db_4(self, fragmentdb):
        """Test psi attribute is read as float."""
        assert isinstance(fragmentdb[0][0].psi, float)

    def test_fragnentdb_db_5(self, fragmentdb):
        """Test omega attribute is read as float."""
        assert isinstance(fragmentdb[0][0].omega, float)

    def test_fragmentdb_db_6(self, fragmentdb):
        """Test .db attribute property."""
        assert fragmentdb._db is fragmentdb.db
  
    def test_fragmentdb_db_7(self, fragmentdb):
        """Test length evaluation of FragmentAngleDB."""
        assert len(fragmentdb) == len(fragmentdb.db)

    def test_fragmentdb_db_8(self, fragmentdb):
        """Test deep copy executes properly."""
        fragmentdb2 = copy.deepcopy(fragmentdb)
        assert fragmentdb2 == fragmentdb
        assert fragmentdb2 is not fragmentdb
    
    def test_fragmentdb_db_9(self, fragmentdb):
        """Test total length of the fragmentdb object."""
        assert fragmentdb.total_len() == 31 

    def test_StopIteration(self, fragmentdb):
        """Test fragmentdb is iterable with StopIteration implemented."""
        with pytest.raises(StopIteration):
            for i in range(100):
                next(fragmentdb)

    def test_for_loop(self, fragmentdb):
        """Test all items in fragmentdb are lists."""
        for i in fragmentdb:
            assert isinstance(i, list)

    def test_get_pure_fragment_1(self, fragmentdb):
        """Test a pure fragment is of list type."""
        fragment = fragmentdb.get_pure_fragment()
        assert isinstance(fragment, list)
        
    def test_get_pure_fragment_2(self, fragmentdb):
        """Test get_pure_fragment API."""
        fragment = fragmentdb.get_pure_fragment()
        assert isinstance(fragment[0], LF.ResidueAngle)
   
    def test_get_angle_fragment_1(self, fragmentdb):
        """Test angle fragment API for fragsize=None."""
        fraglist = fragmentdb.get_angle_fragment()
        assert isinstance(fraglist, list)
        
    def test_get_angle_fragment_2(self, fragmentdb):
        """Test angle fragment API for fragsize=4."""
        fraglist = fragmentdb.get_angle_fragment(fragsize=4)
        assert isinstance(fraglist, list)
        assert len(fraglist) == 4

    @pytest.fixture
    def dest_pickle_path(self):
        """Destination test pickle file."""
        return Path(tcommons.folder_output, 'dummypickle')

    def test_save_text_to_pickle(self, fragmentdb, dest_pickle_path):
        """Test saving a DB to pickle."""
        fragmentdb.save_pickle_db(dest_pickle_path)
    
    def test_reload_pickle(self, fragmentdb, dest_pickle_path):
        """Test if a pickle saved from txt equals txt when reloaded."""
        b = LF.FragmentAngleDB.from_file(dest_pickle_path)
        assert fragmentdb == b
