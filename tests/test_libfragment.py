"""Test lib fragment."""
import copy
from typing import NamedTuple

import pytest

from idpconfgen import Path
from idpconfgen.libs import libfragment as LF

from . import tcommons


class TestResidueAngle:
    """Test ResidueAngle."""
    
    @pytest.fixture
    def residue_angle_a(self):
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
        return LF.ResidueAngle(
            pdbid='XXZ',
            residue='R',
            dssp='L',
            phi=10.0,
            psi=11.0,
            omega=12.0,
            )

    def test_init(self):
        with pytest.raises(TypeError):
            LF.ResidueAngle()
    
    def test_equality(self, residue_angle_a, residue_angle_b):
        assert residue_angle_a == residue_angle_b

    def test_string(self, residue_angle_a):
        s = '  XXZ  R L None None None  572.958  630.254  687.549'
        # this is definitively dangerous!!
        assert s == str(residue_angle_a)


class TestFragmentAngleDB:
    """Test FragmentAngleDB."""
    
    @pytest.fixture
    def file_LVALL_sample(self):
        return Path(tcommons.data_folder, 'LVALL_sample')

    def test_init(self):
        """Test class can be initiated empty."""
        LF.FragmentAngleDB()
    
    @pytest.fixture
    def data_text_file(self, file_LVALL_sample):
        return LF.FragmentAngleDB.read_text_file(file_LVALL_sample)
    
    def test_static_read_text_len_1(self, data_text_file):
        assert len(data_text_file) == 3

    def test_static_read_text_len_2(self, data_text_file):
        assert len(data_text_file[0]) == 12

    def test_static_read_text_len_3(self, data_text_file):
        assert len(data_text_file[1]) == 10

    def test_static_read_text_len_4(self, data_text_file):
        assert len(data_text_file[1]) == 10

    def test_static_read_text_len_5(self, data_text_file):
        assert len(data_text_file[2]) == 9

    def test_static_read_text_6(self, data_text_file):
        assert all(isinstance(i, str) for b in data_text_file for i in b)
    
    def test_static_read_text_7(self, data_text_file):
        assert len(data_text_file[0][0].split()) == 9
    
    @pytest.fixture
    def fragmentdb(self, file_LVALL_sample):
        return LF.FragmentAngleDB.from_file(file_LVALL_sample)
    
    def test_fragnentdb_db_0(self, fragmentdb):
        assert isinstance(fragmentdb, LF.FragmentAngleDB)

    def test_fragnentdb_db_1(self, fragmentdb):
        assert len(fragmentdb[0]) == 12
    
    def test_fragnentdb_db_2(self, fragmentdb):
        assert all(isinstance(i, LF.ResidueAngle) for i in fragmentdb[0])

    def test_fragnentdb_db_3(self, fragmentdb):
        assert isinstance(fragmentdb[0][0].phi, float)
    
    def test_fragnentdb_db_4(self, fragmentdb):
        assert isinstance(fragmentdb[0][0].psi, float)

    def test_fragnentdb_db_5(self, fragmentdb):
        assert isinstance(fragmentdb[0][0].omega, float)

    def test_fragmentdb_db_6(self, fragmentdb):
        assert fragmentdb._db is fragmentdb.db
  
    def test_fragmentdb_db_7(self, fragmentdb):
        assert len(fragmentdb) == len(fragmentdb.db)

    def test_fragmentdb_db_8(self, fragmentdb):
        fragmentdb2 = copy.deepcopy(fragmentdb)
        assert fragmentdb2 == fragmentdb
        assert fragmentdb2 is not fragmentdb
    
    def test_fragmentdb_db_9(self, fragmentdb):
        assert fragmentdb.total_len() == 31 

    def test_StopIteration(self, fragmentdb):
        with pytest.raises(StopIteration):
            for i in range(100):
                next(fragmentdb)

    def test_for_loop(self, fragmentdb):
        for i in fragmentdb:
            assert isinstance(i, list)

    #@pytest.mark.parametrize(
    #    'fname',
    #    [
    #        (Path(tcommons.data_folder, 'LVALL_sample')),
    #        #(Path(tcommons.data_folder, 'fragment_angle_db.pickle')),
    #        ],
    #    )
    #def test_get_pure_fragment(self, fname):
    #    fragdb = LF.FragmentAngleDB.from_file(fname)
    #    
    #    fragment = fragdb.get_pure_fragment()

    def test_get_pure_fragment_1(self, fragmentdb):
        fragment = fragmentdb.get_pure_fragment()
        assert isinstance(fragment, list)
        
    def test_get_pure_fragment_2(self, fragmentdb):
        fragment = fragmentdb.get_pure_fragment()
        assert isinstance(fragment[0], LF.ResidueAngle)
   
    def test_get_angle_fragment_1(self, fragmentdb):
        fraglist = fragmentdb.get_angle_fragment()
        assert isinstance(fraglist, list)
        
    def test_get_angle_fragment_2(self, fragmentdb):
        fraglist = fragmentdb.get_angle_fragment(fragsize=4)
        assert isinstance(fraglist, list)
        assert len(fraglist) == 4

    #@pytest.fixture(
    #    params=[
    #        Path(tcommons.data_folder, 'LVALL_sample'),
    #        #Path(tcommons.data_folder, 'fragment_angle_db.pickle'),
    #        ],
    #    )
    #def frag_file_path(self, request):
    #    return request.param
    #
    #@pytest.fixture(params=range(1, 3))
    #def frag_size(self, request):
    #    return request.param
    #
    #@pytest.fixture()#(params=[(frag_file_path, frag_size)])
    #def angle_fragment(self, frag_file_path, frag_size):#request):
    #    fragdb = LF.FragmentAngleDB.from_file(frag_file_path)
    #    return fragdb.get_angle_fragment(frag_size)
    #
    #def test_get_angle_fragment_type(self, angle_fragment):
    #    assert isinstance(angle_fragment, list)
    #
    #def test_get_angle_fragment_item_type(self, angle_fragment):
    #    assert isinstance(angle_fragment[0], LF.ResidueAngle)

    #def test_get_angle_fragment_item_len(self, angle_fragment):
    #    assert len(angle_fragment) == 1
    #
    #@pytest.fixture(params=[0])
    #def frag_size_null(self, request):
    #    return request.param

    #@pytest.fixture()
    #def angle_fragment_null(self, frag_file_path, frag_size_null):
    #    fragdb = LF.FragmentAngleDB.from_file(frag_file_path)
    #    return fragdb.get_angle_fragment(frag_size_null)

    #def test_get_angle_fragment_null_type(self, angle_fragment_null):
    #    assert isinstance(angle_fragment_null, list)

    #def test_get_angle_fragment_null_len(self, angle_fragment_null):
    #    assert len(angle_fragment_null) == 0
    #@pytest.mark.parametrize(
    #    'fname',
    #    [
    #        (Path(tcommons.data_folder, 'LVALL_sample')),
    #        (Path(tcommons.data_folder, 'fragment_angle_db.pickle')),
    #        ],
    #    )
    #def test_get_angle_fragment(self, fname):
    #    fragdb = LF.FragmentAngleDB.from_file(fname)
    #    
    #    fragment = fragdb.get_angle_fragment(fragsize=1)
    #    assert isinstance(fragment, list)
    #    assert isinstance(fragment[0], LF.ResidueAngle)
    #    assert len(fragment) == 1

    #@pytest.fixture(
    #    params=[
    #        Path(tcommons.data_folder, 'LVALL_sample'),
    #        Path(tcommons.data_folder, 'fragment_angle_db.pickle'),
    #        ],
    #    )
    #def fragment(self, request):
    #    fragdb = LF.FragmentAngleDB.from_file(request.param)
    #    fragment = fragdb.get_angle_fragment(fragsize=0)
    #    return fragment

    #def test_get_angle_fragment_0size_type(self, fragment):
    #    assert isinstance(fragment, list)

    #def test_get_angle_fragment_0size_len(self, fragment):
    #    assert len(fragment) == 0


class TestFragDB2pickle:
    """Test writing DB to pickle and readig it back."""
    
    @pytest.fixture(scope='class')
    def LVALL_sample(self):
        return LF.FragmentAngleDB.from_file(
            Path(tcommons.data_folder, 'LVALL_sample')
            )
  
    @pytest.fixture(scope='class')
    def dest_pickle_path(self):
        return Path(tcommons.folder_output, 'dummypickle')

    def test_save_text_to_pickle(self, LVALL_sample, dest_pickle_path):
        LVALL_sample.save_pickle_db(dest_pickle_path)
    
    def test_reload_pickle(self, LVALL_sample, dest_pickle_path):
        """Test if a pickle saved from txt equals txt when reloaded."""
        b = LF.FragmentAngleDB.from_file(dest_pickle_path)
        assert LVALL_sample == b

    #def test_save_pickle_db(self, dest_pickle_path):
    #    b = LF.FragmentAngleDB.from_file(dest_pickle_path)
    #    c = LF.FragmentAngleDB.from_file(
    #        Path(tcommons.data_folder, 'fragment_angle_db.pickle')
    #        )
    #    assert b == c
