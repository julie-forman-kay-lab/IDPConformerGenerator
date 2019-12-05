"""Test lib fragment."""
import copy

import pytest

from idpconfgen import Path
from idpconfgen.libs import libfragment as LF

from . import tcommons


class TestResidueAngle:
    """Test ResidueAngle."""

    def test_init(self):
        LF.ResidueAngle()

    def test_repr(self):
        result = repr(LF.ResidueAngle())
        expected = (
            'ResidueAngle('
            'pdbid=None, '
            'residue=None, '
            'dssp=None, '
            'pdb_res_1=None, '
            'pdb_res_2=None, '
            'pdb_res_3=None, '
            'phi=None, '
            'psi=None, '
            'omega=None'
            ')'
            )
        assert result == expected

    def test_equality(self):
        a = LF.ResidueAngle(
            pdbid='XXZ',
            residue='R',
            dssp='L',
            phi=10.0,
            psi=11.0,
            omega=12.0,
            )
        b = LF.ResidueAngle(
            pdbid='XXZ',
            residue='R',
            dssp='L',
            phi=10.0,
            psi=11.0,
            omega=12.0,
            )
        assert a == b

    def test_string(self):
        a = LF.ResidueAngle(
            pdbid='XXZ',
            residue='R',
            dssp='L',
            phi=0.0,
            psi=1.0,
            omega=3.141,
            )
        s = '  XXZ  R L None None None    0.000   57.296  179.966'
        # this is definitively dangerous!!
        assert s == str(a)


class TestFragmentABC:
    def test_init(self):
        with pytest.raises(TypeError):
            LF.FragmentDBABC()


class TestFragmentAngleDB:
    def test_init(self):
        LF.FragmentAngleDB()

    def test_static_read_text(self):
        data = LF.FragmentAngleDB.read_text_file(
            Path(tcommons.data_folder, 'LVALL_sample')
            )
        assert len(data) == 3
        assert len(data[0]) == 12 
        assert len(data[1]) == 10 
        assert len(data[2]) == 9
        assert all(isinstance(i, str) for b in data for i in b)
        assert len(data[0][0].split()) == 9

    def tests_static_parse_raw_data(self):
        """Test data parsing to fragment blocks."""
        data = LF.FragmentAngleDB.read_text_file(
            Path(tcommons.data_folder, 'LVALL_sample')
            )

        parsed_data = LF.FragmentAngleDB._parse_raw_data(data)
        assert len(parsed_data[0]) == 12
        assert all(isinstance(i, LF.ResidueAngle) for i in parsed_data[0])
        assert isinstance(parsed_data[0][0].phi, float)
        assert isinstance(parsed_data[0][0].psi, float)
        assert isinstance(parsed_data[0][0].omega, float)
    
    @pytest.mark.parametrize(
        'fname',
        [
            (Path(tcommons.data_folder, 'LVALL_sample')),
            (Path(tcommons.data_folder, 'fragment_angle_db.pickle')),
            ],
        )
    def test_from_file(self, fname):
        """Test read from file."""
        fragdb = LF.FragmentAngleDB.from_file(fname)
        assert isinstance(fragdb, LF.FragmentAngleDB)

    @pytest.mark.parametrize(
        'fname',
        [
            (Path(tcommons.data_folder, 'LVALL_sample')),
            (Path(tcommons.data_folder, 'fragment_angle_db.pickle')),
            ],
        )
    def test_property_db(self, fname):
        fragdb = LF.FragmentAngleDB.from_file(fname)
        assert fragdb._db is fragdb.db
        assert len(fragdb) == len(fragdb.db)
        fragdb2 = copy.deepcopy(fragdb)
        assert fragdb2 == fragdb
        assert fragdb is not fragdb2

    @pytest.mark.parametrize(
        'fname',
        [
            (Path(tcommons.data_folder, 'LVALL_sample')),
            (Path(tcommons.data_folder, 'fragment_angle_db.pickle')),
            ],
        )
    def test_get_pure_fragment(self, fname):
        fragdb = LF.FragmentAngleDB.from_file(fname)
        
        fragment = fragdb.get_pure_fragment()
        assert isinstance(fragment, list)
        assert isinstance(fragment[0], LF.ResidueAngle)

    @pytest.mark.parametrize(
        'fname',
        [
            (Path(tcommons.data_folder, 'LVALL_sample')),
            (Path(tcommons.data_folder, 'fragment_angle_db.pickle')),
            ],
        )
    def test_get_angle_fragment(self, fname):
        fragdb = LF.FragmentAngleDB.from_file(fname)
        
        fragment = fragdb.get_angle_fragment(fragsize=1)
        assert isinstance(fragment, list)
        assert isinstance(fragment[0], LF.ResidueAngle)
        assert len(fragment) == 1

    @pytest.mark.parametrize(
        'fname',
        [
            (Path(tcommons.data_folder, 'LVALL_sample')),
            (Path(tcommons.data_folder, 'fragment_angle_db.pickle')),
            ],
        )
    def test_get_angle_fragment_0size(self, fname):
        fragdb = LF.FragmentAngleDB.from_file(fname)
        
        fragment = fragdb.get_angle_fragment(fragsize=0)
        assert isinstance(fragment, list)
        assert len(fragment) == 0

    def test_save_pickle_db(self):
        a = LF.FragmentAngleDB.from_file(
            Path(tcommons.data_folder, 'LVALL_sample')
            )
        
        dest_pickle = Path(tcommons.folder_output, 'dummypickle')
        a.save_pickle_db(dest_pickle)
        b = LF.FragmentAngleDB.from_file(dest_pickle)
        assert a == b
