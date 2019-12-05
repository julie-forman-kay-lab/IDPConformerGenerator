"""Test lib fragment."""
import pytest

from idpconfgen import Path
from idpconfgen.libs import libfragment as LF

from . import tcommons


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
