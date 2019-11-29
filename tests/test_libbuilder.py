"""Test libbuilder."""
import numpy as np
import pytest

from idpconfgen import Path
from idpconfgen.core import definitions as DEFS
from idpconfgen.libs import libbuilder as LB

from . import tcommons

class TestConformerTemplate:
    
    @pytest.mark.parametrize(
        'in1,expected',
        [
            ('QWERTYIPAS', (
                'GLN',
                'TRP',
                'GLU',
                'ARG',
                'THR',
                'TYR',
                'ILE',
                'PRO',
                'ALA',
                'SER',
                )),
            (['M', 'A', 'V'], ('MET', 'ALA', 'VAL')),
            (['MET', 'ALA', 'VAL'], ('MET', 'ALA', 'VAL')),
            ],
        )
    def test_parse_seq(self, in1, expected):
        assert LB.ConformerTemplate._parse_seq(in1) == expected

    @pytest.mark.parametrize(
        'in1,error',
        [
            ('MAVO', ValueError),
            (['MET', 'AAA'], ValueError),
            ],
        )
    def test_parse_errors(self, in1, error):
        with pytest.raises(error):
            LB.ConformerTemplate._parse_seq(in1)

    def test_property_seq(self):
        ct = LB.ConformerTemplate('MARVEL')
        assert ct.seq == ('MET', 'ALA', 'ARG', 'VAL', 'GLU', 'LEU')
        assert isinstance(ct.seq, tuple)
    
    def test_property_seq_AttributeError(self):
        ct = LB.ConformerTemplate('MARVEL')
        with pytest.raises(AttributeError):
            ct.seq = 1
    
    @pytest.mark.parametrize(
        'in1',
        [
            ('M'),
            ('MARVEL'),
            ('MEGAPRTEIN')
            ],
        )
    def test_property_coords(self, in1):
        ct = LB.ConformerTemplate(in1)
        num_of_residues = len(in1) * len(DEFS.backbone_atoms) + 1
        assert ct.coords.size == num_of_residues * 3
        assert ct.coords.shape == (num_of_residues, 3)

    def test_property_coords_AttributeError(self):
        ct = LB.ConformerTemplate('MSME')
        with pytest.raises(AttributeError):
            ct.coords = 1

    def test_property_atomnames(self):
        ct = LB.ConformerTemplate('MARVEL')
        assert ct.atomnames.shape == (ct.coords.shape[0],)

    def test_property_atomnames_AttributeError(self):
        ct = LB.ConformerTemplate('MARVEL')
        with pytest.raises(AttributeError):
            ct.atomnames = 1
    
    @pytest.mark.parametrize(
        'resind,aname,expected',
        [
            (0, 'N', 0),
            (0, 'CA', 1),
            (10, 'C', 42),
            (321321, DEFS.COO_name, -1),
            ],
        )
    def test_get_index(self, resind, aname, expected):
        result = LB.ConformerTemplate._get_index(resind, aname)
        assert result == expected
   
    @pytest.mark.parametrize(
        'resind,aname,error',
        [
            (0, 'Z', ValueError),
            ],
        )
    def test_get_index_Errors(self, resind, aname, error):
        """
        1: tests against a wrong atom name.
        """
        with pytest.raises(error):
            LB.ConformerTemplate._get_index(resind, aname)

    @pytest.mark.parametrize(
        'resindex,atomname,coords,realindex',
        [
            (0, 'N', [1.000, 1.000, 1.000], 0),
            (5, 'CA', [1.123, 2.134, 0.321], 5 * 4 + 1),
            (3, 'O', [0.001, 0.230, 0.543], 3 * 4 + 3),
            (5, DEFS.COO_name, [1.000, 1.030, 1.040], -1),
            (0, DEFS.COO_name, [1.000, 1.030, 1.040], -1),
            ],
        )
    def test_add_coords(self,resindex, atomname, coords, realindex):
        ct = LB.ConformerTemplate('MARVEL')

        coords = np.array(coords)
        ct.add_atom_coords(resindex, atomname, coords)
        assert all(np.equal(ct.coords[realindex], coords))

    def test_add_coords_resindx_error(self):
        ct = LB.ConformerTemplate('MARVEL')
        with pytest.raises(IndexError):
            ct.add_atom_coords(10, 'N', np.array([1., 1., 1.]))

    @pytest.mark.parametrize(
        'seq',
        [
            ('MARVEL'),
            ],
        )
    def test_is_complete(self, seq):
        ct = LB.ConformerTemplate(seq)
        for i in range(len(seq)):
            for c in DEFS.backbone_atoms:
                ct.add_atom_coords(i, c, np.array([1.000, 1.000, 1.000]))
        else:
            ct.add_atom_coords(0, DEFS.COO_name, np.array([1.00, 1.00, 1.00]))

        assert ct.is_complete()
    
    def test_is_complete_manual(self):
        ct = LB.ConformerTemplate('MAV')
        ct._coords[:] = 1
        assert ct.is_complete()

    def test_is_complete_false_1(self):
        ct = LB.ConformerTemplate('MAV')
        assert not ct.is_complete()

    def test_is_complete_false_2(self):
        ct = LB.ConformerTemplate('MAV')
        ct._coords[:3, :] = 1
        assert not ct.is_complete()
    
    @pytest.mark.parametrize(
        'resind,aname',
        [
            (0, 'N'),
            (0, DEFS.COO_name),
            ],
        )
    def test_get_coords(self, resind, aname):
        ct = LB.ConformerTemplate('MAV')
        coords = np.array([1., 2., 0.7])
        ct.add_atom_coords(resind, aname, coords)
        result = ct.get_coord(resind, aname)
        assert result.shape == coords.shape
        assert np.all(np.equal(result, coords))


class TestBuilder:
    
    def test_init(self):
        LB.ConformerBuilderNeRF(None, None)

    
def test_Conformer_Builder_integration_1():
    conf = LB.ConformerTemplate('MARVEL')
    builder = LB.ConformerBuilderNeRF(conf, None)
    conf.add_atom_coords(0, 'N', np.array([1., 1., 1.,]))
    coords1 = conf.coords
    conf = LB.ConformerTemplate('MARVEL')
    builder = LB.ConformerBuilderNeRF(conf, None)
    assert coords1 is not conf.coords
    assert not np.all(np.equal(coords1, conf.coords))


class TestFragLoopDB:
    def test_init(self):
        LB.FragmentAngleDB()

    def test_static_read_text(self):
        data = LB.FragmentAngleDB.read_text_file(
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
        data = LB.FragmentAngleDB.read_text_file(
            Path(tcommons.data_folder, 'LVALL_sample')
            )

        parsed_data = LB.FragmentAngleDB._parse_raw_data(data)
        assert len(parsed_data[0]) == 12
        assert len(parsed_data[0][0]) == 6
        assert all(isinstance(i, LB.ResidueAngle) for i in parsed_data[0])
        assert isinstance(parsed_data[0][0].phi, float)
        assert isinstance(parsed_data[0][0].psi, float)
        assert isinstance(parsed_data[0][0].omega, float)
    
    @pytest.mark.parametrize(
        'fname',
        [
            (Path(tcommons.data_folder, 'LVALL_sample')),
            ],
        )
    def test_from_file(self, fname):
        """Test read from file."""
        fragdb = LB.FragmentAngleDB.from_file(fname)
        assert isinstance(fragdb, LB.FragmentAngleDB)


class TestResidueAngleTuple:
    """Test ResidueAngle Tuple."""    

    resang = LB.ResidueAngle('1XXX', 'A', 'L', 1.0, 2.0, 3.0)

    def test_attributes(self):
        assert self.resang.pdbid == '1XXX'
        assert self.resang.letter == 'A'
        assert self.resang.dssp == 'L'
        assert self.resang.phi == 1.0
        assert self.resang.psi == 2.0
        assert self.resang.omega == 3.0

