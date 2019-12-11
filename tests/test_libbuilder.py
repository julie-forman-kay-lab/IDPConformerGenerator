"""Test libbuilder."""
import numpy as np
import pytest


from idpconfgen import Path
from idpconfgen.core import definitions as DEFS
from idpconfgen.libs import libbuilder as LB
from idpconfgen.libs import libfragment as LF

from . import tcommons


class TestAtomNeRF:
    """Test AtomNeRF object."""
    def test_Atom1(self):
        a = LB.AtomNeRF('N', 'UPPER', 1, 0, -1, 0)
        assert a.name == 'N'
        assert a.rosetta_name == 'UPPER'
        assert a.poff == 1
        assert a.xoff == 0
        assert a.yoff == -1
        assert a.resoff == 0


    def test_Atom2_AttributeError(self):
        a = LB.AtomNeRF('N', 'UPPER', 1, 0, -1, 0)
        with pytest.raises(AttributeError):
            a.name = 'Z'
        with pytest.raises(AttributeError):
            a.rosetta_name = 'Z'
        with pytest.raises(AttributeError):
            a.poff = 0
        with pytest.raises(AttributeError):
            a.xoff = 1
        with pytest.raises(AttributeError):
            a.yoff = 2
        with pytest.raises(AttributeError):
            a.resoff = 0


    def test_NeRF_building_order_NtoC(self):
        exp = [
            LB.N_atom_NeRF,
            LB.O_atom_NeRF,
            LB.CA_atom_NeRF,
            LB.C_atom_NeRF,
            ]
        assert exp == LB.NeRF_building_order_NtoC


    def test_NeRF_building_order_CtoN(self):
        exp = [
            LB.C_atom_NeRF_b,
            LB.CA_atom_NeRF_b,
            LB.O_atom_NeRF_b,
            LB.N_atom_NeRF_b,
            ]
        assert exp == LB.NeRF_building_order_CtoN


class TestConformerTemplate:
    """Test base Conformer Template."""

    @pytest.fixture
    def conformer_template_marvel(self):
        return LB.ConformerTemplate('MARVEL')

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
        """Test _parse_seq with different inputs."""
        assert LB.ConformerTemplate._parse_seq(in1) == expected

    @pytest.mark.parametrize(
        'in1,error',
        [
            ('MAVO', ValueError),
            (['MET', 'AAA'], ValueError),
            ],
        )
    def test_parse_errors(self, in1, error):
        """Test errors on _parse_seq."""
        with pytest.raises(error):
            LB.ConformerTemplate._parse_seq(in1)

    def test_property_seq(self, conformer_template_marvel):
        """Test ConformerTemplate seq values."""
        ct = conformer_template_marvel
        assert ct.seq == ('MET', 'ALA', 'ARG', 'VAL', 'GLU', 'LEU')

    def test_seq_type(self, conformer_template_marvel):
        """Test ConformerTemplate seq is tuple type."""
        ct = conformer_template_marvel
        assert isinstance(ct.seq, tuple)
    
    def test_property_seq_AttributeError(self, conformer_template_marvel):
        """Test AttributeError on attribute assignement .seq."""
        with pytest.raises(AttributeError):
            conformer_template_marvel.seq = 1
    
    @pytest.mark.parametrize(
        'in1',
        [
            ('M'),
            ('MARVEL'),
            ('MEGAPRTEIN')
            ],
        )
    def test_property_coords(self, in1):
        """Test property coords."""
        ct = LB.ConformerTemplate(in1)
        num_of_residues = len(in1) * len(DEFS.backbone_atoms) + 1
        assert ct.coords.size == num_of_residues * 3
        assert ct.coords.shape == (num_of_residues, 3)

    def test_property_coords_AttributeError(self):
        """Test AttributeError on attribute assignement .coords."""
        ct = LB.ConformerTemplate('MSME')
        with pytest.raises(AttributeError):
            ct.coords = 1

    def test_property_atomnames(self, conformer_template_marvel):
        """Test atomnames."""
        ct = conformer_template_marvel
        assert ct.atomnames.shape == (ct.coords.shape[0],)

    def test_property_atomnames_AttributeError(self, conformer_template_marvel):
        """Test atomnames AttributeError on assignment."""
        with pytest.raises(AttributeError):
            conformer_template_marvel.atomnames = 1
    
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
        """Test _get_index static method."""
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
        """Test adding coords."""
        ct = LB.ConformerTemplate('MARVEL')

        coords = np.array(coords)
        ct.add_atom_coords(resindex, atomname, coords)
        assert all(np.equal(ct.coords[realindex], coords))

    def test_add_coords_resindx_error(self, conformer_template_marvel):
        """Test IndexError on adding coords outside range."""
        with pytest.raises(IndexError):
            conformer_template_marvel.add_atom_coords(
                10, 'N', np.array([1., 1., 1.]))

    @pytest.mark.parametrize(
        'seq',
        [
            ('MARVEL'),
            ],
        )
    def test_is_complete(self, seq):
        """Test conformer is complete."""
        ct = LB.ConformerTemplate(seq)
        for i in range(len(seq)):
            for c in DEFS.backbone_atoms:
                ct.add_atom_coords(i, c, np.array([1.000, 1.000, 1.000]))
        else:
            ct.add_atom_coords(0, DEFS.COO_name, np.array([1.00, 1.00, 1.00]))

        assert ct.is_complete()
    
    def test_is_complete_manual(self):
        """Test is complete assigning manually."""
        ct = LB.ConformerTemplate('MAV')
        ct._coords[:] = 1
        assert ct.is_complete()

    def test_is_complete_false_1(self):
        """Test is not complete assessement no coords. """
        ct = LB.ConformerTemplate('MAV')
        assert not ct.is_complete()

    def test_is_complete_false_2(self):
        """Test is not complete assessement some coords. """
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


class TestFragmentAngleDBNeRF:
    """
    Test dedicated Fragment Angle DB for NeRF protocol.
    
    Tests on FragmentAngleDB are in tests_libfragment.py
    """
    
    fragdb = LF.FragmentAngleDB.from_file(
        Path(tcommons.data_folder, 'LVALL_sample')
        )
    
    def test_init(self):
        LB.FragmentAngleDBNeRF(self.fragdb)

    @pytest.mark.parametrize(
        'in1,expected',
        [
            (
                [
                    LF.ResidueAngle(
                        pdbid='XXX',
                        residue='A',
                        dssp='L',
                        phi=1.0,
                        psi=2.0,
                        omega=3.0,
                        ),
                    LF.ResidueAngle(
                        pdbid='XXZ',
                        residue='R',
                        dssp='L',
                        phi=4.0,
                        psi=5.0,
                        omega=6.0,
                        ),
                    ],
                {
                    0: {
                        'PHI': 4.0,
                        'PSI': 2.0,
                        'OMEGA': 3.0,
                        },
                    }
                ),
            (
                [
                    LF.ResidueAngle(
                        pdbid='XXX',
                        residue='A',
                        dssp='L',
                        phi=1.0,
                        psi=2.0,
                        omega=3.0,
                        ),
                    LF.ResidueAngle(
                        pdbid='XXZ',
                        residue='R',
                        dssp='L',
                        phi=4.0,
                        psi=5.0,
                        omega=6.0,
                        ),
                    LF.ResidueAngle(
                        pdbid='XXX',
                        residue='A',
                        dssp='L',
                        phi=7.0,
                        psi=8.0,
                        omega=9.0,
                        ),
                    LF.ResidueAngle(
                        pdbid='XXZ',
                        residue='R',
                        dssp='L',
                        phi=10.0,
                        psi=11.0,
                        omega=12.0,
                        ),
                    ],
                {
                    0: {
                        'PHI': 4.0,
                        'PSI': 2.0,
                        'OMEGA': 3.0,
                        },
                    1: {
                        'PHI': 7.0,
                        'PSI': 5.0,
                        'OMEGA': 6.0,
                        },
                    2: {
                        'PHI': 10.0,
                        'PSI': 8.0,
                        'OMEGA': 9.0,
                        },
                    }
                ),
            ],
        )
    def test_transform_frag2dict(self, in1, expected):
        fdb = LB.FragmentAngleDBNeRF(self.fragdb)
        result = fdb.transform_frag2dict(in1)
        assert dict(result) == expected

#class TestConformerNeRF:
#   
#    @pytest.mark.parametrize(
#        'seq',
#        [
#            ('MARVEL'),
#            ],
#        )
#    def test_init(self, seq):
#        LB.ConformerNeRF(seq)
#

#class TestBuilder:
#    
#    def test_init(self):
#        LB.ConformerBuilderNeRF(None, None, None)
#
#    
##def test_Conformer_Builder_integration_1():
##    conf = LB.ConformerTemplate('MARVEL')
##    builder = LB.ConformerBuilderNeRF(conf, None, None)
##    conf.add_atom_coords(0, 'N', np.array([1., 1., 1.,]))
##    coords1 = conf.coords
##    conf = LB.ConformerTemplate('MARVEL')
##    builder = LB.ConformerBuilderNeRF(conf, None, None)
##    assert coords1 is not conf.coords
##    assert not np.all(np.equal(coords1, conf.coords))
#
#
#def test_fragment_ABC():
#    with pytest.raises(TypeError):
#        LB.FragmentDBABC()
#
#
#
#
#class TestResidueAngleTuple:
#    """Test ResidueAngle Tuple."""    
#
#    resang = LB.ResidueAngle('1XXX', 'A', 'L', 1.0, 2.0, 3.0)
#
#    def test_attributes(self):
#        assert self.resang.pdbid == '1XXX'
#        assert self.resang.letter == 'A'
#        assert self.resang.dssp == 'L'
#        assert self.resang.phi == 1.0
#        assert self.resang.psi == 2.0
#        assert self.resang.omega == 3.0
#
#class TestRosettaAtom:
#    """Test RosettaAtom class to represent atoms from Rosetta building DB."""
#
#    ra = LB.RosettaAtomData(1, 2, 3, 4, 5, 6)
#
#    @pytest.mark.parametrize(
#        'attr1,expected',
#        [
#            (ra.polar_theta, 1),
#            (ra.polar_phi, 2),
#            (ra.polar_r, 3),
#            (ra.parent_atom, 4),
#            (ra.xaxis_atom, 5),
#            (ra.yaxis_atom, 6),
#            ],
#        )
#    def test_polar_theta(self, attr1, expected):
#        assert attr1 == expected
