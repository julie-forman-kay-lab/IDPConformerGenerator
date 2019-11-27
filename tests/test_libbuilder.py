"""Test libbuilder."""
import pytest

from idpconfgen.core import definitions as DEFS
from idpconfgen.libs import libbuilder as LB


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
        assert ct.coords.size == len(in1) * len(DEFS.backbone_atoms) + 1

    def test_property_coords_AttributeError(self):
        ct = LB.ConformerTemplate('MSME')
        with pytest.raises(AttributeError):
            ct.coords = 1
