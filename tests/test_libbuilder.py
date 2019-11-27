"""Test libbuilder."""
import pytest

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
        ct = LB.ConformerTemplate('MAVERL')
        assert ct.seq == ('MET', 'ALA', 'VAL', 'GLU', 'ARG', 'LEU')
        assert isinstance(ct.seq, tuple)
