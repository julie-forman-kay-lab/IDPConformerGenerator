"""Test core definitions."""
import pytest

from idpconfgen.core import definitions as DEFS



def test_aa_letters():
    assert len(DEFS.aa3to1) == len(DEFS.aa1to3)
    assert list(DEFS.aa3to1.keys()) == list(DEFS.aa1to3.values())
    assert list(DEFS.aa1to3.keys()) == list(DEFS.aa3to1.values())


@pytest.mark.parametrize(
    'in1,in2',
    [
        (DEFS.N_seed, [0.000, 0.000, 0.000]),
        (DEFS.CA_seed, [1.458, 0.000, 0.000]),
        (DEFS.C_seed, [2.009, 1.420, 0.000]),
        ],
    )
def test_bb_seed(in1, in2):
    assert all(i - j < 0.00001 for i, j in zip(list(in1), in2))


def test_bb_atoms_len():
    assert len(DEFS.backbone_atoms) == 4


def test_bb_atoms():
    assert DEFS.backbone_atoms == ('N', 'CA', 'C', 'O')
