from idpconfgen.core import build_definitions as BD
from idpconfgen.core.definitions import aa3to1

import pytest


@pytest.fixture
def ff14SB():
    """The fs14SB dictionary."""
    return BD.read_ff14SB_params()


def test_ff14SB_params_length(ff14SB):
    """Test reads ff14SB parameters."""
    assert len(ff14SB) == 114


def test_ff14SB_params_random_1(ff14SB):
    """Test some random value."""
    assert ff14SB['protein-C']['element'] == 'C'


def test_ff14SB_params_random_2(ff14SB):
    """Test some random value."""
    assert ff14SB['protein-CT']['class'] == 'CT'


def test_ff14SB_params_atom_keys(ff14SB):
    """Test atom type keys are okay."""
    for key in ff14SB.keys():
        if key.startswith('protein-'):
            assert set(ff14SB[key].keys()) == set([
                'mass', 'element', 'class', 'sigma', 'epsilon',
                ])

def test_ff14SB_coloumb_and_lj(ff14SB):
    """Assert presence of keys."""
    assert 'coulomb14scale' in ff14SB
    assert 'lj14scale' in ff14SB


def test_ff14SB_has_all_residues(ff14SB):
    """Assert forcefield contains all residues used."""
    s1 = set(ff14SB.keys())
    s2 = set(aa3to1.keys())
    s2.remove('HIS')
    assert s2.issubset(s1)

