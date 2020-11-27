import warnings
from string import digits

from idpconfgen.core import build_definitions as BD
from idpconfgen.core.definitions import aa3to1
from idpconfgen.libs.libpdb import atom_element, atom_name

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

    required_keys = set([
        'mass',
        'element',
        'class',
        'sigma',
        'epsilon',
        ])

    for key in ff14SB.keys():
        if key.startswith('protein-'):
            assert set(ff14SB[key].keys()) == required_keys, f'Wrong for {key}'


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


def test_ff14SB_atoms_have_params(ff14SB):
    """Test atoms in residues have parameters."""
    required_keys = set(['charge', 'type'])

    for key, residue in ff14SB.items():
        if key[0].isupper():  # this should refer only to residues
            for atom_name, aparams in residue.items():
                assert set(aparams.keys()) == required_keys, f'Wrong for {key}'


def test_ff14SB_residue_atoms_names(ff14SB):
    """Test all atoms in FF are connected to the system."""
    for key, residue in ff14SB.items():
        if key[0].isupper():
            if len(key) == 3 and key[0].isupper():
                try:
                    s1 = set(residue.keys())
                    k1 = aa3to1[key]
                    s2 = set(BD.atom_labels_amber[k1])

                    diff = s1.difference(s2)
                    assert not diff, (key, diff)

                    s2 = set(BD.atom_labels_amber[k1])

                    diff = s1.difference(s2)
                    assert not diff, (key, diff)

                except KeyError:
                    warnings.warn(UserWarning(f'{key} is ignored'))


def test_residue_atoms_pdb_amber():
    """Confirm residue atom names.

    Quantities should be the same for amber and pdb atom names.
    """
    for key, value in BD.atom_labels_pdb.items():
        assert len(value) == len(BD.atom_labels_amber[key]), \
            f'There are different number of atoms for {key}'


def test_amber_templates_element():
    """Test elements in Amber template files are correct."""
    for pdb in BD._amber_pdbs:
        with open(pdb, 'r') as fin:
            for line in fin:
                from_name = line[atom_name].strip()[0]
                element = line[atom_element].strip()
                assert from_name == element, f'Elements differ for {pdb.name}'


def test_pdb_templates_element():
    """Test elements in PDB template files are correct."""
    for pdb in BD._sidechain_template_files:
        with open(pdb, 'r') as fin:
            for line in fin:
                from_name = line[atom_name].strip().lstrip(digits)[0]
                element = line[atom_element].strip()
                assert from_name == element, f'Elements differ for {pdb.name}'
