

import numpy as np
import pytest

from idpconfgen.core.definitions import (
    C_ele_name,
    N_ele_name,
    O_ele_name,
    C_vdw_radius,
    N_vdw_radius,
    O_vdw_radius,
    C_alpha_pdb_label,
    C_carbonyl_pdb_label,
    N_amide_pdb_label,
    O_carbonyl_pdb_label,
    O2_carboxyl_pdb_label,
    C_alpha,
    C_carbonyl,
    N_amide,
    O_carbonyl,
    O2_carboxyl,
    N_seed,
    CA_seed,
    C_seed,
    backbone_atoms,
    num_of_backbone_atoms,
    )


@pytest.mark.parametrize(
    'ele, atoms',
    [
        (C_ele_name, (C_alpha, C_carbonyl)),
        (N_ele_name, (N_amide,)),
        (O_ele_name, (O_carbonyl, O2_carboxyl)),
        ]
    )
def test_atoms_are_same_element(ele, atoms):
    """Test atoms are of correct element."""
    assert all(ele == atom.element_name for atom in atoms)


@pytest.mark.parametrize(
    'vdw, atoms',
    [
        (C_vdw_radius, (C_alpha, C_carbonyl)),
        (N_vdw_radius, (N_amide,)),
        (O_vdw_radius, (O_carbonyl, O2_carboxyl)),
        ]
    )
def test_atoms_VDW_radius(vdw, atoms):
    """Test atoms are of correct element."""
    assert all(vdw == atom.vdw for atom in atoms)


@pytest.mark.parametrize(
    'pdb_label,atom',
    [
        (C_alpha_pdb_label, C_alpha),
        (C_carbonyl_pdb_label, C_carbonyl),
        (N_amide_pdb_label, N_amide),
        (O_carbonyl_pdb_label, O_carbonyl),
        (O2_carboxyl_pdb_label, O2_carboxyl),
        ]
    )
def test_atoms_pdb_label(pdb_label, atom):
    """Test atoms PDB label."""
    assert pdb_label == atom.pdb_label


@pytest.fixture(
    params=[
        N_seed,
        CA_seed,
        C_seed,
        ]
    )
def XYZ_seed_coordinates(request):
    """XYZ seeding coordinates."""
    return request.param


def test_XYZ_seed_coords_type(XYZ_seed_coordinates):
    """Test seed coords are numpy.ndarray."""
    assert isinstance(XYZ_seed_coordinates, np.ndarray)


def test_XYZ_seed_coords_dtype(XYZ_seed_coordinates):
    """Test seed coords numpy.ndarray dtype is np.float32."""
    assert XYZ_seed_coordinates.dtype == np.float32


def test_backbone_atoms():
    """Test tuple."""
    assert isinstance(backbone_atoms, tuple)


def test_num_bb_atoms():
    """Test is integer."""
    assert isinstance(num_of_backbone_atoms, int)
