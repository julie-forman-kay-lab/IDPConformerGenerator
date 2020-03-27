"""Static definitions that serve the whole program infrastructure."""
from argparse import Namespace
from collections import namedtuple

import numpy as np


# Definition of Atoms required for building conformers
BuildingAtom = namedtuple(
    'BuildingAtom',
    [
        'element_name', # according to Periodic Table
        'vdw',  # Van der Waals
        'pdb_label',  # atom label in the PDB file
        ]
    )

# element names
C_ele_name = 'C'
N_ele_name = 'N'
O_ele_name = 'O'

# https://periodictable.com/Properties/A/VanDerWaalsRadius.an.html
# Atomic radius in pm
C_vdw_radius = 170
N_vdw_radius = 155
O_vdw_radius = 152

# PDB atom labels
C_alpha_pdb_label = 'CA'
C_carbonyl_pdb_label = 'C'
N_amide_pdb_label = 'N'
O_carbonyl_pdb_label = 'O'
O2_carboxyl_pdb_label = 'OXT'

# Atoms used in the building process
# Different atoms may be of the same element type
C_alpha = BuildingAtom(C_ele_name, C_vdw_radius, C_alpha_pdb_label)
C_carbonyl = BuildingAtom(C_ele_name, C_vdw_radius, C_carbonyl_pdb_label)
N_amide = BuildingAtom(N_ele_name, N_vdw_radius, N_amide_pdb_label)
O_carbonyl = BuildingAtom(O_ele_name, O_vdw_radius, O_carbonyl_pdb_label)
O2_carboxyl = BuildingAtom(O_ele_name, O_vdw_radius, O2_carboxyl_pdb_label)

# Building seed coordinates for the backbone of the first residue
N_seed = np.array([0.000, 0.000, 0.000], dtype=np.float32)
CA_seed = np.array([1.458, 0.000, 0.000], dtype=np.float32)
C_seed = np.array([2.009, 1.420, 0.000], dtype=np.float32)

# other general definitions
backbone_atoms = (N_amide, C_alpha, C_carbonyl, O_carbonyl)
num_of_backbone_atoms = len(backbone_atoms)

# keys from https://github.com/cmbi/dssp/blob/7c2942773cd37d47b3e4611597d5e1eb886d95ba/src/dssp.cpp#L66-L74  # noqa:
dssp_ss_keys = Namespace(
    ahelix='H',
    helix_3='G',
    helix_5='I',
    bbridge='B',
    strand='E',
    turn='T',
    bend='S',
    loop=' ',
    )

dssp_ss_keys.all_helix = (
    dssp_ss_keys.ahelix,
    dssp_ss_keys.helix_3,
    dssp_ss_keys.helix_5,
    )

dssp_ss_keys.all_strand = (
    dssp_ss_keys.bbridge,
    dssp_ss_keys.strand,
    )

dssp_ss_keys.all_loops = (
    dssp_ss_keys.turn,
    dssp_ss_keys.bend,
    dssp_ss_keys.loop,
    )
