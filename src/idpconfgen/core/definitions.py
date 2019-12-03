"""Static definitions that serve the whole program infrastructure."""
from argparse import Namespace
from pathlib import Path

import numpy as np

_file_path = Path(__file__).resolve()
data_folder = Path(_file_path.parents[1], 'data')


class ReprClean:
    """Baseclass for repr method."""
    def __repr__(self):
        kwargs = \
            ', '.join(f'{key}={val!r}' for key, val in self.__dict__.items())
        rpr = '{}({})'.format(__class__.__name__, kwargs)
        return rpr



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


aa3to1 = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'}

aa1to3 = {v:k for k, v in aa3to1.items()}

CA_radius = 1.7
C_radius = 1.7
N_radius = 1.55
O_radius = 1.52

#back bone seed coordinates
N_seed = np.array([0.000, 0.000, 0.000], dtype=np.float32)
CA_seed = np.array([1.458, 0.000, 0.000], dtype=np.float32)
C_seed = np.array([2.009, 1.420, 0.000], dtype=np.float32)

N_name = 'N'
CA_name = 'CA'
C_name = 'C'
O_name = 'O'
COO_name = 'X'


backbone_atoms = (N_name, CA_name, C_name, O_name)
num_bb_atoms = len(backbone_atoms)
