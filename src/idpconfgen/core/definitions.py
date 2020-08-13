"""Static definitions that serve the whole program infrastructure."""
from argparse import Namespace
from collections import namedtuple
from itertools import chain
# does not import the Path from IDPConfgen to avoid circular imports
from pathlib import Path


core_folder = Path(__file__).parent

# Amino-acid 3 to 1 letter code dictionary
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

# Amino-acid 1 to 3 letter code dictionary
aa1to3 = {v: k for k, v in aa3to1.items()}

# amino-acids atom labels
# from: http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
# PDB column
# Taken from PDB entry 6I1B REVDAT 15-OCT-92.
atom_labels = {
    'A': ('N', 'CA', 'C', 'O', 'CB', 'H', 'HA', 'HB1', 'HB2', 'HB3'),
    'C': ('N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', '1HB', '2HB', 'HG'),
    'D': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', 'H', 'HA', '1HB', '2HB'),
    'E': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'H', 'HA', '1HB', '2HB', '1HG', '2HG'),
    'F': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'H', 'HA', '1HB', '2HB', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'),
    'G': ('N', 'CA', 'C', 'O', 'H', '1HA', '2HA'),
    'H': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD2', 'CE1', 'ND1', 'NE2', 'H', 'HA', '1HB', '2HB', 'HD1', 'HD2', 'HE1'),
    'I': ('N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'H', 'HA', 'HB', '1HG1', '2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1'),
    'K': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', 'H', 'HA', '1HB', '2HB', '1HG', '2HG', '1HD','2HD', '1HE', '2HE', '1HZ', '2HZ', '3HZ'),
    'L': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'H', 'HA', '1HB', '2HB', 'HG', '1HD1', '2HD1', '3HD1', '1HD2', '2HD2', '3HD2'),
    'M': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CE', 'SD', 'H', 'HA', '1HB', '2HB', '1HG', '2HG', '1HE', '2HE', '3HE'),
    'N': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'ND2', 'OD1', 'H', 'HA', '1HB', '2HB', '1HD2', '2HD2'),
    'P': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'H2', 'H1', 'HA', '1HB', '2HB', '1HG', '2HG', '1HD', '2HD'),
    'Q': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE2', 'OE1', 'H', 'HA', '1HB', '2HB', '1HG', '2HG', '1HE2', '2HE2'),
    'R': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CZ', 'NE', 'NH1', 'NH2', 'H', 'HA', '1HB', '2HB', '1HG', '2HG', '1HD', '2HD', 'HE', '1HH1', '2HH1', '1HH2', '2HH2'),
    'S': ('N', 'CA', 'C', 'O', 'CB', 'OG', 'H', 'HA', '1HB', '2HB', 'HG'),
    'T': ('N', 'CA', 'C', 'O', 'CB', 'CG2', 'OG1', 'H', 'HA', 'HB', 'HG1', '1HG2', '2HG2', '3HG2'),
    'V': ('N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', '1HG1', '2HG1', '3HG1', '1HG2', '2HG2', '3HG2'),
    'W': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'NE1', 'H', 'HA', '1HB', '2HB', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2'),
    'Y': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'H', 'HA', '1HB', '2HB', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'),
    }

# heavy atoms
heavy_atoms = {'C', 'O', 'N', 'S', 'P'}

#
# https://www.cgl.ucsf.edu/chimerax/docs/user/radii.html
vdW_radii_tsai_1999 = {
'C': 1.7,
'H': 1.0,
'N': 1.625,
'O': 1.480,
'P': 1.871,
'S': 1.782,
}

# Bondi 1964,
# https://en.wikipedia.org/wiki/Van_der_Waals_radius
vdW_radii_bondi_1964 = {
'C': 1.7,
'H': 1.09,
'N': 1.55,
'O': 1.52,
'P': 1.8,
'S': 1.8,
}

vdW_radii_dict = {
    'tsai1999': vdW_radii_tsai_1999,
    'bondi1964': vdW_radii_bondi_1964,
    }


# JSON structure parameter keys
JsonParameters = namedtuple('JsonParameters', 'ss fasta resids')
jsonparameters = JsonParameters(
    ss='dssp',
    fasta='fasta',
    resids='resids',
    )

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
    # dssp_ss_keys.bbridge,
    dssp_ss_keys.strand,
    )

dssp_ss_keys.all_loops = (
    dssp_ss_keys.turn,
    dssp_ss_keys.bend,
    dssp_ss_keys.loop,
    dssp_ss_keys.bbridge,  # convention break!
    )

dssp_ss_keys.all = \
    dssp_ss_keys.all_helix \
    + dssp_ss_keys.all_strand \
    + dssp_ss_keys.all_loops

dssp_ss_keys.valid = dssp_ss_keys.all + ('L',)

dssp_trans = str.maketrans(
    ''.join(dssp_ss_keys.all),
    '{}{}{}'.format(
        'H' * len(dssp_ss_keys.all_helix),
        'E' * len(dssp_ss_keys.all_strand),
        'L' * len(dssp_ss_keys.all_loops),
        )
    )


dssp_trans_bytes = bytes.maketrans(
    b''.join(c.encode() for c in dssp_ss_keys.all),
    b'H' * len(dssp_ss_keys.all_helix)
    + b'E' * len(dssp_ss_keys.all_strand)
    + b'L' * len(dssp_ss_keys.all_loops),
    )


# considers solvent and DNA/RNA
# http://www.wwpdb.org/documentation/file-format-content/format33/sect4.html#HET
# _discarded_residues = (
# 'I', 'C', 'G', 'A', 'U', 'I', 'DC', 'DG', 'DA', 'DU', 'DT', 'DI', 'N',
# )
pdb_ligand_codes_file = Path(core_folder, 'chem_comp_parsed.txt')
pdb_lig_codes_manual = Path(core_folder, 'chem_comp_added.txt')
pdb_ligand_codes = set(
    i.strip()
    for i in chain(
        pdb_ligand_codes_file.read_text().split('\n'),
        pdb_lig_codes_manual.read_text().split('\n'),
        )
    if not i.startswith('#')
    )

blocked_ids_file = Path(core_folder, 'discarded_ids.txt')
blocked_ids = [
    i for i in blocked_ids_file.read_text().split('\n')
    if i and not i.startswith('#')
    ]

residue_elements = {'C', 'O', 'N', 'H', 'S', 'Se', 'D'}
minimal_bb_atoms = ['N', 'CA', 'C']  # ordered!

# Builder Definitions  ###
# average values of the backbone angles calculated from
# Dunbrack PISCES
# cull_d200611/200611/cullpdb_pc90_res1.6_R0.25_d200611_chains8807
# float values are represented as ratio of integers
# https://docs.python.org/3/tutorial/floatingpoint.html
average_N_CA_C = 8731046790257777 / 4503599627370496  # +- 0.04375239960584633
average_CA_C_Np1 = 4587708133805365 / 2251799813685248  # +- 0.022904896537130497
average_CA_C_O = 4733796466948169 / 2251799813685248  # +- 0.019050491268134375
average_Np1_C_O = 4825315589323725 / 2251799813685248  # +- 0.017982788310237034
average_Cm1_N_CA = 2385749441983237 / 1125899906842624  # +- 0.029039312259214314

distance_N_CA = 6576479998126497 / 4503599627370496  # 1.46027 +- 0.013036
distance_CA_C = 6861872558247717 / 4503599627370496  # 1.52364 +- 0.012599
distance_C_Np1 = 2996436734567847 / 2251799813685248  # 1.33068 +- 0.009621
distance_C_O = 5556993099130213 / 4503599627370496  # 1.234 +- 0.0121

distance_N_CA_std = 0.013036529567238726
distance_CA_C_std = 0.012599655969373144
distance_C_Np1_std = 0.009621596711934686
