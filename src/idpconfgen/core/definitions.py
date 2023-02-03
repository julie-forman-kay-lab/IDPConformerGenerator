"""Static definitions that serve the whole program infrastructure."""
from argparse import Namespace
from collections import namedtuple
from itertools import chain
# does not import the Path from IDPConfgen to avoid circular imports
from pathlib import Path


core_folder = Path(__file__).parent


# Bond Geometry definitions
# Keys in library:
bgeo_Cm1NCa = 'Cm1_N_Ca'
bgeo_NCaC = 'N_Ca_C'
bgeo_CaCNp1 = 'Ca_C_Np1'
bgeo_CaCO = 'Ca_C_O'

bgeo_NCa = 'N_Ca'
bgeo_CaC = 'Ca_C'
bgeo_CNp1 = 'C_Np1'
bgeo_CO = 'C_O'


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
    'HIP': 'p',  # other name for double protonated histidine
    'HIE': 'e',  # epsilon protonated histidine
    'HID': 'd',  # gamma protonated histidine
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'SEP': 's',  # phosphorylated serine
    'THR': 'T',
    'TPO': 't',  # phosphorylated threonine
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'}

# Amino-acid 1 to 3 letter code dictionary
aa1to3 = {v: k for k, v in aa3to1.items()}

aa1set = set(aa1to3.keys())
aa3set = set(aa3to1.keys())

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
    # adding PPII definition
    polypII="P",
    )

dssp_ss_keys.all_helix = (
    dssp_ss_keys.ahelix,
    dssp_ss_keys.helix_3,
    # dssp_ss_keys.helix_5,
    )

# dssp_ss_keys.helix_3 = (dssp_ss_keys.helix_3,)
# dssp_ss_keys.helix_5 = (dssp_ss_keys.helix_5,)

dssp_ss_keys.all_strand = (
    # dssp_ss_keys.bbridge,
    dssp_ss_keys.strand,
    )

dssp_ss_keys.all_loops = (
    dssp_ss_keys.turn,
    dssp_ss_keys.bend,
    dssp_ss_keys.loop,
    dssp_ss_keys.bbridge,  # convention break!
    # helix_5 added as loops, following
    # Balasco, N. et al. BioMed Research International vol. 2017 e2617629 (2017)
    dssp_ss_keys.helix_5,
    # polyproline II helix added as loops, following
    # Mansiaux, Y., Joseph, A. P., Gelly, J.-C. & Brevern, A. G.
    # PLOS ONE 6, e18401 (2011)
    dssp_ss_keys.polypII,
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
        'L' * len(dssp_ss_keys.all_loops)
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

#  Builder Definitions  ###
#  average values of the backbone angles calculated from
#  Dunbrack PISCES
#  cull_d200611/200611/cullpdb_pc90_res1.6_R0.25_d200611_chains8807
#  float values are represented as ratio of integers
#  https://docs.python.org/3/tutorial/floatingpoint.html
# average_N_CA_C = 8731046790257777 / 4503599627370496  # +- 0.04375239960584633
# average_CA_C_Np1 = 4587708133805365 / 2251799813685248  # +- 0.022904896537130
# average_Np1_C_O = 4733796466948169 / 2251799813685248  # +- 0.0190504912681343
# average_CA_C_O = 4825315589323725 / 2251799813685248  # +- 0.01798278831023703
# average_Cm1_N_CA = 2385749441983237 / 1125899906842624  # +- 0.029039312259214
# bend_CA_C_OXT = 2 * pi / 3
#
# # pi corrected angles needed for the building algorithm
# build_bend_CA_C_Np1 = (pi - average_CA_C_Np1) / 2
# build_bend_Cm1_N_CA = (pi - average_Cm1_N_CA) / 2
# build_bend_N_CA_C = (pi - average_N_CA_C) / 2
# build_bend_CA_C_OXT = (pi - bend_CA_C_OXT) / 2
# build_bend_CA_C_O = average_CA_C_O / 2  # this angle does not require `pi -`
#
# distance_N_CA = 6576479998126497 / 4503599627370496  # 1.46027 +- 0.013036
# distance_CA_C = 6861872558247717 / 4503599627370496  # 1.52364 +- 0.012599
# distance_C_Np1 = 2996436734567847 / 2251799813685248  # 1.33068 +- 0.009621
# distance_C_O = 5556993099130213 / 4503599627370496  # 1.234 +- 0.0121
# distance_C_OXT = 1.27
#
# distance_N_CA_std = 0.013036529567238726
# distance_CA_C_std = 0.012599655969373144
# distance_C_Np1_std = 0.009621596711934686
