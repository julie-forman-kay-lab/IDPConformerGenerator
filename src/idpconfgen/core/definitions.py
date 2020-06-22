"""Static definitions that serve the whole program infrastructure."""
from argparse import Namespace
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
aa1to3 = {v:k for k, v in aa3to1.items()}


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
    #dssp_ss_keys.bbridge,
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
    b'H' * len(dssp_ss_keys.all_helix) \
    + b'E' * len(dssp_ss_keys.all_strand) \
    + b'L' * len(dssp_ss_keys.all_loops),
        )




# considers solvent and DNA/RNA
# http://www.wwpdb.org/documentation/file-format-content/format33/sect4.html#HET
#_discarded_residues = (
    #'I', 'C', 'G', 'A', 'U', 'I', 'DC', 'DG', 'DA', 'DU', 'DT', 'DI', 'N',
    #)
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
