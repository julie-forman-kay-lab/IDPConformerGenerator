"""
Temporary python script used for testing functions.
Will be deleted upon integration completion.
"""

from Bio.PDB import PDBParser as parser
from idpconfgen.libs.libstructure import Structure

pdb = "/home/nemoliu/Documents/eisdTest/asyn_CSSS_nosub/conformer_2_mcsce.pdb"
p = parser()
struct = p.get_structure(id='p', file=pdb)
print(struct[0]['A'][1]['CA'] - struct[0]['A'][2]['CA'])

# jul 12, error, doesn't recognize MC-SCE PDB files
s = Structure(pdb)
s.build()
print(s[0]['A'][1]['CA'] - s[0]['A'][2]['CA'])