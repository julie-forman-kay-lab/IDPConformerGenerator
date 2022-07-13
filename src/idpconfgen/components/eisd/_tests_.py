"""
Temporary python script used for testing functions.
Will be deleted upon integration completion.
"""
import numpy as np

from pathlib import Path
from Bio.PDB import PDBParser as parser
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
)


def get_scalar(x, y, z):
    return (x**2 + y**2 + z**2)**0.5

pdb = "/home/nemoliu/Documents/eisdTest/asyn_CSSS_nosub/conformer_2_mcsce.pdb"
p = parser()
struct = p.get_structure(id='p', file=pdb)
print(struct[0]['A'][22]['CA'] - struct[0]['A'][23]['CA'])

# jul 12, error, doesn't recognize MC-SCE PDB files
path_pdb = Path(pdb)
s = Structure(path_pdb)
s.build()
s.add_filter(lambda x: x[col_name] == 'CA')
da = s.coords[22, :] - s.coords[20, :]
assert da.shape == (3,)
dist = get_scalar(da[0], da[1], da[2])
print(dist)