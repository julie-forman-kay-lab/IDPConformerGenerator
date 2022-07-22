"""
Temporary python script used for testing functions.
Will be deleted upon integration completion.
"""
import numpy as np
import pandas as pd

from pathlib import Path
from Bio.PDB import PDBParser as parser
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    )

atoms = []
values = []
errors = []

lines = [
        "_Atom_chem_shift.ID",
        "_Atom_chem_shift.Assembly_atom_ID",
        "_Atom_chem_shift.Entity_assembly_ID",
        "_Atom_chem_shift.Entity_assembly_asym_ID",
        "_Atom_chem_shift.Entity_ID",
        "_Atom_chem_shift.Comp_index_ID",
        "_Atom_chem_shift.Seq_ID",
        "_Atom_chem_shift.Comp_ID",
        "_Atom_chem_shift.Atom_ID",
        "_Atom_chem_shift.Atom_type",
        "_Atom_chem_shift.Atom_isotope_number",
        "_Atom_chem_shift.Val",
        "_Atom_chem_shift.Val_err",
        "1     .   2   .   1   2    2    MET   C    C   13   176.341   .",
        "2     .   2   .   1   2    2    MET   CA   C   13   55.634    .",
        "3     .   2   .   1   2    2    MET   CB   C   13   32.613    .",
        "4     .   2   .   1   3    3    GLU   H    H   1    8.781     .",
        "5     .   2   .   1   3    3    GLU   C    C   13   176.158   .",
        "6     .   2   .   1   3    3    GLU   CA   C   13   56.8      .",
        "7     .   2   .   1   3    3    GLU   CB   C   13   29.905    .",
        "8     .   2   .   1   3    3    GLU   N    N   15   122.712   .",
        "9     .   2   .   1   4    4    ALA   H    H   1    8.419     .",
        "10    .   2   .   1   4    4    ALA   C    C   13   177.602   .",
        "11    .   2   .   1   4    4    ALA   CA   C   13   52.395    .",
    ]
    
for i, line in enumerate(lines):
    if "_" in line:
        dtype = line.split(".")[1]
        if dtype == "Val": data_idx = i
        elif dtype == "Val_err": error_idx = i
        elif "Atom_ID" in dtype: atom_idx = i
    else:
        start_idx = i
        break
    
for idx in range(start_idx, len(lines)):
    splitted = lines[idx].split()
    
    atoms.append(splitted[atom_idx])
    values.append(splitted[data_idx])
    errors.append(splitted[error_idx])
    
data = pd.DataFrame({"atom": atoms, "values": values, "errors": errors})
print(data_idx, error_idx, atom_idx, start_idx)
print(data)



'''
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
'''
