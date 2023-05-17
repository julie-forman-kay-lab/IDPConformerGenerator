import glob

from idpconfgen.libs.libstructure import structure_to_pdb, write_PDB 
from idpconfgen.ldrs_helper import psurgeon

from idpconfgen import Path

# Feel free to change the variables to other values for other protein
# systems of interest
case = 'slc26a9'
lower = 49
upper = 135

fld = "./7CH1_SLC26A9.pdb"
combos_list = glob.glob(f"./{case}_sidechains/*.pdb")
output_stitched = f"./{case}_results/"

# Execute stitching protocol
for i, idr in enumerate(combos_list):
    combined_struc = psurgeon([[Path(idr)]], Path(fld), ["Break-IDR"], [(lower, upper)])
    combined_pdb = structure_to_pdb(combined_struc)
    write_PDB(combined_pdb, output_stitched + f"conformer_{i + 1}.pdb")
