import re
import glob
import numpy as np

from functools import partial
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resSeq,
    cols_coords,
    parse_pdb_to_array,
    write_PDB,
    structure_to_pdb,
    )
from idpconfgen.ldrs_helper import (
    align_coords,
    count_clashes,
    next_seeker,
    )

from idpconfgen import Path


# The following function takes elements from `cli_ldrs.py` in IDPConformerGenerator
# Please note that the clash-tolerances are NOT hard-coded in, rather they use the
# default values in IDPConformerGenerator.
def break_idr_attach(idr, fld, upper, lower, output_upper=False, output_lower=False):
    fld_struc = Structure(Path(fld))
    fld_struc.build()
    atom_names = fld_struc.data_array[:, col_name]
    
    fld_seq = fld_struc.data_array[:, col_resSeq].astype(int)
    first_seq = fld_seq[0]

    if first_seq > 1:
        fld_seq -= first_seq - 1
        first_seq = fld_seq[0]
    
    fld_term_idx = {}
    for s, aa in enumerate(fld_seq):
        if aa == lower - 1 and atom_names[s] == 'C':
            fld_term_idx["CL"] = s
        elif aa == lower:
            if atom_names[s] == 'N':
                fld_term_idx["NL"] = s
            elif atom_names[s] == 'CA':
                fld_term_idx["CAL"] = s
        elif aa == upper and atom_names[s] == 'C':
            fld_term_idx['CU'] = s
        elif aa == upper + 1:
            if atom_names[s] == 'N':
                fld_term_idx['NU'] = s
            elif atom_names[s] == 'CA':
                fld_term_idx['CAU'] = s
        if aa > upper + 1:
            break
    
    fld_CLxyz = fld_struc.data_array[fld_term_idx["CL"]][cols_coords].astype(float)
    fld_NLxyz = fld_struc.data_array[fld_term_idx["NL"]][cols_coords].astype(float)
    fld_CALxyz = fld_struc.data_array[fld_term_idx["CAL"]][cols_coords].astype(float)
    fld_coords_L = np.array([fld_CLxyz, fld_NLxyz, fld_CALxyz])
    
    fld_CUxyz = fld_struc.data_array[fld_term_idx["CU"]][cols_coords].astype(float)
    fld_NUxyz = fld_struc.data_array[fld_term_idx["NU"]][cols_coords].astype(float)
    fld_CAUxyz = fld_struc.data_array[fld_term_idx["CAU"]][cols_coords].astype(float)
    fld_coords_U = np.array([fld_CUxyz, fld_NUxyz, fld_CAUxyz])
    
    with open(idr) as f:
        idr_data = f.read()
    try:
        idr_arr = parse_pdb_to_array(idr_data)
    except AssertionError:
        return
    
    if output_lower:
        aligned_break_L = align_coords(idr_arr, fld_coords_L, "C-IDR")
        clashesL, fragmentL = count_clashes(
            aligned_break_L,
            fld_struc,
            "C-IDR",
            max_clash=80,
            tolerance=0.8
            )
        if type(clashesL) is int:
            aligned_struc = structure_to_pdb(fragmentL)
            idr_split = idr.split("/")
            idr_int = int(re.search(r'\d+', idr_split[-1]).group())
            write_PDB(aligned_struc, output_lower + f"aligned_conf_L_{idr_int}.pdb")
    
    elif output_upper:
        aligned_break_U = align_coords(idr_arr, fld_coords_U, "N-IDR")
        clashesU, fragmentU = count_clashes(
            aligned_break_U,
            fld_struc,
            "N-IDR",
            max_clash=80,
            tolerance=0.8
            )
        if type(clashesU) is int:
            aligned_struc = structure_to_pdb(fragmentU)
            idr_split = idr.split("/")
            idr_int = int(re.search(r'\d+', idr_split[-1]).group())
            write_PDB(aligned_struc, output_upper + f"aligned_conf_U_{idr_int}.pdb")


# Feel free to change the variables to other values for other protein
# systems of interest
case = 'slc26a9'
lower = 49
upper = 135
ncores = 32  # MUST BE CHANGED to number of actual CPUs for multiprocessing

# Set-up paths to directories
idr_lst = glob.glob(f".{case}_idr_ANY_bb/*.pdb")
idr_upper = glob.glob(f"./{case}_nterm/*.pdb")
idr_lower = glob.glob(f"./{case}_cterm/*.pdb")
fld = "./7CH1_SLC26A9.pdb"
output_upper = f"./{case}_nterm/"
output_lower = f"./{case}_cterm/"
output_combos = f"./{case}_matches/"

# Execute alignment to N-terminal part of break
execute = partial(
    break_idr_attach,
    fld=fld,
    upper=upper,
    lower=lower,
    output_upper=output_upper,
    output_lower=False,
    )
execute_pool = pool_function(execute, idr_lst, ncores=ncores)
for _ in execute_pool:
    pass

# Execute alignment to C-terminal part of break
execute = partial(
    break_idr_attach,
    fld=fld,
    upper=upper,
    lower=lower,
    output_upper=False,
    output_lower=output_lower,
    )
execute_pool = pool_function(execute, idr_lst, ncores=ncores)
for _ in execute_pool:
    pass

# Execute next-seeker protocol
execute = partial(
    next_seeker,
    nterm_idr_lib=idr_upper,
    max_clash=40,
    tolerance=0.4,
    output_folder=output_combos,
    )
execute_pool = pool_function(execute, idr_lower, ncores=ncores)
for _ in execute_pool:
    pass
