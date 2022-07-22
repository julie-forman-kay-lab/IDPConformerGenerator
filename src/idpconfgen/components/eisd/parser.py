"""
Main logic for parsing experimental and back-calculated raw outputs.
Accepts primarily NMR-STAR formatted files for solution data.

USAGE EXAMPLE:
    For NMR-STAR files, it must be the data between the `loop_` and 
    `stop_` codes without new lines. The easiest way is to just
    copy the block of data from NMR-STAR file from the BMRB
    to a new text file.
    
    For example an exerpt from entry 4155:
    _Atom_chem_shift.ID
    _Atom_chem_shift.Assembly_atom_ID
    _Atom_chem_shift.Entity_assembly_ID
    _Atom_chem_shift.Entity_assembly_asym_ID
    _Atom_chem_shift.Entity_ID
    _Atom_chem_shift.Comp_index_ID
    _Atom_chem_shift.Seq_ID
    _Atom_chem_shift.Comp_ID
    _Atom_chem_shift.Atom_ID
    _Atom_chem_shift.Atom_type
    _Atom_chem_shift.Atom_isotope_number
    _Atom_chem_shift.Val
    _Atom_chem_shift.Val_err
    1     .   2   .   1   2    2    MET   C    C   13   176.341   .
    2     .   2   .   1   2    2    MET   CA   C   13   55.634    .
    3     .   2   .   1   2    2    MET   CB   C   13   32.613    .
    4     .   2   .   1   3    3    GLU   H    H   1    8.781     .
    5     .   2   .   1   3    3    GLU   C    C   13   176.158   .
    6     .   2   .   1   3    3    GLU   CA   C   13   56.8      .
    7     .   2   .   1   3    3    GLU   CB   C   13   29.905    .
    8     .   2   .   1   3    3    GLU   N    N   15   122.712   .
    9     .   2   .   1   4    4    ALA   H    H   1    8.419     .
    10    .   2   .   1   4    4    ALA   C    C   13   177.602   .
    11    .   2   .   1   4    4    ALA   CA   C   13   52.395    .

Currently accepting the following formats:
    SAXS: SASBDB, CRYSOL
    Chemical shift: NMR-STAR
    FRET: (ask greg/claudiu for the "standard")
    J-Couplings: NMR-STAR (coupling constants)
    NOE: NMR-STAR (homonuclear - distances NMR-STAR)
    PRE: (changes in R2 or as distance (use NOE) - NMR STAR)
    RDC: NMR-STAR
    Rh: (single values)
"""
import numpy as np
import pandas as pd

from idpconfgen.components.eisd import (
    star_val,
    star_err,
    star_atmID,
    star_max,
    star_min,
    parse_mode_exp,
    saxs_name,
    cs_name,
    fret_name,
    jc_name,
    noe_name,
    pre_name,
    rdc_name,
    rh_name,
    )

class Stack():
    def __init__(self, name, data, sigma=None, mu=None):
        """
        Initialize stack of experimental or back-calculated
        data.

        Parameters
        ----------
        name : string
            Name of experimental datatype.
        
        pd.DataFrame
            Containing atoms, values, and errors for the 3 columns
        """
        self.name = name
        self.data = data
        self.sigma = sigma
        self.mu = mu


def parse_cs_data(fpath):
    atoms = []
    values = []
    errors = []
    
    with open(fpath) as f:
        lines = f.readlines()
        
        for i, line in enumerate(lines):
            if "_" in line:
                dtype = line.split(".")[1].strip()
                if dtype == "Val": data_idx = i
                elif dtype == "Val_err": error_idx = i
                elif dtype == "Atom_ID": atom_idx = i
            else:
                start_idx = i
                break
            
        for idx in range(start_idx, len(lines)):
            splitted = lines[idx].split()
            
            atoms.append(splitted[atom_idx])
            values.append(float(splitted[data_idx]))
            errors.append(float(splitted[error_idx]))
            
    return pd.DataFrame({star_atmID: atoms, star_val: values, star_err: errors})


# Good for NMR-STAR formatted JC, RDC, NOE, PRE
def parse_nmrstar_data(fpath, type=None):
    values = []
    upper = []
    lower = []
    errors = []
    
    with open(fpath) as f:
        lines = f.readlines()
        
        for i, line in enumerate(lines):
            if "_" in line:
                dtype = line.split(".")[1].strip()
                if dtype == "Val": data_idx = i
                elif dtype == "Val_err": error_idx = i
                if type == noe_name or type == pre_name:
                    if dtype == "Val_max": max_idx = i
                    elif dtype == "Val_min": min_idx = i
            else:
                start_idx = i
                break
            
        for idx in range(start_idx, len(lines)):
            splitted = lines[idx].split()
            
            values.append(float(splitted[data_idx]))
            errors.append(float(splitted[error_idx]))
            
            if type == noe_name or type == pre_name:
                max = splitted[max_idx]
                min = splitted[min_idx]
                if max == ".": max = 0.0
                if min == ".": min = 0.0
                upper.append(float(max))
                lower.append(float(min))
    
    if type == noe_name or type == pre_name:
        return pd.DataFrame({star_val: values, star_max: upper, star_min: lower, star_err: errors})
    
    return pd.DataFrame({star_val: values, star_err: errors})

def parse_data(filenames, mode, bc_errors={}):
    """
    Main function to read/parse all experimental and back-calculated files.

    Parameters
    ----------
    filenames : dict
        Dictionary of modules with their relative path
        to the data file. First key-layer has all of the
        different experimental modules eisd can handel.
    
    mode : str
        Parameter must be one of the following:
        - 'exp' = experimental data
        - 'bc' = back-calculated data
        
    Returns
    -------
    parsed : dict
        Dictionary of properties with their pandas dataframe.
    
    errlogs : list
        List of possible messages to display for the user.    
    """
    parsed = {}
    errlogs = []
    data = pd.DataFrame(dtype=np.float)
    
    # Parse experimental data
    if mode == parse_mode_exp:
        for module in filenames:
            
            if module == saxs_name:
                # do saxs parsing
                None
            elif module == cs_name:
                # do cs parsing
                data = parse_cs_data(filenames[module])
            elif module == fret_name:
                # do fret parsing
                None
            elif module == jc_name:
                # do JC parsing
                data = parse_nmrstar_data(filenames[module])
            elif module == noe_name:
                # do NOE parsing
                data = parse_nmrstar_data(filenames[module], noe_name)
            elif module == pre_name:
                # no PRE parsing
                data = parse_nmrstar_data(filenames[module], pre_name)
            elif module == rdc_name:
                # do RDC parsing
                data = parse_nmrstar_data(filenames[module])
            elif module == rh_name:
                # do Rh parsing
                None
                     
            parsed[module] = Stack(module, data, None, None)
    
    # Parse back calculated data
    else:
        None

    return parsed, errlogs