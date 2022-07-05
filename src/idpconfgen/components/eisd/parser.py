"""
Main logic for parsing experimental and back-calculated raw outputs.

Currently accepting the following formats:
    SAXS: SASBDB, CRYSOL
    Chemical shift: NMR-STAR, SHIFTX2
    FRET: (ask greg/claudiu for the "standard")
    J-Couplings:
    NOE: NMR-STAR
    PRE:
    RDC: NMR-STAR
    Rh: (single value)
"""
import numpy as np
import pandas as pd

from idpconfgen.components.eisd import (
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
        self.name = name
        self.data = data
        self.sigma = sigma
        self.mu = mu
    
    def get_index(self):
        pass

def parse_data(filenames, mode):
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
    
    # Parse experimental data
    if mode == parse_mode_exp:
        for module in filenames:
            p = pd.read_csv(filenames[module])
            
            if module == saxs_name:
                # do saxs parsing
                None
            elif module == cs_name:
                # do cs parsing
                None
            elif module == fret_name:
                # do fret parsing
                None
            elif module == jc_name:
                # do JC parsing
                None
            elif module == noe_name:
                # do NOE parsing
                None
            elif module == pre_name:
                # no PRE parsing
                None
            elif module == rdc_name:
                # do RDC parsing
                None
            elif module == rh_name:
                # do Rh parsing
                None
                     
            parsed[module] = Stack(module, p, None, None)
    # Parse back calculated data
    else:
        None

    return parsed, errlogs