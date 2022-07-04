"""
Main logic for parsing experimental and back-calculated raw outputs.

Currently accepting the following formats:
    SAXS:
    Chemical shift: NMR_STAR, SHIFTX2
    FRET:
    J-Couplings:
    NOE:
    PRE:
    RDC:
    Rh:
"""
import numpy as np
import pandas as pd

from idpconfgen.components.eisd import (
    parse_mode_exp,
    parse_mode_back,
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
    
    if mode == parse_mode_exp:
        try:
            saxs = pd.read_csv(filenames[saxs_name])
            
            parsed[saxs_name] = Stack(saxs_name, saxs, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {saxs_name} experimental data not found."
            "Skipping this module"
            )
        
        try:
            cs = pd.read_csv(filenames[cs_name])
            
            parsed[cs_name] = Stack(cs_name, cs, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {cs_name} experimental data not found."
            "Skipping this module"
            )
        
        try:
            fret = pd.read_csv(filenames[fret_name])
            
            parsed[fret_name] = Stack(fret_name, fret, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {fret_name} experimental data not found."
            "Skipping this module"
            )
        
        try:
            jc = pd.read_csv(filenames[jc_name])
            
            parsed[jc_name] = Stack(jc_name, jc, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {jc_name} experimental data not found."
            "Skipping this module"
            )
        
        try:
            noe = pd.read_csv(filenames[noe_name])
            
            parsed[noe_name] = Stack(noe_name, noe, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {noe_name} experimental data not found."
            "Skipping this module"
            )
        try:
            pre = pd.read_csv(filenames[pre_name])
            
            parsed[pre_name] = Stack(pre_name, pre, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {pre_name} experimental data not found."
            "Skipping this module"
            )
        try:
            rdc = pd.read_csv(filenames[rdc_name])
            
            parsed[rdc_name] = Stack(rdc_name, rdc, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {rdc_name} experimental data not found."
            "Skipping this module"
            )
        
        try:
            rh = pd.read_csv(filenames[rh_name])
            
            parsed[rh_name] = Stack(rh_name, rh, None, None)
        except FileNotFoundError:
            errlogs.append
            (f"Note: {rh_name} experimental data not found."
            "Skipping this module"
            )
            

    return parsed, errlogs