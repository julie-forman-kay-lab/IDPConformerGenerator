# Note to self and devs: calculators will be done on a separate repository
# to maintain modularity and flexibility, this script will be removed soon
"""
Contains functions to back-calculate experimental
data from generated conformer pools.

Partially inspired from: https://github.com/Oufan75/X-EISD/blob/master/EISD_back_calc.ipynb
"""
import numpy as np
import pandas as pd

from idpconfgen.libs.libstructure import Structure, col_name
from idpconfgen.libs.libcalc import calc_torsion_angles

# chemical shift back calculations using UCBShifts (Li et al. 2020
# https://github.com/THGLab/CSpred.git); reading outputs
UCBSHIFT_MAP = {'H':0, 'HA':1, 'C':2, 'CA':3, 'CB':4, 'N':5}


def get_scalar(x, y, z):
    return (x ** 2 + y ** 2 + z ** 2) ** 0.5


def fret_calculator(exp_file, struct_path):
    fret_bc = []
    
    exp = pd.read_csv(exp_file)
    res1 = exp.res1.values.astype(np.int)
    res2 = exp.res2.values.astype(np.int)
    assert abs(res1 - res2) == len(s.residues) - 1
    
    scaler = exp.scale.values
    
    s = Structure(struct_path)
    s.build()
    s.add_filter(lambda x: x[col_name] == 'CA')
    #assumes CA as atom labeled
    for j in range(exp.shape[0]):
        r1 = np.int(res1[j])
        r2 = np.int(res2[j])
        # distance between 2 CA atoms
        _d = s.coords[r1, :] - s.coords[r2, :]
        assert _d.shape == (3,)
        d = get_scalar(_d[0], _d[1], _d[2])
        
        # scale_factor to adjust for dye size and CA to label distances
        scale_factor = ((np.abs(r1 - r2) + 7) / np.abs(r1 - r2)) ** 0.5
        d = d * scale_factor
        eff = 1.0 / (1.0 + (d / scaler[j]) ** 6.0)

        fret_bc.append(eff)
    
    fret_bc = np.reshape(fret_bc, (-1, exp.shape[0]))
    return fret_bc
