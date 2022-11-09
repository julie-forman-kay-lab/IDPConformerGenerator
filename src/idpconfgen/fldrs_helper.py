"""
Client for building IDRs on PDB files in the cartesian coordinate space.

Methodology deviates from traditional IDP or beads-on-a-string FLDR/S approach.

Name: FLDR/S (Folded region/structure sampling)

Plan of action:
1. Create internal database of start/end point teathers
    - Case 1: N-IDR - build backwards from N-terminal break point
    - Case 2: C-IDR - build forwards by populating coordinate array
    with PDB structure
    - Case 3: Break-IDR - construct internal dictionary of starting
    and end-point teathers
2. IDP generation will use original philosophies to retain sequence/secondary
    structure identities in context of protein of interest.
    - Case 1: append empty cells to array and build backwards by inverting torsions
    - Case 2: populate coordinate array then build CIDR as usual
    - Case 3: apply loop-closure algorithms
    (inverse kinematics: CCD, RCD, KIC, RLG, etc.)

Important note: do not perform clash check within folded region
"""
# Definitions of cases (subject to change)
n_idr = 1
c_idr = 2
break_idr = 3

from functools import partial

import numpy as np

from idpconfgen import Path, log
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resName,
    col_resSeq,
    cols_coords,
)
from idpconfgen.logger import S, T, init_files, report_on_crash


def psurgeon(case, fld_struc, idp_coords):
    """
    Protein surgeon grafts disordered regions onto folded structures.

    Parameters
    ----------
    case : int
        Case number 1 = N-IDR, 2 = C-IDR, 3 = break-IDR
    
    fld_struc : Structure
        IDPConformerGenerator Structure class of folded structure
    
    idp_coords : np.ndarray
        Array of coordinates for IDP generated
    
    Returns
    -------
    pdb_coords : np.ndarray
        Array of final coordinates 
    """
    if case == n_idr:
        None
    elif case == c_idr:
        None
    elif case == break_idr:
        None
    else:
        # code should not reach here
        return
    