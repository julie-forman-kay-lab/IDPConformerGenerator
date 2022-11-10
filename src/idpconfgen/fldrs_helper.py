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

from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen import Path, log
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resName,
    col_resSeq,
    cols_coords,
)
from idpconfgen.logger import S, T, init_files, report_on_crash


def disorder_bounds(seq):
    """
    Find out where the disordered residues are given a special sequence.

    Parameters
    ----------
    seq : string
        Special sequence where folded residues are marked as 'X'
    
    Return
    ------
    dis_bounds : list
        List of ranges for boundaries of disordered sequences
    """


def break_check(fdata):
    """
    Calculate where breaks are in the backbone.
    
    Inspired from `get_torsions` in `libhigherlevel`.

    Parameters
    ----------
    fdata : str, bytes or Path
        A path to the structure file, or the string representing
        the file.
        In fact, accepts any type `:class:libstructure.Structure` would
        accept.
    
    Return
    ------
    breaks : list
        List of ranges of breaks in the sequence of fdata.
    """
    structure = Structure(fdata)
    structure.build()
    structure.add_filter_backbone(minimal=True)

    data = structure.filtered_atoms
    names = data[:, col_name]
    coords_raw = structure.coords

    n_mask = names == 'N'
    ca_mask = names == 'CA'
    c_mask = names == 'C'

    n = coords_raw[n_mask, :]
    ca = coords_raw[ca_mask, :]
    c = coords_raw[c_mask, :]

    try:
        coords = np.empty((n.shape[0] * 3, 3), dtype=np.float64)
        coords[0::3, :] = n
        coords[1::3, :] = ca
        coords[2::3, :] = c
    except ValueError as err:
        errmsg = (
            'Coordinates do not match expectation. '
            'Some possibly missing.'
            )
        raise IDPConfGenException(errmsg) from err

    coords_distances = np.linalg.norm(coords[:-1, :] - coords[1:, :], axis=1)
    assert coords_distances.size == coords.shape[0] - 1
    
    if np.any(coords_distances > 2.1):
        breaks = []
        for i, dist in enumerate(coords_distances):
            if dist > 2.1:
                # The following normalizes residue numbers in case PDB doesn't start with 1
                # then returns the index of the ranges of disordered regions
                first = int(data[:, col_resSeq][0])
                current = int(data[:, col_resSeq][i])
                next = int(data[:, col_resSeq][i + 1])
                # we add 1 because "first" is inclusive
                breaks.append((current - first + 1, next - first))
        return breaks
    
    return False
        

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
    