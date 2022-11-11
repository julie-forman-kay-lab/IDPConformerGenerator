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
    - Case 3: generate disordered loops from tether regions and try to connect them

Important note: do not perform clash check within folded region
"""
# Definitions of cases (subject to change)
n_idr = 1
c_idr = 2
break_idr = 3

from functools import partial

import numpy as np

from idpconfgen.core.definitions import aa3to1
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


def consecutive_grouper(seq):
    """
    Use negative indexing to group together consecutive numbers.

    References:
    https://stackoverflow.com/questions/70363072/group-together-consecutive-numbers-in-a-list
    
    Parameters
    ----------
    seq : string
        Special sequence where indices of disordered residues are stored.
    
    Return
    ------
    bounds : list
        List of ranges for boundaries of disordered sequences.
    """
    grouped = [[seq[0]]]
    for x in seq[1:]:
        if x == grouped[-1][-1] + 1:
            grouped[-1].append(x)
        else:
            grouped.append([x])
            
    bounds=[]
    for group in grouped:
        first = group[0]
        last = group[len(group) - 1]
        bounds.append((first, last))
    
    return bounds


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
    fld_seqs : list
        List of FASTA sequence of folded regions in the sequence of fdata.
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
        whole = []
        for i, dist in enumerate(coords_distances):
            if dist < 2.1:
                whole.append(i)
        
        whole = consecutive_grouper(whole)
        fld_seqs = []
        for idx in whole:
            fld_idx = list(range(idx[0], idx[1], 3))
            fld_seqs.append(''.join(aa3to1.get(f) for f in data[:, col_resName][fld_idx].tolist()))
        
        return fld_seqs
    
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
    