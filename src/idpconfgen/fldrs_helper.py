"""
Client for building IDRs on PDB files in the cartesian coordinate space.

Methodology deviates from traditional IDP or beads-on-a-string FLDR/S approach.

Name: FLDR/S (Folded disordered region/structure sampling)

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

disorder_cases = {
    0: "N-IDR",
    1: "Break-IDR",
    2: "C-IDR",
    }

import random

import numpy as np

from idpconfgen import Path

from idpconfgen.core.definitions import aa3to1, vdW_radii_tsai_1999
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resName,
    col_resSeq,
    cols_coords,
    col_x,
    col_y,
    col_z,
)


def calculate_distance(coords1, coords2):
    """
    Returns the distance between two 3D coordinates.
    
    Calculates the distance between 2 coordinates using Euclidean distance
    formula.
    
    Parameters
    ----------
    coords1 : np.ndarray
    
    coords2 : np.ndarray
    
    Return
    ------
    float distance
    """
    return np.sqrt(np.sum((coords1 - coords2) ** 2))


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
        bounds.append((first, last + 1))
    
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
    
    return


def pmover(case, fld_xyz, idp_path):
    """
    Protein cartesian space mover.
    
    Shifts entire protein chain based on one point.

    Parameters
    ----------
    case : string
        Case could be `nidr`, `cidr`, or `break` as defined above.
    
    fld_xyz : tuple
        Backbone N(x, y, z) float coordinates of interest
        where we want to move the IDP chain relative to.
    
    idp_path : Path
        Path to the IDP conformer we want to move.
    
    Returns
    -------
    Overwrites PDB of IDP conformer with new coordinates.
    """
    Nx = fld_xyz[0]
    Ny = fld_xyz[1]
    Nz = fld_xyz[2]
    
    structure = Structure(idp_path)
    structure.build()
    atom_names = structure.data_array[:, col_name]
    
    if case == disorder_cases[0]:  # N-IDR
        # In the case of N-IDR, we want to move relative to C-term
        # A little bit complicated, need to calculate difference between
        # C-term Nitrogen on IDP and N(x,y,z)
        for i, atom in enumerate(atom_names):
            if atom_names[len(atom_names) - 1 - i] == "N":
                index = len(atom_names) - 1 - i
                break
    elif case == disorder_cases[1]:  # break
        pass
    elif case == disorder_cases[2]:  # C-IDR
        # In the case of C-IDR, we want to move relative to carbonyl C
        for i, atom in enumerate(atom_names):
            if atom == "C":
                index = i
                break
    
    idp_xyz = structure.data_array[index][cols_coords]
    dx = Nx - float(idp_xyz[0])
    dy = Ny - float(idp_xyz[1])
    dz = Nz - float(idp_xyz[2])
    
    for i, coords in enumerate(structure.data_array[:, cols_coords]):
        x = str(round(dx + float(coords[0]), 3))
        y = str(round(dy + float(coords[1]), 3))
        z = str(round(dz + float(coords[2]), 3))
        
        structure.data_array[i][col_x] = x
        structure.data_array[i][col_y] = y
        structure.data_array[i][col_z] = z
    
    structure.write_PDB(idp_path)
    
    return


def rotator(chain, case):
    """
    Rotation function that rotates the protein chain randomly.
    
    The point of rotation depends on the case of IDR.
    For example, with N-IDR, the point of rotation is about
    the CA atom at the C-term of the IDR chain.
    
    Parameters
    ----------
    chain : Path
        Path to the chain of interest we want to rotate.
    
    case : string
        Disordered case of interest determines which point
        to rotate about.

    Returns
    -------
    rotated : IDPConformerGenerator.Structure
        Structure object with the new rotated coordinates.
    """
    minrad = 0
    maxrad = 2 * np.pi
    # Select random angle to rotate
    angle = random.uniform(minrad, maxrad)
    
    idp = Structure(chain)
    idp.build()
    atom_names = idp.data_array[:, col_name]
    
    if case == disorder_cases[0]:  # N-IDR
        for i, atom in enumerate(atom_names):
            if atom_names[len(atom_names) - 1 - i] == "N":
                index = len(atom_names) - 1 - i
                break
    elif case == disorder_cases[2]:  # C-IDR
        for i, atom in enumerate(atom_names):
            if atom == "C":
                index = i
                break
    
    point_xyz = idp.data_array[index][cols_coords]
    px = float(point_xyz[0])
    py = float(point_xyz[1])
    pz = float(point_xyz[2])
    
    coords = idp.data_array[:, cols_coords]  # needs confirming
    
    # Create a 3x3 rotation matrix
    rotation_matrix = np.array(
        [[np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0, 0, 1]]
    )
    
    # Create an array of the coordinates to be rotated
    # ones appended to the end to make a 4xN matrix for multiplication
    coords_with_ones = np.hstack([coords, np.ones((coords.shape[0], 1))])
    
    rotated_coords = np.dot(rotation_matrix, coords_with_ones.T).T
    rotated_coords = rotated_coords[:, :-1]
    
    # Translate rotated coordinates so they are rotated about the desired point
    rotated_coords[:, 0] += px - np.dot(rotation_matrix[0, :2], [px, py])
    rotated_coords[:, 1] += py - np.dot(rotation_matrix[1, :2], [px, py])
    rotated_coords[:, 2] += pz - np.dot(rotation_matrix[2, :2], [px, py])

    idp.data_array[:, cols_coords] = rotated_coords
    
    return idp


def count_clashes(chain1, chain2):
    """
    Checks for steric clashes between two protein chains.
    
    Input should be a list of tuples (each has atom symbol and 3D coordinates)
    Function returns number of steric clashes detected between two chains.

    Parameters
    ----------
    chain1 : List[Tuple[str, np.ndarray]]
    
    chain2 : List[Tuple[str, np.ndarray]]
    
    Returns
    -------
    num_clashes : int
        Number of steric clashes determined using vdW radii
    """
    num_clashes = 0
    
    # Loop through all pairs of atoms in the 2 protein chains
    for atom1, coords1 in chain1:
        for atom2, coords2 in chain2:
            # calculate distance between atoms
            distance = calculate_distance(coords1, coords2)
            
            # get vdW radii for each atom
            vdw_radius1 = vdW_radii_tsai_1999[atom1]
            vdw_radius2 = vdW_radii_tsai_1999[atom2]
            
            # Check if a steric clash is detected by comparing
            # distance between atoms to the sum of their vdW radii
            if distance < vdw_radius1 + vdw_radius2:
                num_clashes += 1
    
    return num_clashes
    


def psurgeon(idp_struc, case, fl_seq, bounds, fld_struc):
    """
    Protein surgeon grafts disordered regions onto folded structures.

    Parameters
    ----------
    case : string
        Case could be `nidr`, `cidr`, or `break` as defined above.
    
    fl_seq : string
        Full length sequence of requested protein to build.
    
    bounds : tuple
        Starting and ending indicies of where disordered region is.
    
    fld_struc : Path
        Path to the folded structure .PDB file.
    
    idp_struc : Path
        Path to the donor conformer to graft.
    
    Returns
    -------
    pdb_coords : np.ndarray
        Array of final coordinates 
    """
    fld = Structure(fld_struc)
    idr = Structure(idp_struc)
    fld.build()
    idr.build()
    
    # For break and cidr cases, need to check if grafting
    # was already completed for N-IDR and previous breaks
    if case == disorder_cases[0]:
        """
        Plan (Nov 18):
        1. Remove last residue (C-term) on idr
        2. Initialize a new structure array with total length
        of idr and fld after C-term removal
        3. Populate new structure array with idr first then fld
        """
            
        pass
    elif case == disorder_cases[1]:  # break
        pass
    elif case == disorder_cases[2]:  # C-IDR
        pass
    
    return
