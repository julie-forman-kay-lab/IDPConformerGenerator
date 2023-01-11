"""
Client for building IDRs on PDB files in the cartesian coordinate space.

Methodology deviates from traditional IDP or beads-on-a-string FLDR/S approach.

Name: FLDR/S (Folded disordered region/structure sampling)
"""

disorder_cases = {
    0: "N-IDR",
    1: "Break-IDR",
    2: "C-IDR",
}

import random
import os

import numpy as np

from itertools import product
from idpconfgen import Path

from idpconfgen.core.definitions import aa3to1, vdW_radii_tsai_1999
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs.libstructure import (
    Structure,
    col_chainID,
    col_element,
    col_name,
    col_resName,
    col_resSeq,
    col_serial,
    col_segid,
    cols_coords,
    col_x,
    col_y,
    col_z,
)


def tolerance_calculator(tolerance):
    """
    Calculates the max number of spherical clashes and distance tolerance
    given a floating point tolerance between 0-1.

    Parameter
    ---------
    tolerance : float
    
    Returns
    -------
    max_clash : int
    dist_tolerance : float
    """
    if tolerance > 1.0:
        tolerance = 1.0
    elif tolerance < 0.0:
        tolerance = 0.0
        
    max_clash = int(tolerance * 100)
    dist_tolerance = tolerance
    
    return max_clash, dist_tolerance


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


def store_idp_paths(folder, temp_dir):
    """
    Stores all of the paths for different cases of IDRs in a dictionary.

    Parameters
    ----------
    folder : Path
        Output folder of interest
    
    temp_dir : str
        Name of the temporary directory of interest
    
    Returns
    -------
    disorder_case : dict
        Dictionary of paths to conformers associated with each disorder case
    """
    case_path = {}
    
    if os.path.exists(folder.joinpath(temp_dir + disorder_cases[0])):
        fpath = folder.joinpath(temp_dir + disorder_cases[0])
        idr_confs = os.listdir(fpath)
        case_path[disorder_cases[0]] = [Path(fpath.joinpath(cpath)) for cpath in idr_confs]
    if os.path.exists(folder.joinpath(temp_dir + disorder_cases[1])):
        fpath = folder.joinpath(temp_dir + disorder_cases[1])
        idr_confs = os.listdir(fpath)
        case_path[disorder_cases[1]] = [Path(fpath.joinpath(cpath)) for cpath in idr_confs]
        # What to do if we have multiple breaks? Maybe split to subdirs
    if os.path.exists(folder.joinpath(temp_dir + disorder_cases[2])):
        fpath = folder.joinpath(temp_dir + disorder_cases[2])
        idr_confs = os.listdir(fpath)
        case_path[disorder_cases[2]] = [Path(fpath.joinpath(cpath)) for cpath in idr_confs]

    return case_path


def create_combinations(list1, list2, num_combinations):
    """
    Create unique combinations between two lists.
    
    Made for N-IDR and C-IDR combinations. Where list1 = N-IDR paths,
    and list2 = C-IDR paths. Itertools product is used here because
    order matters.
    
    Parameters
    ----------
    list1 : list
    
    list2 : list
    
    num_combinations : int
    
    Return
    ------
    selected_combinations : list
        List of tuples of different combinations as follows:
        [(item from list 1, item from list 2), ...]
    """
    all_combinations = list(product(list1, list2))
    max_combinations = len(all_combinations)

    selected_combinations = random.sample(all_combinations, min(num_combinations, max_combinations))

    return selected_combinations


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


def align_coords(sample, target, case):
    """
    Translates and rotates coordinates based on the IDR case.
    
    Set of `target` coordinates should be for [[C], [N], [CA]].
    Where the `C` position is of the previous residue for N-IDR
    and next residue `C` for C-IDR.

    Parameters
    ----------
    sample : np.array
        Array format of the PDB of IDR in question. The return
        of `parse_pdb_to_array` from `libstructure`.
    
    target : np.array
        Set of 3D coordinates representing positions of C, N, CA
        fixed points to align to.
    
    case : str
        IDR case as mentioned above (N-IDR, C-IDR, Break-IDR).
    
    Returns
    -------
    sample : np.array
        New array format of the PDB now with coordinates rotated
        and translated about target.
    """
    atom_names = sample[:, col_name]
    idr_term_idx = {}
    counter = 0
    if case == disorder_cases[0]:  # N-IDR
        # In the case of N-IDR, we want to move relative to C-term
        for i, atom in enumerate(atom_names):
            # Use last residues of N-IDR for alignment
            j = len(atom_names) - 1 - i
            if counter == 1 and atom_names[j] == "N":
                idr_term_idx["N"] = j
                counter += 1
            elif counter == 0 and atom_names[j] == "CA":
                idr_term_idx["CA"] = j
                counter += 1
            elif counter == 2 and atom_names[j] == "C":
                idr_term_idx["C"] = j
                break
    elif case == disorder_cases[1]:  # break
        pass
    elif case == disorder_cases[2]:  # C-IDR
        # In the case of C-IDR, we want to move relative to N-term
        for i, atom in enumerate(atom_names):
            # Use first residues of C-IDR for alignment
            if counter == 5: break
            if atom == "N":
                idr_term_idx["N"] = i
                counter += 1
            elif atom == "CA":
                idr_term_idx["CA"] = i
                counter += 1
            elif atom == "C":
                idr_term_idx["C"] = i
                counter += 1

    idr_Cxyz = sample[idr_term_idx["C"]][cols_coords].astype(float).tolist()
    idr_Nxyz = sample[idr_term_idx["N"]][cols_coords].astype(float).tolist()
    idr_CAxyz = sample[idr_term_idx["CA"]][cols_coords].astype(float).tolist()
    idr_xyz = sample[:, cols_coords].astype(float)
    
    idr_coords = np.array([idr_Cxyz, idr_Nxyz, idr_CAxyz])
    
    centered_idr = idr_coords - idr_coords.mean(axis=0)
    centered_fld = target - target.mean(axis=0)
    
    covariance_matrix = np.dot(centered_idr.T, centered_fld)
    U, S, Vt = np.linalg.svd(covariance_matrix)
    rotation_matrix = np.dot(U, Vt)

    rotated_points = np.dot(idr_xyz, rotation_matrix)
    sample[:, cols_coords] = rotated_points.astype(str)
    
    translation_vector = \
        target[0] - sample[idr_term_idx["C"]][cols_coords].astype(float)

    for i, coords in enumerate(sample[:, cols_coords]):
        x = str(round(translation_vector[0] + float(coords[0]), 3))
        y = str(round(translation_vector[1] + float(coords[1]), 3))
        z = str(round(translation_vector[2] + float(coords[2]), 3))
        
        sample[i][col_x] = x
        sample[i][col_y] = y
        sample[i][col_z] = z
    
    return sample


def count_clashes(
    fragment,
    parent,
    case=None,
    max_clash=60,
    tolerance=0.6,      
    ):
    """
    Checks for steric clashes between two protein chains using vdW radii.

    Parameters
    ----------    
    fragment : np.array
        Array of the IDR fragment of interest
    
    parent : IDPConformerGenerator.Structure
        Structure of static protein chain of interest.
        Must already be built using `.build()`
    
    case : str, optional
        Disorder case of interest will change how clash is calculated
    
    max_clash : int, optional
        Integer number for maximum number of allowed clashes
    
    tolerance : float, optional
        Tolerance applicable to vdW clash validation in Angstroms
    
    Returns
    -------
    num_clashes : int
        Number of steric clashes determined using vdW radii
    
    True : Bool
        Too many clashes observed, not worthwhile continuing
    """  
    num_clashes = 0
    
    parent_atoms = parent.data_array[:, col_element]
    fragment_atoms = fragment[:, col_element]
    fragment_seq = fragment[:, col_resSeq]
    parent_coords = parent.data_array[:, cols_coords].astype(float)
    fragment_coords = fragment[:, cols_coords].astype(float)
    
    if case == disorder_cases[0]:
        # N-IDR, remove last 4 resiudes of fragment from consideration
        last_r = int(fragment_seq[len(fragment_seq) - 1])
        for i, _ in enumerate(fragment_seq):
            j = len(fragment_seq) - 1 - i
            try:  # In case the user wants to build less than 3 residues
                prev = int(fragment_seq[j - 1])
            except IndexError: 
                continue
            fragment_atoms = np.delete(fragment_atoms, j, axis=0)
            fragment_coords = np.delete(fragment_coords, j, axis=0)
            if last_r - prev == 5: break
    elif case == disorder_cases[1]:
        pass
    elif case == disorder_cases[2]:
        # C-IDR, remove first 4 residues of fragment from consideration
        first_r = int(fragment_seq[0])
        for i, _ in enumerate(fragment_seq):
            try:  # In case the user wants to build less than 3 residues
                next = int(fragment_seq[i + 1])
            except IndexError: 
                continue
            fragment_atoms = np.delete(fragment_atoms, i, axis=0)
            fragment_coords = np.delete(fragment_coords, i, axis=0)
            if next - first_r == 5: break
    
    # Loop through all pairs of atoms in the 2 protein chains
    for i, atom1 in enumerate(parent_atoms):
        for j, atom2 in enumerate(fragment_atoms):
            # calculate distance between atoms
            distance = calculate_distance(parent_coords[i], fragment_coords[j])
            
            # get vdW radii for each atom
            vdw_radius1 = vdW_radii_tsai_1999[atom1]
            vdw_radius2 = vdW_radii_tsai_1999[atom2]
            
            # Check if a steric clash is detected by comparing
            # distance between atoms to the sum of their vdW radii
            if num_clashes > 150:
                return True, fragment
            if distance < vdw_radius1 + vdw_radius2 + tolerance:
                num_clashes += 1
    
    return num_clashes, fragment


def psurgeon(idp_struc, fld_struc, case):
    """
    Protein surgeon grafts disordered regions onto folded structures.

    Parameters
    ----------
    idp_struc : Path or IDPConformerGenerator.Structure
        Donor conformer to graft
    
    case : str
        Case could be `nidr`, `cidr`, or `break` as defined above

    fld_struc : Path or IDPConformerGenerator.Structure
        Folded structure to be grafted on

    Returns
    -------
    new_struc_arr : np.ndarray (N x 16)
        Array format of the PDB file.
    """
    if type(fld_struc) is Path:
        fld = Structure(fld_struc)
        fld.build()
    if type(idp_struc) is Path:
        idr = Structure(idp_struc)
        idr.build()
    elif type(idp_struc) is tuple:
        nidr = Structure(idp_struc[0])
        cidr = Structure(idp_struc[1])
        nidr.build()
        cidr.build()
        
    fld_seq = fld.data_array[:, col_resSeq]
    fld_data_array = fld.data_array
    
    if case == disorder_cases[0]:
        idr_seq = idr.data_array[:, col_resSeq]
        idr_data_array = idr.data_array
        # N-IDR, remove last resiude of fragment
        # and remove the first residue of folded-domain
        for i, seq in enumerate(idr_seq):
            j = len(idr_seq) - 1 - i
            curr = idr_seq[j]
            prev = idr_seq[j - 1]
            idr_data_array = idr_data_array[:-1]
            if prev != curr:
                break
        
        idr._data_array = idr_data_array
        
        for i, seq in enumerate(fld_seq):
            curr = fld_seq[i]
            next = fld_seq[i + 1]
            fld_data_array = fld_data_array[1:]
            if next != curr:
                break
        
        fld._data_array = fld_data_array
        
        new_struc_arr = np.append(idr._data_array, fld._data_array, axis=0)
        
        new_serial = [str(i) for i in range(1, len(new_struc_arr) + 1)]
        new_struc_arr[:, col_serial] = new_serial
        new_struc_arr[:, col_chainID] = "A"
        new_struc_arr[:, col_segid] = "A"
    elif case == disorder_cases[1]:  # internal loop TBD
        pass
    elif case == disorder_cases[2]:
        idr_seq = idr.data_array[:, col_resSeq]
        idr_data_array = idr.data_array
        # C-IDR, remove last resiude of folded protein
        # and remove the first residue of C-IDR
        for i, seq in enumerate(fld_seq):
            j = len(fld_seq) - 1 - i
            curr = fld_seq[j]
            prev = fld_seq[j - 1]
            fld_data_array = fld_data_array[:-1]
            if prev != curr:
                break
        
        fld._data_array = fld_data_array
        
        for i, seq in enumerate(idr_seq):
            curr = idr_seq[i]
            next = idr_seq[i + 1]
            idr_data_array = idr_data_array[1:]
            if next != curr:
                break
        
        idr._data_array = idr_data_array

        # Fix residue connectivity issue
        last_residue_fld = int(fld.data_array[:, col_resSeq][-1])        
        curr_residue = last_residue_fld + 1
        idr_seq = idr.data_array[:, col_resSeq]
        for i, seq in enumerate(idr_seq):
            curr = seq
            idr._data_array[:, col_resSeq][i] = str(curr_residue)
            try:
                if idr_seq[i + 1] != curr:
                    curr_residue += 1
            except IndexError:
                break
        
        # Initialize and clean new structure
        new_struc_arr = np.append(fld._data_array, idr._data_array, axis=0)
        
        new_serial = [str(i) for i in range(1, len(new_struc_arr) + 1)]
        new_struc_arr[:, col_serial] = new_serial
        new_struc_arr[:, col_chainID] = "A"
        new_struc_arr[:, col_segid] = "A"
        
    elif case == disorder_cases[0] + disorder_cases[2]:
        nidr_seq = nidr.data_array[:, col_resSeq]
        cidr_seq = cidr.data_array[:, col_resSeq]
        nidr_data_array = nidr.data_array
        cidr_data_array = cidr.data_array
        
        # For cases where we have both C-IDR and N-IDR
        # idp_struc should be a list of tuple (N-IDR, C-IDR) paths
        # N-IDR, remove last resiude of fragment
        for i, seq in enumerate(nidr_seq):
            j = len(nidr_seq) - 1 - i
            curr = nidr_seq[j]
            prev = nidr_seq[j - 1]
            nidr_data_array = nidr_data_array[:-1]
            if prev != curr:
                break
        
        nidr._data_array = nidr_data_array
        
        # Remove first residue of folded domain
        for i, seq in enumerate(fld_seq):
            curr = fld_seq[i]
            next = fld_seq[i + 1]
            fld_data_array = fld_data_array[1:]
            if next != curr:
                break
        
        fld._data_array = fld_data_array
        
        new_struc_arr = np.append(nidr._data_array, fld._data_array, axis=0)
        new_struc_seq = new_struc_arr[:, col_resSeq]
        
        # Remove last residue of folded domain (now protein)
        for i, seq in enumerate(new_struc_seq):
            j = len(new_struc_seq) - 1 - i
            curr = new_struc_seq[j]
            prev = new_struc_seq[j - 1]
            new_struc_arr = new_struc_arr[:-1]
            if prev != curr:
                break
        
        # Remove first residue of C-IDR
        for i, seq in enumerate(cidr_seq):
            curr = cidr_seq[i]
            next = cidr_seq[i + 1]
            cidr_data_array = cidr_data_array[1:]
            if next != curr:
                break
            
        cidr._data_array = cidr_data_array
        
        # Fix residue connectivity issue
        last_residue_fld = int(new_struc_arr[:, col_resSeq][-1])        
        curr_residue = last_residue_fld
        cidr_seq = cidr.data_array[:, col_resSeq]
        for i, seq in enumerate(cidr_seq):
            curr = seq
            cidr._data_array[:, col_resSeq][i] = str(curr_residue)
            try:
                if cidr_seq[i + 1] != curr:
                    curr_residue += 1
            except IndexError:
                break
        
        new_struc_arr = np.append(new_struc_arr, cidr.data_array, axis=0)
        
        new_serial = [str(i) for i in range(1, len(new_struc_arr) + 1)]
        new_struc_arr[:, col_serial] = new_serial
        new_struc_arr[:, col_chainID] = "A"
        new_struc_arr[:, col_segid] = "A"
    
    return new_struc_arr


''' Deprecated code.
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
        # In the case of C-IDR, we want to move relative to first N
        for i, atom in enumerate(atom_names):
            if atom == "N":
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
    chain : Path or IDPConformerGenerator.Structure
        Chain of interest we want to rotate.
    
    case : string
        Disordered case of interest determines which point
        to rotate about.

    Returns
    -------
    idp : IDPConformerGenerator.Structure
        Structure object with the new rotated coordinates.
    """
    minrad = 0
    maxrad = 2 * np.pi
    # Select random angle to rotate
    angle = random.uniform(minrad, maxrad)
    
    if type(chain) is Path:
        idp = Structure(chain)
        idp.build()
    else:
        idp = chain
        
    atom_names = idp.data_array[:, col_name]
    
    if case == disorder_cases[0]:  # N-IDR
        for i, atom in enumerate(atom_names):
            if atom_names[len(atom_names) - 1 - i] == "N":
                index = len(atom_names) - 1 - i
                break
    elif case == disorder_cases[2]:  # C-IDR
        for i, atom in enumerate(atom_names):
            if atom == "N":
                index = i
                break
    
    # Create a 3x3 rotation matrix
    rotation_matrix = np.array(
        [[np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0, 0, 1]]
    )
    
    # Mark the fixed point and translate point as the origin
    point_xyz = idp.data_array[index][cols_coords].astype(float)
    coords = idp.data_array[:, cols_coords].astype(float) - point_xyz
    
    # Rotate coordinates about a random angle and translate point back
    rotated_coords = np.dot(rotation_matrix, coords.T).T
    rotated_coords += point_xyz
    
    idp.data_array[:, cols_coords] = rotated_coords.astype(str)
    
    return idp, chain  # returns original chain/path also
'''
