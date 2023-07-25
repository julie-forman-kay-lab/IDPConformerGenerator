"""
Client for building IDRs on PDB files in the cartesian coordinate space.

Methodology deviates from traditional IDP or beads-on-a-string FLDR/S approach.

Name: FLDR/S (Folded disordered region/structure sampling)
"""
import os
import random
from itertools import product

import numpy as np

from idpconfgen import Path
from idpconfgen.core.definitions import aa3to1, vdW_radii_tsai_1999
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs.libcalc import calc_torsion_angles
from idpconfgen.libs.libstructure import (
    Structure,
    col_chainID,
    col_name,
    col_resName,
    col_resSeq,
    col_segid,
    col_serial,
    col_x,
    col_y,
    col_z,
    cols_coords,
    structure_to_pdb,
    write_PDB,
    )


disorder_cases = {
    0: "N-IDR",
    1: "Linker-IDR",
    2: "C-IDR",
    }


def tolerance_calculator(tolerance):
    """
    Calculate the max number of tolerated spherical clashes and distance.

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
        
    max_clash = int(tolerance * 80)
    dist_tolerance = tolerance * 0.8
    
    return max_clash, dist_tolerance


def calculate_distance(coords1, coords2):
    """
    Calculate the distance between two 3D coordinates.
    
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


def calculate_angle(a, b, c):
    """
    Calculate angle between three 3D coordinates.

    Parameters
    ----------
    a : np.ndarray
    
    b : np.ndarray
    
    c : np.ndarray
    
    Return
    ------
    float angle in radians
    """
    ab = a - b
    cb = c - b

    # Calculate the dot product
    dot_product = np.dot(ab, cb)

    # Calculate the magnitudes of the vectors
    mag_ab = np.linalg.norm(ab)
    mag_cb = np.linalg.norm(cb)

    # Calculate the cosine of the angle between the vectors
    cos_angle = dot_product / (mag_ab * mag_cb)
    # Calculate the angle in radians
    angle = np.arccos(cos_angle)

    return angle


def consecutive_grouper(seq):
    """
    Use negative indexing to group together consecutive numbers.

    Reference
    ---------
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
            
    bounds = []
    for group in grouped:
        first = group[0]
        last = group[len(group) - 1]
        bounds.append((first, last + 1))
    
    return bounds


def create_combinations(lst, num_combinations):
    """
    Create unique combinations between list of lists.
    
    Made for N-IDR and C-IDR combinations. Where list1 = N-IDR paths,
    and list2 = C-IDR paths. Itertools product is used here because
    order matters.
    
    Parameters
    ----------
    lst : list
        Can be a list of lists, but normally in the order of
        N-IDR, Linker-IDR, C-IDR paths
    
    num_combinations : int
    
    Return
    ------
    selected_combinations : list
        List of lists of different combinations as follows:
        [[item from list 1, item from list 2], ...]
    """
    all_combinations = list(product(*lst))
    max_combinations = len(all_combinations)

    selected_combinations = \
        random.sample(all_combinations, min(num_combinations, max_combinations))

    return np.array(selected_combinations).tolist()


def create_all_combinations(folder, nconfs):
    """
    Generate combinations of all cases.

    Parameters
    ----------
    folder : Path
        Temporary folder with all IDR cases.
    
    nconfs : int
        Desired number of list elements based on number of conformers.
    
    Returns
    -------
    combinations : list
        List of all combinations with first element being path for N-IDR,
        then Linker-IDR, and last element C-IDR if necessary.
    """
    nidr_path = folder.joinpath(disorder_cases[0])
    lidr_path = folder.joinpath(disorder_cases[1])
    cidr_path = folder.joinpath(disorder_cases[2])
    nidr_files = []
    lidr_combinations = []
    cidr_files = []
    
    if os.path.isdir(nidr_path):
        nidr_confs = os.listdir(nidr_path)
        nidr_files = [Path(nidr_path.joinpath(cpath)) for cpath in nidr_confs]
    if os.path.isdir(lidr_path):
        lidr_cases_dir = os.listdir(lidr_path)
        lidr_confs_lst = []
        for c, cpath in enumerate(lidr_cases_dir):
            lidr_matches = lidr_path.joinpath(cpath + f"/{c}_match")
            lidr_confs = os.listdir(lidr_matches)
            lidr_files = [Path(lidr_matches.joinpath(fpath)) for fpath in lidr_confs]  # noqa: E501
            lidr_confs_lst.append(lidr_files)
        lidr_combinations = create_combinations(lidr_confs_lst, nconfs)
    if os.path.isdir(cidr_path):
        cidr_confs = os.listdir(cidr_path)
        cidr_files = [Path(cidr_path.joinpath(cpath)) for cpath in cidr_confs]
    
    if len(nidr_files) and len(cidr_files) and len(lidr_combinations) > 0:
        return create_combinations([nidr_files, lidr_combinations, cidr_files], nconfs)  # noqa: E501
    elif len(nidr_files) and len(lidr_combinations) > 0:
        return create_combinations([nidr_files, lidr_combinations], nconfs)
    elif len(cidr_files) and len(lidr_combinations) > 0:
        return create_combinations([lidr_combinations, cidr_files], nconfs)
    elif len(nidr_files) and len(cidr_files) > 0:
        return create_combinations([nidr_files, cidr_files], nconfs)
    elif len(nidr_files) > 0:
        return nidr_files
    elif len(lidr_combinations) > 0:
        return lidr_combinations
    elif len(cidr_files) > 0:
        return cidr_files


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
            fld_seqs.append(''.join(aa3to1.get(f) for f in data[:, col_resName][fld_idx].tolist()))  # noqa: E501
        
        return fld_seqs
    
    return


def align_coords(sample, target, case):
    """
    Translate and rotate coordinates based on the IDR case.
    
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
        
        For Linker-IDR it would be (F, L), (CA, N) where "F" and "L"
        are the positions of "C" for the first and last bit of the
        chain break respectively and "CA" is for alignment to
        first part of break while "N" is for the last part.
    
    case : str
        IDR case as mentioned above (N-IDR, C-IDR, Linker-IDR).
    
    Returns
    -------
    sample : np.array
        New array format of the PDB now with coordinates rotated
        and translated about target.
    """
    atom_names = sample[:, col_name]
    res_seq = sample[:, col_resSeq].astype(int)

    first_seq = res_seq[0]
    last_seq = res_seq[-1]

    idr_term_idx = {}

    if case == disorder_cases[0]:  # N-IDR
        # In the case of N-IDR, we want to move relative to C-term
        for i, _atom in enumerate(atom_names):
            # Use last residues of N-IDR for alignment
            j = len(atom_names) - 1 - i
            seq = res_seq[j]

            if seq == last_seq and atom_names[j] == "N":
                idr_term_idx["N"] = j
            elif seq == last_seq and atom_names[j] == "CA":
                idr_term_idx["CA"] = j
            elif seq == last_seq - 1 and atom_names[j] == "C":
                idr_term_idx["C"] = j
            elif seq == last_seq - 2:
                break
    elif case == disorder_cases[2]:  # C-IDR
        # We want to move relative to N-term of fragment
        for i, atom in enumerate(atom_names):
            seq = res_seq[i]

            # Use first residues of C-IDR for alignment
            if seq == first_seq + 1 and atom == "N":
                idr_term_idx["N"] = i
            elif seq == first_seq + 1 and atom == "CA":
                idr_term_idx["CA"] = i
            elif seq == first_seq and atom == "C":
                idr_term_idx["C"] = i
            elif seq == first_seq + 2:
                break

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


def next_seeker(
        cterm_idr,
        nterm_idr_lib,
        max_clash,
        tolerance,
        output_folder,
        ):
    """
    Next-seeker protocol to find possible matches that will close the break.

    Parameters
    ----------
    cterm_idr : str or Path
        Path to the IDR chain of interest from C-term side of break
    
    nterm_idr_lib : list
        List of paths of IDR chains from N-term side of break
    
    max_clash : int
        Integer number for maximum number of allowed clashes
    
    tolerance : float
        Tolerance applicable to vdW clash validation in Angstroms
    
    output_folder : str
        Output folder to store conformers with matches
    
    Returns
    -------
    False
        If no matches have been found. Otherwise write them to output folder.
    """
    matches = 0
    idr_struc = Structure(Path(cterm_idr))
    idr_struc.build()
    idr_arr = idr_struc.data_array
    idr_name = idr_arr[:, col_name]
    idr_coords = idr_arr[:, cols_coords].astype(float)
    idr_resseq = idr_arr[:, col_resSeq].astype(int)
    
    idr_C = []
    idr_CA = []
    idr_O = []
    idr_res = []
    for n, name in enumerate(idr_name):
        if name == 'C':
            idr_C.append(idr_coords[n])
            idr_res.append(idr_resseq[n])
        elif name == 'CA':
            idr_CA.append(idr_coords[n])
        elif name == 'O':
            idr_O.append(idr_coords[n])
    
    for nterm_idr in nterm_idr_lib:
        nterm_idr_struc = Structure(Path(nterm_idr))
        nterm_idr_struc.build()
        nterm_idr_arr = nterm_idr_struc.data_array
        nterm_idr_name = nterm_idr_arr[:, col_name]
        nterm_idr_coords = nterm_idr_arr[:, cols_coords].astype(float)
        
        nterm_idr_N = []
        nterm_idr_CA = []
        for n, name in enumerate(nterm_idr_name):
            if name == 'N':
                nterm_idr_N.append(nterm_idr_coords[n])
            elif name == 'CA':
                nterm_idr_CA.append(nterm_idr_coords[n])
        
        for i, curr_c in enumerate(idr_C):
            try:
                next_n = nterm_idr_N[i + 1]
                next_ca = nterm_idr_CA[i + 1]
            except IndexError:
                break
            CN_dist = calculate_distance(curr_c, next_n)
            CCA_dist = calculate_distance(curr_c, next_ca)
            CACN_ang = calculate_angle(idr_CA[i], curr_c, next_n)
            CACNCA_coords = np.array([idr_CA[i], curr_c, next_n, next_ca])
            omega = calc_torsion_angles(CACNCA_coords)
            # Here is a set of geometric checks to ensure we have closure
            # Refer to distances and angles in `core/build_definitions.py`
            # |omega| angle must be greater than 150 deg
            if 1.32 <= CN_dist <= 1.56 and 1.91 <= CACN_ang <= 2.15 and 2.2 <= CCA_dist <= 2.7 and np.abs(omega) >= 2.61:  # noqa: E501
                term_residue = idr_res[i]
                
                idr_list = []
                for p, _ in enumerate(idr_resseq):
                    next = idr_resseq[p + 1]
                    idr_list.append(idr_arr[p])
                    nterm_idr_arr = nterm_idr_arr[1:]
                    if next == term_residue + 1:
                        break
                idr_arr = np.array(idr_list)
                    
                nterm_idr_struc._data_array = nterm_idr_arr
                
                clashes, _ = count_clashes(
                    idr_arr,
                    nterm_idr_struc,
                    disorder_cases[1],
                    max_clash,
                    tolerance,
                    )
                
                if type(clashes) is int:
                    nterm_idr_list = nterm_idr_arr.tolist()
                    final_struc_arr = np.array(idr_list + nterm_idr_list)
                    final_struc_name = final_struc_arr[:, col_name]
                    final_struc_res = final_struc_arr[:, col_resSeq].astype(int)
                    H_idx = -1  # for cases like Proline without "H"
                    for idx, name in enumerate(final_struc_name):
                        if final_struc_res[idx] == term_residue:
                            if name == 'O':
                                O_idx = idx
                        elif final_struc_res[idx] == term_residue + 1:
                            if name == 'H':
                                H_idx = idx
                        elif final_struc_res[idx] == term_residue + 2:
                            break
                    
                    # Fix the position of the Carbonyl O and Nitrogen H
                    CO_length = calculate_distance(curr_c, idr_O[i])
                    CAC_O_vec = idr_CA[i] - curr_c
                    NC_O_vec = next_n - curr_c
                    O_angle = np.arccos(np.dot(CAC_O_vec, NC_O_vec) / (np.linalg.norm(CAC_O_vec) * np.linalg.norm(NC_O_vec)))  # noqa: E501
                    O_vector = CO_length * np.sin(O_angle / 2) * (CAC_O_vec / np.linalg.norm(CAC_O_vec)) + CO_length * np.sin(O_angle / 2) * (NC_O_vec / np.linalg.norm(NC_O_vec))  # noqa: E501
                    new_O_xyz = curr_c - O_vector
  
                    final_struc_arr[:, col_x][O_idx] = str(new_O_xyz[0])
                    final_struc_arr[:, col_y][O_idx] = str(new_O_xyz[1])
                    final_struc_arr[:, col_z][O_idx] = str(new_O_xyz[2])
                    
                    if H_idx >= 0:
                        # Bond length also taken from
                        # `core/build_definitions.py`
                        NH_length = 1.0
                        CN_H_vec = curr_c - next_n
                        CAN_H_vec = next_ca - next_n
                        H_angle = np.arccos(np.dot(CN_H_vec, CAN_H_vec) / (np.linalg.norm(CN_H_vec) * np.linalg.norm(CAN_H_vec)))  # noqa: E501
                        H_vector = NH_length * np.sin(H_angle / 2) * (CN_H_vec / np.linalg.norm(CN_H_vec)) + NH_length * np.sin(H_angle / 2) * (CAN_H_vec / np.linalg.norm(CAN_H_vec))  # noqa: E501
                        new_H_xyz = next_n - H_vector
                        
                        final_struc_arr[:, col_x][H_idx] = str(new_H_xyz[0])
                        final_struc_arr[:, col_y][H_idx] = str(new_H_xyz[1])
                        final_struc_arr[:, col_z][H_idx] = str(new_H_xyz[2])
                    
                    final_struc = structure_to_pdb(final_struc_arr)
                    cterm_idr_stem = Path(cterm_idr).stem
                    nterm_idr_stem = Path(nterm_idr).stem
                    matches += 1
                    write_PDB(final_struc, str(output_folder) + f"/{cterm_idr_stem}+{nterm_idr_stem}.pdb")  # noqa: E501

    return matches


def count_clashes(
        fragment,
        parent,
        case=None,
        max_clash=40,
        tolerance=0.4,
        ):
    """
    Check for steric clashes between two protein chains using vdW radii.

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
    num_clashes : int or Bool
        Number of steric clashes determined using vdW radii
    
    fragment : np.array
        Array of the IDR fragment of interest
    """
    num_clashes = 0
    
    # Not using `col_element` because some PDB files may not have that column
    # take the first character of `col_name` instead
    parent_atoms = parent.data_array[:, col_name]
    fragment_atoms = fragment[:, col_name]
    fragment_seq = fragment[:, col_resSeq].astype(int)
    parent_coords = parent.data_array[:, cols_coords].astype(float)
    fragment_coords = fragment[:, cols_coords].astype(float)
    
    first_r = fragment_seq[0]
    last_r = fragment_seq[-1]
    
    if case == disorder_cases[0] or case == disorder_cases[1]:
        # N-IDR or Linker-IDR, remove last 2 resiudes of fragment
        # from consideration
        for i, _ in enumerate(fragment_seq):
            j = len(fragment_seq) - 1 - i
            try:  # In case the user wants to build less than 3 residues
                prev = fragment_seq[j - 1]
            except IndexError:
                continue
            fragment_atoms = fragment_atoms[:-1]
            fragment_coords = fragment_coords[:-1]
            if last_r - prev == 3:
                break
    elif case == disorder_cases[2]:
        # C-IDR, remove first 2 residues of fragment from consideration
        for i, _ in enumerate(fragment_seq):
            try:  # In case the user wants to build less than 3 residues
                next = fragment_seq[i + 1]
            except IndexError:
                continue
            fragment_atoms = fragment_atoms[1:]
            fragment_coords = fragment_coords[1:]
            if next - first_r == 3:
                break

    # Loop through all pairs of atoms in the 2 protein chains
    # Take the first character of atom name as the element ID
    for i, atom1 in enumerate(parent_atoms):
        a1 = atom1[0]
        for j, atom2 in enumerate(fragment_atoms):
            a2 = atom2[0]
            # calculate distance between atoms
            distance = calculate_distance(parent_coords[i], fragment_coords[j])
            
            # get vdW radii for each atom
            vdw_radius1 = vdW_radii_tsai_1999[a1]
            vdw_radius2 = vdW_radii_tsai_1999[a2]
            
            # Check if a steric clash is detected by comparing
            # distance between atoms to the sum of their vdW radii
            if num_clashes > max_clash:
                return True, fragment
            if distance < vdw_radius1 + vdw_radius2 + tolerance:
                num_clashes += 1
    
    return num_clashes, fragment


def psurgeon(idp_lst, fld_struc, case, ranges):
    """
    Protein surgeon grafts disordered regions onto folded structures.

    Parameters
    ----------
    idp_lst : list of Path
        Always in the order of [N-IDR, Linker-IDR, C-IDR]
    
    case : list of str
        Case could be `N-IDR`, `C-IDR`, or `Linker-IDR` as defined above

    fld_struc : Path or IDPConformerGenerator.Structure
        Folded structure to be grafted on
    
    ranges : list of tuple of int
        For Linker-IDR, what residue ranges for the chain break

    Returns
    -------
    new_struc_arr : np.ndarray (N x 16)
        Array format of the PDB file.
    """
    if type(fld_struc) is Path:
        fld = Structure(fld_struc)
        fld.build()
        fld_seq = fld.data_array[:, col_resSeq]
        fld_data_array = fld.data_array
    else:
        fld_seq = fld_struc.data_array[:, col_resSeq]
        fld_data_array = fld_struc.data_array
    
    if disorder_cases[0] in case:
        # N-IDR, remove last resiude of fragment
        # and remove the first residue of folded-domain
        nidr = Structure(idp_lst[0])
        nidr.build()
        nidr_seq = nidr.data_array[:, col_resSeq]
        nidr_data_array = nidr.data_array
        
        for i, _ in enumerate(nidr_seq):
            j = len(nidr_seq) - 1 - i
            curr = nidr_seq[j]
            prev = nidr_seq[j - 1]
            nidr_data_array = nidr_data_array[:-1]
            if prev != curr:
                break
        
        nidr._data_array = nidr_data_array
        
        for i, seq in enumerate(fld_seq):
            next = fld_seq[i + 1]
            fld_data_array = fld_data_array[1:]
            if next != seq:
                break
        
        fld._data_array = fld_data_array
        idp_lst.pop(0)
        ranges.pop(0)
        new_struc_arr = np.append(nidr.data_array, fld.data_array, axis=0)

    if disorder_cases[2] in case:
        # C-IDR, remove last resiude of folded protein
        # and remove the first residue of C-IDR
        cidr = Structure(idp_lst[-1])
        cidr.build()
        cidr_seq = cidr.data_array[:, col_resSeq]
        cidr_data_array = cidr.data_array
        
        for i, seq in enumerate(cidr_seq):
            next = cidr_seq[i + 1]
            cidr_data_array = cidr_data_array[1:]
            if next != seq:
                break
        
        cidr._data_array = cidr_data_array
        
        for i, _ in enumerate(fld_seq):
            j = len(fld_seq) - 1 - i
            curr = fld_seq[j]
            prev = fld_seq[j - 1]
            fld_data_array = fld_data_array[:-1]
            if prev != curr:
                break
        
        fld._data_array = fld_data_array
        
        # Fix residue connectivity issue
        last_residue_fld = int(fld.data_array[:, col_resSeq][-1])
        curr_residue = last_residue_fld + 1
        cidr_seq = cidr.data_array[:, col_resSeq]
        for i, seq in enumerate(cidr_seq):
            curr = seq
            cidr._data_array[:, col_resSeq][i] = str(curr_residue)
            try:
                if cidr_seq[i + 1] != curr:
                    curr_residue += 1
            except IndexError:
                break
        
        if disorder_cases[0] in case:
            new_struc_arr = np.append(nidr.data_array, fld.data_array, axis=0)
            new_struc_arr = np.append(new_struc_arr, cidr.data_array, axis=0)
        else:
            new_struc_arr = np.append(fld.data_array, cidr.data_array, axis=0)
        idp_lst.pop(-1)
        ranges.pop(-1)
    
    if disorder_cases[1] in case:
        if len(case) == 1:
            new_struc_arr = fld_data_array
        
        new_struc_seq = new_struc_arr[:, col_resSeq].astype(int)
        first_struc_seq = new_struc_seq[0]
        
        for idx, minmax in enumerate(ranges):
            lower = minmax[0]
            upper = minmax[1]
            new_struc_lst = new_struc_arr.tolist()
            
            idr = Structure(idp_lst[0][idx])
            idr.build()
            
            idr_seq = idr.data_array[:, col_resSeq].astype(int)
            idr_data_array = idr.data_array
            last_res = idr_seq[-1]
            # Linker-IDR, remove first and last residue of fragment
            # and the first residues on either side of the chain break
            for i, seq in enumerate(idr_seq):
                next = idr_seq[i + 1]
                j = len(idr_seq) - 1 - i
                rev = idr_seq[j]
                
                if seq == 1:
                    idr_data_array = idr_data_array[1:]
                if rev == last_res:
                    idr_data_array = idr_data_array[:-1]
                if seq > 2 and rev < last_res - 1:
                    break
                
            idr_data_lst = idr_data_array.tolist()
            surrounding_data_list = []
            if disorder_cases[0] not in case:
                actual_lower = lower + first_struc_seq - 1
                actual_upper = upper + first_struc_seq
            else:
                actual_lower = lower
                actual_upper = upper + 1
            found = False
            for s, seq in enumerate(new_struc_seq):
                if seq == actual_lower or seq == actual_upper:
                    if found is False:
                        found = True
                        insert_idx = s
                    continue
                else:
                    surrounding_data_list.append(new_struc_lst[s])
            # Needs to be reinitialized because we changed the array
            for r, row in enumerate(idr_data_lst):
                idr_seq = int(row[col_resSeq])
                row[col_resSeq] = str(idr_seq + actual_lower - 2)
                surrounding_data_list.insert(insert_idx + r, row)
            
            new_struc_arr = np.array(surrounding_data_list)
            new_struc_seq = new_struc_arr[:, col_resSeq].astype(int)
    
    # Fix serial numbers and chain IDs
    new_serial = [str(i) for i in range(1, len(new_struc_arr) + 1)]
    new_struc_arr[:, col_serial] = new_serial
    new_struc_arr[:, col_chainID] = "A"
    new_struc_arr[:, col_segid] = "A"
    
    return new_struc_arr
