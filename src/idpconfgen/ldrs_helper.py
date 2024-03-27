"""
Client for building IDRs on PDB files in the cartesian coordinate space.

Methodology deviates from traditional IDP or beads-on-a-string FLDR/S approach.

Name: FLDR/S (Folded disordered region/structure sampling)
"""
import os
import random
from itertools import combinations, product

import numpy as np

from idpconfgen import Path
from idpconfgen.core.definitions import (
    aa3to1,
    vdW_radii_ionic_CRC82,
    vdW_radii_tsai_1999,
    )
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs.libcalc import calc_torsion_angles
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libparse import convert_tuples_to_lists
from idpconfgen.libs.libstructure import (
    Structure,
    col_chainID,
    col_element,
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


def combinations_clash_check(selected):
    """
    Clash-check IDRs within a combination against each other.

    Parameters
    ----------
    selected : list
        List of lists of different combinations as follows:
        [[item from list 1, item from list 2], ...]
    """
    if len(selected) == 2:
        first = selected[0]
        second = selected[1]
        # [L-IDR, C-IDR] case
        if type(first) is tuple:
            second_struc = Structure(second)
            second_struc.build()
            if len(first) > 1:
                # Check all L-IDR against C-IDR
                for l in first:    # noqa: E741
                    l_struc = Structure(l)
                    l_struc.build()
                    l_arr = l_struc.data_array
                    l_c, _ = count_clashes(l_arr, second_struc, tolerance=1.0)
                    if type(l_c) is bool:
                        return False
                l_combos = list(combinations(first, 2))
                # Check all L-IDR against each other
                for combo in l_combos:
                    l1 = Structure(combo[0])
                    l1.build()
                    l2 = Structure(combo[1])
                    l2.build()
                    lx2, _ = count_clashes(l1.data_array, l2, tolerance=1.0)
                    if type(lx2) is bool:
                        return False
            else:
                lidr = first[0]
                lidr_struc = Structure(lidr)
                lidr_arr = lidr_struc.data_array
                l_c, _ = count_clashes(lidr_arr, second_struc, tolerance=1.0)
                if type(l_c) is bool:
                    return False
        # [N-IDR, L-IDR] case
        elif type(second) is tuple:
            first_struc = Structure(second)
            first_struc.build()
            if len(second) > 1:
                # Check all L-IDR against N-IDR
                for l in second:    # noqa: E741
                    l_struc = Structure(l)
                    l_struc.build()
                    l_arr = l_struc.data_array
                    n_l, _ = count_clashes(l_arr, first_struc, tolerance=1.0)
                    if type(n_l) is bool:
                        return False
                l_combos = list(combinations(second, 2))
                # Check all L-IDR against each other
                for combo in l_combos:
                    l1 = Structure(combo[0])
                    l1.build()
                    l2 = Structure(combo[1])
                    l2.build()
                    lx2, _ = count_clashes(l1.data_array, l2, tolerance=1.0)
                    if type(lx2) is bool:
                        return False
            else:
                lidr = second[0]
                lidr_struc = Structure(lidr)
                lidr_arr = lidr_struc.data_array
                n_l, _ = count_clashes(lidr_arr, first_struc, tolerance=1.0)
                if type(n_l) is bool:
                    return False
        # [N-IDR, C-IDR] case
        elif type(first) is Path and type(second) is Path:
            first_struc = Structure(first)
            first_struc.build()
            first_arr = first_struc.data_array
            second_struc = Structure(second)
            second_struc.build()
            n_c, _ = count_clashes(first_arr, second_struc, tolerance=1.0)
            if type(n_c) is bool:
                return False
    # [N-IDR, L-IDR, C-IDR] case
    elif len(selected) == 3:
        nidr = selected[0]
        nidr_struc = Structure(nidr)
        nidr_struc.build()
        cidr = selected[2]
        cidr_struc = Structure(cidr)
        cidr_struc.build()
        
        lidr = selected[1]
        if len(lidr) > 1:
            # Check all L-IDR against N-IDR and C-IDR
            for l in lidr:    # noqa: E741
                l_struc = Structure(l)
                l_struc.build()
                l_arr = l_struc.data_array
                n_l, _ = count_clashes(l_arr, nidr_struc, tolerance=1.0)
                l_c, _ = count_clashes(l_arr, cidr_struc, tolerance=1.0)
                if type(n_l) is bool or type(l_c) is bool:
                    return False
            l_combos = list(combinations(lidr, 2))
            # Check all L-IDR against each other
            for combo in l_combos:
                l1 = Structure(combo[0])
                l1.build()
                l2 = Structure(combo[1])
                l2.build()
                lx2, _ = count_clashes(l1.data_array, l2, tolerance=1.0)
                if type(lx2) is bool:
                    return False
        # We just have 1 L-IDR
        else:
            lidr_struc = Structure(lidr[0])
            lidr_struc.build()
            lidr_arr = lidr_struc.data_array
            nidr_arr = nidr_struc.data_array
            n_i, _ = count_clashes(nidr_arr, lidr_struc, tolerance=1.0)
            i_c, _ = count_clashes(lidr_arr, cidr_struc, tolerance=1.0)
            n_c, _ = count_clashes(nidr_arr, cidr_struc, tolerance=1.0)
            
            if type(n_i) is bool or type(i_c) is bool or type(n_c) is bool:
                return False
    
    return selected


def create_combinations(lst, num_combinations, ncores=1):
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
        Number of combinations to create.
    
    ncores : int
        Number of cores to use for clash-checking against multiple IDR
        combinations.
    
    Return
    ------
    selected_combinations : list
        List of lists of different combinations as follows:
        [[item from list 1, item from list 2], ...]
    """
    all_combinations = list(product(*lst))
    max_combinations = len(all_combinations)

    selected_combinations = \
        random.sample(all_combinations, min(num_combinations, max_combinations))  # noqa: E501
    
    if len(lst) == 1:
        return convert_tuples_to_lists(selected_combinations)
    else:
        execute_clash = pool_function(combinations_clash_check, selected_combinations, ncores=ncores)  # noqa: E501
        passed = []
        for result in execute_clash:
            if result is not False:
                passed.append(result)
        
        while len(passed) < num_combinations:
            selected_combinations = \
                random.sample(all_combinations, min(num_combinations, max_combinations))  # noqa: E501
            execute_clash = pool_function(combinations_clash_check, selected_combinations, ncores=ncores)  # noqa: E501
            for result in execute_clash:
                if result is not False:
                    passed.append(result)
            try:
                passed = list(set(passed))
            except TypeError:
                continue
        
        if len(passed) > num_combinations:
            passed = passed[0:num_combinations]

        return convert_tuples_to_lists(passed)


def create_all_combinations(folder, chains, nconfs, ncores=1):
    """
    Generate combinations of all cases.

    Parameters
    ----------
    folder : Path
        Temporary folder with all IDR cases.
    
    chains : list
        List of all chain IDs modeled.
    
    nconfs : int
        Desired number of list elements based on number of conformers.
    
    ncores : int
        Number of cores to perform multiprocessing with.
        Exists for clash-checking step for all combinations.
    
    Returns
    -------
    combinations : dict
        Dictionary of all combinations for each chain being the key and
        values as lists of combinations with first element being path for
        N-IDR, then Linker-IDR, and last element C-IDR if necessary.
    """
    combinations = {}
    for c in chains:
        chain_path = folder.joinpath(f"chain{c}")
        nidr_path = chain_path.joinpath(disorder_cases[0])
        lidr_path = chain_path.joinpath(disorder_cases[1])
        cidr_path = chain_path.joinpath(disorder_cases[2])
        nidr_files = []
        lidr_combinations = []
        cidr_files = []

        if os.path.isdir(nidr_path):
            nidr_confs = os.listdir(nidr_path)
            nidr_files = [Path(nidr_path.joinpath(cpath)) for cpath in nidr_confs]  # noqa: E501
        if os.path.isdir(lidr_path):
            lidr_cases_dir = os.listdir(lidr_path)
            lidr_confs_lst = []
            for i, cpath in enumerate(lidr_cases_dir):
                lidr_matches = lidr_path.joinpath(cpath + f"/{i}_match")
                lidr_confs = os.listdir(lidr_matches)
                lidr_files = [Path(lidr_matches.joinpath(fpath)) for fpath in lidr_confs]  # noqa: E501
                lidr_confs_lst.append(lidr_files)
            lidr_combinations = create_combinations(lidr_confs_lst, nconfs)
        if os.path.isdir(cidr_path):
            cidr_confs = os.listdir(cidr_path)
            cidr_files = [Path(cidr_path.joinpath(cpath)) for cpath in cidr_confs]  # noqa: E501

        if len(nidr_files) and len(cidr_files) and len(lidr_combinations) > 0:
            combinations[c] = create_combinations([nidr_files, lidr_combinations, cidr_files], nconfs, ncores)  # noqa: E501
        elif len(nidr_files) and len(lidr_combinations) > 0:
            combinations[c] = create_combinations([nidr_files, lidr_combinations], nconfs, ncores)  # noqa: E501
        elif len(cidr_files) and len(lidr_combinations) > 0:
            combinations[c] = create_combinations([lidr_combinations, cidr_files], nconfs, ncores)  # noqa: E501
        elif len(nidr_files) and len(cidr_files) > 0:
            combinations[c] = create_combinations([nidr_files, cidr_files], nconfs, ncores)  # noqa: E501
        elif len(nidr_files) > 0:
            combinations[c] = nidr_files[0:nconfs]
        elif len(lidr_combinations) > 0:
            combinations[c] = lidr_combinations
        elif len(cidr_files) > 0:
            combinations[c] = cidr_files[0:nconfs]
    
    return combinations


def break_check(fdata, membrane=False):
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
    
    membrane : bool
        If the PDB file containes a membrane built using CHARMM-GUI's
        bilayer builder or not.
        Defaults to False. No membrane.
    
    Return
    ------
    fld_seqs : list
        List of FASTA sequence of folded regions in the sequence of fdata.
    """
    structure = Structure(fdata)
    structure.build()
    if membrane:
        fld_pro = []
        fld_arr = structure.data_array.tolist()
        fld_segid = structure.data_array[:, col_segid]
        
        for i, segid in enumerate(fld_segid):
            if "PRO" in segid:
                fld_pro.append(fld_arr[i])
        fld_pro = np.array(fld_pro)
        structure._data_array = fld_pro
        
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
            
            res3 = data[:, col_resName][fld_idx].tolist()
            res1 = []
            for rn in res3:
                # Most likely exists a better solution than this.
                # To account for residues with different 3 letter codes.
                if rn == "HSD" or rn == "HIP":
                    rn = "HIS"
                res1.append(aa3to1.get(rn))
            
            fld_seqs.append(''.join(res1))
        
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
    
    # Check for reflection and correct if necessary
    if np.linalg.det(Vt.T @ U.T) < 0:
        Vt[-1, :] *= -1
    # Compute the rotation matrix
    rotation_matrix = np.dot(Vt.T, U.T)

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
        dtype=np.float32,
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
    
    dtype : data type, optional
        Data type for the array for clash-check matrix.
        Defaults to np.float32
        Can be np.float16, np.float32, np.float64, np.long, etc.
    
    Returns
    -------
    num_clashes : int or Bool
        Number of steric clashes determined using vdW radii
    
    fragment : np.array
        Array of the IDR fragment of interest
    """
    num_clashes = 0
    
    # PDB must have element column
    parent_atoms = parent.data_array[:, col_element]
    fragment_atoms = fragment[:, col_element]
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
    
    # Use matrix solution instead of nested for-loop to improve
    # clash checking speed at least 200-fold
    
    # Create dictionary with all possible radii
    all_vdw_radii = {**vdW_radii_tsai_1999, **vdW_radii_ionic_CRC82}
    
    # Calculate all distances
    # In case there's too many atoms, we cannot use float64 or else we will
    # run into a OOM (out of memory) error
    distances = np.sqrt(
        np.sum(
            (parent_coords[:, np.newaxis, :] - fragment_coords) ** 2,
            axis=-1,
            dtype=dtype
            ),
        dtype=dtype
        )
    
    # Get all radii
    vdw_radii1 = np.array([all_vdw_radii[atom] for atom in parent_atoms])
    vdw_radii2 = np.array([all_vdw_radii[atom] for atom in fragment_atoms])
    vdw_radii1 = vdw_radii1[:, np.newaxis]
    vdw_radii2 = vdw_radii2[np.newaxis, :]
    
    # Get maximum radii and compare with distances
    sum_radii = vdw_radii1 + vdw_radii2 + tolerance
    clashes = distances < sum_radii
    num_clashes = np.sum(clashes)
    
    if num_clashes > max_clash:
        return True, fragment
    
    return int(num_clashes), fragment


def psurgeon(
        idp_lst,
        fld_struc,
        case,
        ranges,
        membrane=False,
        skipped_chains=None
        ):
    """
    Protein surgeon grafts disordered regions onto folded structures.

    Parameters
    ----------
    idp_lst : dict of list of Path
        Always in the order of [N-IDR, Linker-IDR, C-IDR]
        Keys are the chain ID
    
    case : dict list of str
        Case could be `N-IDR`, `C-IDR`, or `Linker-IDR` as defined above
        Keys are the chain ID

    fld_struc : Path or IDPConformerGenerator.Structure
        Folded structure to be grafted on
    
    ranges : dict of list of tuple of int
        For Linker-IDR, what residue ranges for the chain break
        Keys are the chain ID
    
    membrane : Bool
        Whether or not a membrane exists within the template structure.
        Carried over from user-input in `ldrs` sub-client.
        Defaults to False.
    
    skipped_chains : list
        Any skipped chains we're building can be appended
        at the end of our new structure.
        Optional. Defaults to None.

    Returns
    -------
    new_struc_arr : np.ndarray (N x 16)
        Array format of the PDB file.
    """
    if type(fld_struc) is Path:
        fld = Structure(fld_struc)
        fld.build()
        fld_seq = fld.data_array[:, col_resSeq]
        fld_chain = fld.data_array[:, col_chainID]
        fld_data_array = fld.data_array
    else:
        fld_seq = fld_struc.data_array[:, col_resSeq]
        fld_chain = fld_struc.data_array[:, col_chainID]
        fld_data_array = fld_struc.data_array
    
    if membrane:
        fld_pro = []
        fld_arr = fld_data_array.tolist()
        fld_segid = fld_data_array[:, col_segid]
        
        for i, segid in enumerate(fld_segid):
            if "PRO" in segid:
                fld_pro.append(fld_arr[i])
        
        # Split the protein and membrane
        fld_pro = np.array(fld_pro)
        fld_data_array = fld_pro
        fld_seq = fld_pro[:, col_resSeq]
        fld_chain = fld_pro[:, col_chainID]
    
    # If only lists are given instead of dictionary
    # convert them to dictionary and process accordingly
    # This is for scripts made for v0.7.2 and earlier
    if type(idp_lst) is list and type(case) is list and type(ranges) is list:  # noqa: E501
        chain = fld_chain[0]
        # Accounting for emtpy chain case (should not happen)
        if chain is None or chain == '':
            chain = "A"
        case = {chain: case}
        ranges = {chain: ranges}
        idp_lst = {chain: idp_lst}
    kCases = set(case.keys())
    kRanges = set(ranges.keys())
    assert kCases == kRanges

    new_struc_arr_complete = []
    for chain in case:
        c = case[chain]
        r = ranges[chain]
        fld_data_lst = fld_data_array.tolist()
        fld_seq_lst = fld_seq.tolist()
        fld_data_seg_lst = []
        fld_seq_seg_lst = []

        for i, ch in enumerate(fld_chain):
            if ch == chain or ch is None or ch == '':
                fld_data_seg_lst.append(fld_data_lst[i])
                fld_seq_seg_lst.append(fld_seq_lst[i])
        
        fld_data_seg = np.array(fld_data_seg_lst)
        fld_seq_seg = np.array(fld_seq_seg_lst)
        if disorder_cases[0] in c:
            # N-IDR, remove last resiude of fragment
            # and remove the first residue of folded-domain
            nidr = Structure(idp_lst[chain][0])
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

            for i, seq in enumerate(fld_seq_seg):
                next = fld_seq_seg[i + 1]
                fld_data_seg = fld_data_seg[1:]
                if next != seq:
                    break

            idp_lst[chain].pop(0)
            r.pop(0)
            new_struc_arr = np.array(nidr_data_array.tolist() + fld_data_seg.tolist())  # noqa: E501

        if disorder_cases[2] in c:
            # C-IDR, remove last resiude of folded protein
            # and remove the first residue of C-IDR
            cidr = Structure(idp_lst[chain][-1])
            cidr.build()
            cidr_seq = cidr.data_array[:, col_resSeq]
            cidr_data_array = cidr.data_array

            for i, seq in enumerate(cidr_seq):
                next = cidr_seq[i + 1]
                cidr_data_array = cidr_data_array[1:]
                if next != seq:
                    break

            for i, _ in enumerate(fld_seq_seg):
                j = len(fld_seq_seg) - 1 - i
                curr = fld_seq_seg[j]
                prev = fld_seq_seg[j - 1]
                fld_data_seg = fld_data_seg[:-1]
                if prev != curr:
                    break

            # Fix residue connectivity issue
            last_residue_fld = int(fld_data_seg[:, col_resSeq][-1])
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
                
            if disorder_cases[0] in c:
                new_struc_arr = np.array(nidr_data_array.tolist() + fld_data_seg.tolist())  # noqa: E501
                new_struc_arr = np.array(new_struc_arr.tolist() + cidr_data_array.tolist())  # noqa: E501
            else:
                new_struc_arr = np.array(fld_data_seg.tolist() + cidr_data_array.tolist())  # noqa: E501
            idp_lst[chain].pop(-1)
            r.pop(-1)

        if disorder_cases[1] in c:
            if len(c) == 1:
                new_struc_arr = fld_data_seg

            new_struc_seq = new_struc_arr[:, col_resSeq].astype(int)
            first_struc_seq = new_struc_seq[0]

            for idx, minmax in enumerate(r):
                lower = minmax[0]
                upper = minmax[1]
                new_struc_lst = new_struc_arr.tolist()

                try:
                    idr = Structure(idp_lst[chain][0][idx])
                except TypeError:
                    # For cases where there's only 1 L-IDR
                    idr = Structure(idp_lst[chain][0])
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
                if disorder_cases[0] not in c:
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

        new_struc_arr[:, col_chainID] = chain
        new_struc_arr[:, col_segid] = chain
        if len(new_struc_arr_complete) == 0:
            new_struc_arr_complete = new_struc_arr
        else:
            new_struc_arr_complete = np.array(new_struc_arr_complete.tolist() + new_struc_arr.tolist())  # noqa: E501
    if skipped_chains is not None:
        new_struc_arr_complete = np.array(new_struc_arr_complete.tolist() + skipped_chains)  # noqa: E501
    new_serial = [str(i) for i in range(1, len(new_struc_arr_complete) + 1)]
    new_struc_arr_complete[:, col_serial] = new_serial
    
    return new_struc_arr_complete


def inter_chain_cc(all_files, nconfs):
    """
    Check for clashes between all chains of IDRs.

    Parameters
    ----------
    all_files : dict
        Keys are chains and values are lists of Paths
        for IDRs selected that don't clash with each other.
    
    nconfs : int
        Number of final conformers.
    
    Returns
    -------
    final_combos : list of dict
        Each value is a dictionary where key is the chain and value
        are the list containing IDRs.
    """
    final_combos = []
    total_chains = len(all_files)
    all_chains = all_files.keys()

    while len(final_combos) < nconfs:
        reset = False
        random_combo = [random.randint(0, nconfs - 1) for _ in range(total_chains)]  # noqa: E501
        candidate_combo = []
        
        for i, selections in enumerate(all_files.values()):
            sele = selections[random_combo[i]]
            if type(sele) is list:
                candidate_combo.append(sele)
            else:
                candidate_combo.append([sele])
            
        for i in range(len(candidate_combo)):
            for j in range(i + 1, len(candidate_combo)):
                idr_set1 = candidate_combo[i]
                idr_set2 = candidate_combo[j]
                idr_pairs = product(idr_set1, idr_set2)
                for pair in idr_pairs:
                    frag_s = Structure(pair[0])
                    struc_s = Structure(pair[1])
                    frag_s.build()
                    struc_s.build()
                    clashes, _ = count_clashes(frag_s.data_array, struc_s)
                    if clashes:
                        reset = True
                        break
                if reset:
                    break
            if reset:
                break
        if reset:
            continue
        else:
            temp_dict = {}
            for i, c in enumerate(all_chains):
                temp_dict[c] = candidate_combo[i]
            final_combos.append(temp_dict)

    return final_combos
