"""
Client for building IDRs on PDB files in the cartesian coordinate space.

Methodology deviates from traditional IDP or beads-on-a-string FLDR/S approach.

Name: FLDR/S (Folded disordered region/structure sampling)
"""
import math
import os
import random
from functools import partial
from itertools import combinations, islice, product
from types import SimpleNamespace

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core.build_definitions import (
    build_bend_CA_C_O,
    build_bend_H_N_C,
    distance_C_O,
    distance_H_N,
    )
from idpconfgen.core.definitions import (
    aa3to1,
    vdW_radii_ionic_CRC82,
    vdW_radii_tsai_1999,
    )
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs.libcalc import calc_torsion_angles, make_coord_Q_planar
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


def _combination_fragments(selected):
    """
    Return the IDR fragment paths of a combination as a flat list.

    Combinations nest the Linker-IDRs of a chain in their own list,
    for example ``[N-IDR, [Linker-IDR 0, Linker-IDR 1], C-IDR]``.

    Parameters
    ----------
    selected : Path, list or tuple
        A fragment path or a (nested) list or tuple of them.

    Returns
    -------
    list of Path
    """
    if isinstance(selected, (list, tuple)):
        fragments = []
        for item in selected:
            fragments.extend(_combination_fragments(item))
        return fragments
    return [selected]


def _idr_rows(rows, case=None):
    """
    Return the rows of the IDR residues of an LDRS fragment.

    LDRS fragments carry two folded residues bonded to the folded
    domain: after an N-IDR, before a C-IDR, and on both sides of a
    Linker-IDR. Those are removed.

    Parameters
    ----------
    rows : np.ndarray
        Data array of the fragment.

    case : str, optional
        Disorder case of the fragment. None keeps every residue.

    Returns
    -------
    np.ndarray
    """
    resseq = rows[:, col_resSeq].astype(int)
    keep = np.unique(resseq)
    if case in (disorder_cases[1], disorder_cases[2]):
        keep = keep[2:]
    if case in (disorder_cases[0], disorder_cases[1]):
        keep = keep[:-2]
    return rows[np.isin(resseq, keep)]


def combinations_clash_check(selected, cases=None):
    """
    Clash-check IDRs within a combination against each other.

    Every IDR fragment of the combination is checked against every
    other one, whatever the mix of N-IDR, Linker-IDRs and C-IDR.

    Parameters
    ----------
    selected : list or tuple
        One combination: one item from each list given to
        :func:`create_combinations`, for example
        ``(N-IDR, [Linker-IDR 0, Linker-IDR 1], C-IDR)``.

    cases : list of str, optional
        Disorder case of each item of ``selected``, for example
        ``['N-IDR', 'Linker-IDR', 'C-IDR']``. When given, the folded
        residues carried by the fragments are not checked: IDRs a few
        folded residues apart carry the same folded residues.

    Returns
    -------
    list or tuple, or False
        ``selected`` when no pair of IDRs has more than 40 contacts
        (within the vdW radii sum + 1.0 A), False otherwise.
    """
    if cases is None:
        cases = [None] * len(selected)

    idrs = []
    for item, case in zip(selected, cases):
        for path in _combination_fragments(item):
            structure = Structure(path)
            structure.build()
            rows = _idr_rows(structure.data_array, case)
            if len(rows):
                idrs.append(rows)

    for first, second in combinations(idrs, 2):
        num_clashes, _ = count_clashes(
            first,
            SimpleNamespace(data_array=second),
            tolerance=1.0,
            )
        if type(num_clashes) is bool:
            return False

    return selected


def _random_combinations(lst):
    """
    Yield the combinations of one item from each list in random order.

    Each combination is yielded once. Large products are sampled
    index by index instead of being built in memory.

    Parameters
    ----------
    lst : list of lists

    Yields
    ------
    tuple
        One item from each list, in the order of ``lst``.
    """
    sizes = [len(items) for items in lst]
    total = math.prod(sizes)
    if total <= 100_000:
        all_combinations = list(product(*lst))
        random.shuffle(all_combinations)
        yield from all_combinations
        return

    seen = set()
    while len(seen) < total:
        indexes = tuple(random.randrange(size) for size in sizes)
        if indexes in seen:
            continue
        seen.add(indexes)
        yield tuple(items[i] for items, i in zip(lst, indexes))


def create_combinations(
        lst,
        num_combinations,
        ncores=1,
        max_trials=None,
        cases=None,
        ):
    """
    Create unique combinations between list of lists.
    
    Made for N-IDR and C-IDR combinations. Where list1 = N-IDR paths,
    and list2 = C-IDR paths. Itertools product is used here because
    order matters.

    Combinations are drawn at random without repetition and clash-checked
    with :func:`combinations_clash_check` until ``num_combinations`` pass,
    every combination was tried, or ``max_trials`` were tried.
    
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

    max_trials : int, optional
        Maximum number of combinations to clash-check.
        Defaults to ``max(100, 20 * num_combinations)``.

    cases : list of str, optional
        Disorder case of each list in ``lst``, given to
        :func:`combinations_clash_check`.
    
    Return
    ------
    selected_combinations : list
        List of lists of different combinations as follows:
        [[item from list 1, item from list 2], ...]
        Fewer than ``num_combinations`` if fewer passed (a warning
        is logged).

    Raises
    ------
    IDPConfGenException
        If no combination passed the clash check.
    """
    candidates = _random_combinations(lst)

    if len(lst) == 1:
        return convert_tuples_to_lists(islice(candidates, num_combinations))

    if max_trials is None:
        max_trials = max(100, 20 * num_combinations)

    check = partial(combinations_clash_check, cases=cases)
    passed = []
    tried = 0
    while len(passed) < num_combinations and tried < max_trials:
        batch_size = min(
            max(num_combinations - len(passed), ncores),
            max_trials - tried,
            )
        batch = list(islice(candidates, batch_size))
        if not batch:
            break
        tried += len(batch)
        execute_clash = pool_function(check, batch, ncores=ncores)
        for result in execute_clash:
            if result is not False:
                passed.append(result)

    if not passed:
        raise IDPConfGenException(errmsg=(
            f"None of the {tried} combinations of IDR conformers tried "
            "passed the clash check between IDRs (more than 40 contacts "
            "within the vdW radii sum + 1.0 A). Build more conformers "
            "for each IDR."
            ))

    if len(passed) < num_combinations:
        log.warning(
            f"Only {len(passed)} of the {num_combinations} requested "
            f"combinations of IDR conformers passed the clash check "
            f"between IDRs ({tried} tried)."
            )

    return convert_tuples_to_lists(passed[:num_combinations])


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
    nidr, lidr, cidr = disorder_cases[0], disorder_cases[1], disorder_cases[2]  # noqa: E501
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
            # keep Linker-IDRs in sequence order: linker_0, linker_1, ...
            lidr_cases_dir = sorted(
                os.listdir(lidr_path),
                key=lambda name: int(name.rsplit('_', 1)[-1]),
                )
            lidr_confs_lst = []
            for cpath in lidr_cases_dir:
                num = cpath.rsplit('_', 1)[-1]
                lidr_matches = lidr_path.joinpath(cpath, f"{num}_match")
                lidr_confs = os.listdir(lidr_matches)
                lidr_files = [Path(lidr_matches.joinpath(fpath)) for fpath in lidr_confs]  # noqa: E501
                lidr_confs_lst.append(lidr_files)
            lidr_combinations = create_combinations(
                lidr_confs_lst,
                nconfs,
                ncores,
                cases=[lidr] * len(lidr_confs_lst),
                )
        if os.path.isdir(cidr_path):
            cidr_confs = os.listdir(cidr_path)
            cidr_files = [Path(cidr_path.joinpath(cpath)) for cpath in cidr_confs]  # noqa: E501

        if len(nidr_files) and len(cidr_files) and len(lidr_combinations) > 0:
            combinations[c] = create_combinations([nidr_files, lidr_combinations, cidr_files], nconfs, ncores, cases=[nidr, lidr, cidr])  # noqa: E501
        elif len(nidr_files) and len(lidr_combinations) > 0:
            combinations[c] = create_combinations([nidr_files, lidr_combinations], nconfs, ncores, cases=[nidr, lidr])  # noqa: E501
        elif len(cidr_files) and len(lidr_combinations) > 0:
            combinations[c] = create_combinations([lidr_combinations, cidr_files], nconfs, ncores, cases=[lidr, cidr])  # noqa: E501
        elif len(nidr_files) and len(cidr_files) > 0:
            combinations[c] = create_combinations([nidr_files, cidr_files], nconfs, ncores, cases=[nidr, cidr])  # noqa: E501
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

    The C-side conformer (from the folded region before the break)
    and each N-side conformer (from the folded region after it)
    cover the same residues. Wherever the C of a residue of the first and
    the N and CA of the next residue of the second form a peptide bond,
    the two are joined there: residues up to the junction come from the
    C-side conformer and the others from the N-side one. A pair of
    conformers can meet at several residues; every junction that passes
    the clash check is written, as ``<C-side>+<N-side>_j<residue>.pdb``.

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
    int
        Number of joined conformers written to the output folder.
    """
    matches = 0
    idr_struc = Structure(Path(cterm_idr))
    idr_struc.build()
    idr_arr = idr_struc.data_array
    idr_name = idr_arr[:, col_name]
    idr_coords = idr_arr[:, cols_coords].astype(float)
    idr_resseq = idr_arr[:, col_resSeq].astype(int)
    cterm_idr_stem = Path(cterm_idr).stem
    
    idr_C = []
    idr_CA = []
    idr_res = []
    for n, name in enumerate(idr_name):
        if name == 'C':
            idr_C.append(idr_coords[n])
            idr_res.append(idr_resseq[n])
        elif name == 'CA':
            idr_CA.append(idr_coords[n])
    
    for nterm_idr in nterm_idr_lib:
        nterm_idr_struc = Structure(Path(nterm_idr))
        nterm_idr_struc.build()
        nterm_idr_full = nterm_idr_struc.data_array
        nterm_idr_name = nterm_idr_full[:, col_name]
        nterm_idr_coords = nterm_idr_full[:, cols_coords].astype(float)
        nterm_idr_stem = Path(nterm_idr).stem
        
        nterm_idr_N = []
        nterm_idr_CA = []
        for n, name in enumerate(nterm_idr_name):
            if name == 'N':
                nterm_idr_N.append(nterm_idr_coords[n])
            elif name == 'CA':
                nterm_idr_CA.append(nterm_idr_coords[n])
        
        # every residue where the pair meets gives a joined conformer
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
                
                # Both conformers share one atom layout: join the C-side
                # conformer up to `term_residue` with the N-side one after
                # it. Slice the full arrays at every junction; they serve
                # the other junctions and the other N-side conformers.
                cut = int(np.argmax(idr_resseq > term_residue))
                idr_list = list(idr_arr[:cut])
                nterm_idr_arr = nterm_idr_full[cut:]

                nterm_idr_struc._data_array = nterm_idr_arr
                
                clashes, _ = count_clashes(
                    np.array(idr_list),
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
                    
                    # Place the carbonyl O and the amide H of the new
                    # peptide bond in its plane
                    new_O_xyz = make_coord_Q_planar(
                        idr_CA[i],
                        curr_c,
                        next_n,
                        distance=distance_C_O,
                        bend=build_bend_CA_C_O,
                        )
                    for col, value in zip(cols_coords, new_O_xyz):
                        final_struc_arr[O_idx, col] = f'{value:.3f}'
                    
                    if H_idx >= 0:
                        new_H_xyz = make_coord_Q_planar(
                            next_ca,
                            next_n,
                            curr_c,
                            distance=distance_H_N,
                            bend=build_bend_H_N_C,
                            )
                        for col, value in zip(cols_coords, new_H_xyz):
                            final_struc_arr[H_idx, col] = f'{value:.3f}'
                    
                    final_struc = structure_to_pdb(final_struc_arr)
                    matches += 1
                    write_PDB(final_struc, str(output_folder) + f"/{cterm_idr_stem}+{nterm_idr_stem}_j{term_residue}.pdb")  # noqa: E501

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

            nidr_resseq = nidr_data_array[:, col_resSeq].astype(int)
            offset = int(fld_seq_seg[0]) - nidr_resseq[-1]
            nidr_data_array[:, col_resSeq] = (nidr_resseq + offset).astype(str)  # noqa: E501

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
            if disorder_cases[0] not in c and disorder_cases[2] not in c:
                # Linker-IDRs only
                new_struc_arr = fld_data_seg

            new_struc_seq = new_struc_arr[:, col_resSeq].astype(int)
            first_struc_seq = new_struc_seq[0]

            # The Linker-IDRs come nested, ``[[L0, L1, ...]]``, or,
            # for chains with Linker-IDRs only, flat, ``[L0, L1, ...]``
            linker_files = idp_lst[chain][0]
            if not isinstance(linker_files, (list, tuple)):
                linker_files = idp_lst[chain]

            for idx, minmax in enumerate(r):
                lower = minmax[0]
                upper = minmax[1]
                new_struc_lst = new_struc_arr.tolist()

                idr = Structure(linker_files[idx])
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
                # residues on either side of the chain break; the first
                # residue of the chain is position 0 of the sequence
                actual_lower = lower + first_struc_seq - 1
                actual_upper = upper + first_struc_seq
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


def inter_chain_cc(all_files, nconfs, max_trials=None):
    """
    Check for clashes between all chains of IDRs.

    Combinations of one entry per chain are drawn at random without
    repetition until ``nconfs`` pass, every combination was tried, or
    ``max_trials`` were tried. A combination passes when no IDR of a
    chain has a contact (within the vdW radii sum + 0.4 A) with an IDR
    of another chain.

    Parameters
    ----------
    all_files : dict
        Keys are chains and values are lists of Paths
        for IDRs selected that don't clash with each other.
    
    nconfs : int
        Number of final conformers.

    max_trials : int, optional
        Maximum number of combinations across chains to check.
        Defaults to ``max(1000, 100 * nconfs)``.
    
    Returns
    -------
    final_combos : list of dict
        Each value is a dictionary where key is the chain and value
        are the list containing IDRs.
        Fewer than ``nconfs`` if fewer passed (a warning is logged).

    Raises
    ------
    IDPConfGenException
        If no combination passed the clash check.
    """
    if max_trials is None:
        max_trials = max(1000, 100 * nconfs)

    chains = list(all_files)
    selections = []
    for chain in chains:
        selections.append([
            sele if type(sele) is list else [sele]
            for sele in all_files[chain]
            ])

    structures = {}

    def _built(path):
        if path not in structures:
            structure = Structure(path)
            structure.build()
            structures[path] = structure
        return structures[path]

    final_combos = []
    tried = 0
    for candidate_combo in _random_combinations(selections):
        if len(final_combos) == nconfs or tried == max_trials:
            break
        tried += 1

        clash = False
        for idr_set1, idr_set2 in combinations(candidate_combo, 2):
            idr_pairs = product(
                _combination_fragments(idr_set1),
                _combination_fragments(idr_set2),
                )
            for path1, path2 in idr_pairs:
                clashes, _ = count_clashes(
                    _built(path1).data_array,
                    _built(path2),
                    )
                if clashes:
                    clash = True
                    break
            if clash:
                break

        if not clash:
            final_combos.append(dict(zip(chains, candidate_combo)))

    if not final_combos:
        raise IDPConfGenException(errmsg=(
            f"None of the {tried} combinations of IDR conformers across "
            "chains tried passed the clash check between chains. Build "
            "more conformers for each IDR."
            ))

    if len(final_combos) < nconfs:
        log.warning(
            f"Only {len(final_combos)} of the {nconfs} requested "
            "combinations of IDR conformers passed the clash check "
            f"between chains ({tried} tried)."
            )

    return final_combos


def atom_vdw_radii(elements, names=None, default=1.7):
    """
    Return the van der Waals radius of each atom as in `count_clashes`.

    Parameters
    ----------
    elements : array-like of str
        Element symbols. Missing or unknown symbols are inferred from the
        first letter of `names`.

    names : array-like of str, optional
        Atom names.

    default : float
        Radius used when the element cannot be resolved.

    Returns
    -------
    np.ndarray of float
    """
    all_vdw_radii = {**vdW_radii_tsai_1999, **vdW_radii_ionic_CRC82}
    elements = [str(e).strip().upper() for e in elements]
    if names is None:
        names = [''] * len(elements)
    radii = np.empty(len(elements), dtype=np.float64)
    for i, (element, name) in enumerate(zip(elements, names)):
        if element not in all_vdw_radii:
            element = str(name).strip().lstrip('0123456789')[:1].upper()
        radii[i] = all_vdw_radii.get(element, default)
    return radii


def is_hydrogen(elements, names):
    """Return a boolean mask of hydrogen atoms."""
    out = []
    for element, name in zip(elements, names):
        element = str(element).strip().upper()
        if not element:
            element = str(name).strip().lstrip('0123456789')[:1].upper()
        out.append(element in ('H', 'D'))
    return np.array(out, dtype=bool)


def backbone_anchor(rows):
    """
    Return N, CA, C, O coordinates of a residue.

    Parameters
    ----------
    rows : np.ndarray
        Structure data array rows of one residue.

    Returns
    -------
    np.ndarray, shape (4, 3)
        NaN for missing atoms.
    """
    anchor = np.full((4, 3), np.nan)
    names = list(rows[:, col_name])
    for k, atom in enumerate(('N', 'CA', 'C', 'O')):
        if atom in names:
            anchor[k] = rows[names.index(atom), cols_coords].astype(float)
    return anchor


PACKING_CODES = {
    **{key: value for key, value in aa3to1.items() if value.isupper()},
    'HIP': 'H',
    'HIE': 'H',
    'HID': 'H',
    'HSD': 'H',
    'HSE': 'H',
    'HSP': 'H',
    'SEP': 'S',
    'TPO': 'T',
    'CYX': 'C',
    'MSE': 'M',
    }


def packing_context(rows, exclude=None):
    """
    Collect template residues to pack side chains next to.

    Parameters
    ----------
    rows : np.ndarray
        Structure data array; residues are consecutive rows with the same
        chain, segment, number and name.

    exclude : np.ndarray of bool, optional
        Rows of residues to leave out.

    Returns
    -------
    dict
        ``sequence`` (one-letter codes of the protein residues with
        complete N, CA, C and O), ``backbone`` (their N, CA, C, O,
        shape (R, 4, 3)), ``atoms`` (their heavy atoms) and ``owner``
        (the residue of each atom).
    """
    exclude = (
        np.zeros(len(rows), dtype=bool) if exclude is None
        else np.asarray(exclude, dtype=bool)
        )
    keys = [
        (c, s, r, n) for c, s, r, n in zip(
            rows[:, col_chainID],
            rows[:, col_segid],
            rows[:, col_resSeq],
            rows[:, col_resName],
            )
        ]
    hydrogen = is_hydrogen(rows[:, col_element], rows[:, col_name])
    sequence, backbone, atoms, owner = [], [], [], []
    start = 0
    for stop in range(1, len(rows) + 1):
        if stop < len(rows) and keys[stop] == keys[start]:
            continue
        code = PACKING_CODES.get(str(keys[start][3]).strip())
        part = slice(start, stop)
        start = stop
        if code is None or exclude[part].any():
            continue
        anchor = backbone_anchor(rows[part])
        if not np.all(np.isfinite(anchor)):
            continue
        heavy = ~hydrogen[part]
        sequence.append(code)
        backbone.append(anchor)
        atoms.append(rows[part][heavy][:, cols_coords].astype(np.float64))
        owner.append(np.full(int(heavy.sum()), len(owner)))
    return {
        'sequence': ''.join(sequence),
        'backbone': np.array(backbone).reshape(-1, 4, 3),
        'atoms': np.concatenate(atoms) if atoms else np.zeros((0, 3)),
        'owner': np.concatenate(owner) if owner else np.zeros(0, dtype=int),
        }


def prepare_ccd_context(chain_rows, environment_rows, lower, upper):
    """
    Collect what CCD needs to close a Linker-IDR onto its folded anchors.

    Parameters
    ----------
    chain_rows : np.ndarray
        Structure data array of the protein atoms of the modeled chain.

    environment_rows : np.ndarray
        Structure data array of everything the loop must not clash with:
        all chains, cofactors and membrane of the template.

    lower, upper : int
        Residue numbers of the last folded residue before the IDR and of
        the first folded residue after it.

    Returns
    -------
    dict
        ``lower_anchor`` and ``upper_anchor`` (N, CA, C, O coordinates),
        ``lower_outer`` and ``upper_outer`` (the same for the residues
        before `lower` and after `upper`, NaN if missing), ``lower_rows``
        and ``upper_rows`` (template rows of the residues bonded to the
        IDR), the environment arrays ``env_coords``, ``env_radii``,
        ``env_hydrogen``, ``env_lower`` and ``env_upper`` (the last two
        mark the atoms of the residues bonded to the IDR), and
        ``packing`` (:func:`packing_context` of the environment without
        the four overlap residues).
    """
    resseq = chain_rows[:, col_resSeq].astype(int)
    lower_rows = chain_rows[resseq == int(lower)]
    upper_rows = chain_rows[resseq == int(upper)]
    lower_anchor = backbone_anchor(lower_rows)
    upper_anchor = backbone_anchor(upper_rows)
    if not (
            np.all(np.isfinite(lower_anchor[:3]))
            and np.all(np.isfinite(upper_anchor[:3]))
            ):
        raise IDPConfGenException(errmsg=(
            'Loop closure needs N, CA and C of residues '
            f'{lower} and {upper} in the folded structure.'
            ))
    overlap_rows = chain_rows[np.isin(
        resseq,
        (int(lower) - 1, int(lower), int(upper), int(upper) + 1),
        )]

    env_names = environment_rows[:, col_name]
    env_elements = environment_rows[:, col_element]
    keys = [
        (str(n).strip(), str(x).strip(), str(y).strip(), str(z).strip())
        for n, x, y, z in zip(
            env_names,
            environment_rows[:, col_x],
            environment_rows[:, col_y],
            environment_rows[:, col_z],
            )
        ]

    def _mask(rows):
        wanted = {
            (str(n).strip(), str(x).strip(), str(y).strip(), str(z).strip())
            for n, x, y, z in zip(
                rows[:, col_name],
                rows[:, col_x],
                rows[:, col_y],
                rows[:, col_z],
                )
            }
        return np.array([k in wanted for k in keys], dtype=bool)

    return {
        'lower_anchor': lower_anchor,
        'upper_anchor': upper_anchor,
        'lower_outer': backbone_anchor(chain_rows[resseq == int(lower) - 1]),
        'upper_outer': backbone_anchor(chain_rows[resseq == int(upper) + 1]),
        'lower_rows': np.array(lower_rows, copy=True),
        'upper_rows': np.array(upper_rows, copy=True),
        'env_coords': environment_rows[:, cols_coords].astype(np.float64),
        'env_radii': atom_vdw_radii(env_elements, env_names),
        'env_hydrogen': is_hydrogen(env_elements, env_names),
        'env_lower': _mask(lower_rows),
        'env_upper': _mask(upper_rows),
        'packing': packing_context(environment_rows, _mask(overlap_rows)),
        }


def splice_folded_rows(fragment_rows, lower_rows, upper_rows):
    """
    Replace the junction residues of a closed loop with the folded ones.

    Parameters
    ----------
    fragment_rows : np.ndarray
        Structure data array of the closed loop, residues numbered 1 to L.

    lower_rows, upper_rows : np.ndarray
        Template rows of the folded residues before and after the IDR.

    Returns
    -------
    np.ndarray
    """
    resseq = fragment_rows[:, col_resSeq].astype(int)
    last = int(resseq.max())

    def _renumbered(rows, number):
        rows = np.array(rows, copy=True)
        rows[:, col_resSeq] = str(number)
        rows[:, col_chainID] = fragment_rows[0, col_chainID]
        rows[:, col_segid] = ''
        return rows

    parts = [
        fragment_rows[resseq == 1],
        _renumbered(lower_rows, 2),
        fragment_rows[(resseq > 2) & (resseq < last - 1)],
        _renumbered(upper_rows, last - 1),
        fragment_rows[resseq == last],
        ]
    out = np.concatenate(parts, axis=0)
    out[:, col_serial] = [str(i) for i in range(1, len(out) + 1)]
    return out
