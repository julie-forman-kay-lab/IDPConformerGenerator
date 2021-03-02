"""
Tools to validate conformers.

Recognition of this module should be grated to:
    - @AlaaShamandy
    - @joaomcteixeira

    Please see:
    https://github.com/julie-forman-kay-lab/IDPConformerGenerator/pull/23
"""
import numpy as np
from scipy.spatial import distance

from idpconfgen import Path
from idpconfgen.core.build_definitions import (
    average_distance_CA_C,
    average_distance_C_Np1,
    average_distance_N_CA,
    )
from idpconfgen.core.definitions import vdW_radii_dict
from idpconfgen.libs.libstructure import (
    Structure,
    col_element,
    col_name,
    col_resSeq,
    cols_coords,
    generate_backbone_pairs_labels,
    )


# Notes on improving distance calculation and vdW energy function
# https://stackoverflow.com/questions/37298710/fastest-way-to-compute-distance-beetween-each-points-in-python
# https://stackoverflow.com/questions/37794849/efficient-and-precise-calculation-of-the-euclidean-distance


def vdw_clash_by_threshold_common_preparation(
        protein_atoms,
        protein_elements,
        residue_numbers,
        atoms_to_consider=False,
        elements_to_consider=False,
        residues_apart=2,
        vdW_radii='tsai1999',
        ):
    """
    Prepare masks for vdW clash calculation.

    Masks are prepared considering all-to-all distances will be computed
    using scipy.distance.cdist, so a (N, 3) array originates a (N, N)
    distance result.

    `atoms_to_consider` and `elements_to_consider` are evaluated with
    logical AND, that is, only entries that satisfy both are considered.

    .. see-also::
        vdw_clash_by_threshold_calc

    Parameters
    ----------
    protein_atoms : np.array, of shape (N,), dtype string
        A sequence of the protein atom names.

    protein_elements : np.array, of shape (N,), dytpe compatible with '<U2'
        A sequence of the protein atom elements aligned with `protein_atoms`.

    residue_numbers : np.array of shape (N,), dtype=np.int

    atoms_to_consider : sequence, optional
        A tuple of the atom names to consider in the calculation.
        Defaults to FALSE, considers all atoms.

    elements_to_consider : sequence, optional
        A tuple of the element types to consider in the calculation.
        Defaults to FALSE, considers all elements.
    """
    # uses all options if values evaluate to FALSE
    atoms_to_consider = atoms_to_consider or protein_atoms
    elements_to_consider = elements_to_consider or protein_elements

    atoms_to_analyze = np.logical_and(
        np.isin(protein_elements, elements_to_consider),
        np.isin(protein_atoms, atoms_to_consider),
        )

    # computes only pairs of atoms where both atoms are in atoms_to_consider
    atc_mask = np.logical_and.outer(atoms_to_analyze, atoms_to_analyze)

    vdW_vector = np.zeros(protein_atoms.size)
    for atom, radius in vdW_radii_dict[vdW_radii].items():
        vdW_vector[protein_elements == atom] = radius

    pure_radii_sum = np.add.outer(vdW_vector, vdW_vector)

    residue_distances = np.subtract.outer(residue_numbers, residue_numbers)

    return atc_mask, pure_radii_sum, residue_distances >= residues_apart


def vdw_clash_by_threshold_calc(
        coords,
        atc_mask,
        pure_radii_sum,
        distances_apart,
        vdW_overlap=0.0,
        ):
    """
    Calculate van der Waals clashes from a pure sphere overlap.

    Other masks used as parameters will be applied to the result of:

    scipy.distance.cdist(coords, coords, 'euclidean')

    Parameters
    ----------
    coords : np.array, dtype=np.float, shape (N, 3)
        The protein XYZ coordinates.

    atc_mask : np.array, dtype=np.bool, shape (N, N)
        A boolean masks to filter only the atoms relevant to report.
        Usually this mask is prepared beforehand and can contain different
        considerations, such as residues apart and especific atom types.
        If `coords` contain only the coordinates desired to compute,
        then `atc_mask` should contain only TRUE entries.

    pure_radii_sum : np.array, dtype=float, shape (N, N)
        An all-to-all sum of the vdW radii. In other words, the threeshold
        after which a clash is considered to exist for each atom pair,
        before applying `vdW_overlap` allowance.

    distances_apart : np.array, dtype=int, shape (N, N)
        An all-to-all atom-to-atom residue to residue distance matrix.

    vdW_overlap : float
        The vdW overlap tolerance to apply.

    Returns
    -------
    tuple
        rows, cols : of cdist applied to `coords` where clashes where found.
        distances found for those clashes
        computed distance threshold
        overlap, computed overlap distance between threshold and distance
    """
    distances = distance.cdist(coords, coords, 'euclidean')

    # computes clahses based on overlap, so we can apply a overlap
    # threshold
    overlap = pure_radii_sum - distances
    clashes_raw = overlap >= vdW_overlap

    rows_cols = np.logical_and(
        np.logical_and(clashes_raw, atc_mask),
        distances_apart,
        ).nonzero()

    rows, cols = rows_cols

    return (
        rows,
        cols,
        distances[rows_cols],
        pure_radii_sum[rows_cols],
        overlap[rows_cols],
        )


def vdw_clash_by_threshold(
        coords,
        protein_atoms,
        protein_elements,
        atoms_to_consider,
        elements_to_consider,
        residue_numbers,
        residues_apart=2,
        vdW_radii='tsai1999',
        vdW_overlap=0.0,
        # 07/ago/2020, I have decided not to use **kwargs here
        # even if that requires a duplication of parameters in this function
        ):
    """
    Calculate vdW clashes from XYZ coordinates and identity masks.

    Parameters
    ----------
    coordinates : numpy array, dtype=float, shape (N, 3)
        The atom XYZ coordinates.

    protein_atoms : numpy array, dtype=str, shape (N,)
        The protein atom names.

    protein_elements : numpy array, dtype=str, shape(N,)
        The protein atom elements.

    atoms_to_consider : list-like
        The atoms in `protein_atoms` to consider in the vdW clash
        analysis.

    elements_to_consider : list-like
        The elements in `protein_elements` to consider in the vdW clash
        analysis.

    residue_number : numpy array, dtype=int, shape (N,)
        The residue number corresponding to each atom.

    residues_apart : int, optional
        The minimum number of residues apart to consider for a clash.
        Defaults to 2.

    vdW_radii : str, optional
        The VDW radii set to consider. Defaults to 'tsai1999'.


    vdW_overlap : float, optional
        An overlap allowance in Angstroms.
        Defaults to 0.0, any distance less than vdW+vdW is considered
        a clash.

    Returns
    -------
    Same :func:`vdw_clash_by_threshold_calc` returns.
    """
    # The VDW clash analysis is separated in two steps,
    # 1) the preparation of the data and masks
    # 2) the actual clash calculation

    atc_mask, pure_radii_sum, distances_apart = \
        vdw_clash_by_threshold_common_preparation(
            protein_atoms,
            protein_elements,
            residue_numbers,
            atoms_to_consider=atoms_to_consider,
            elements_to_consider=elements_to_consider,
            residues_apart=residues_apart,
            vdW_radii=vdW_radii,
            )

    return vdw_clash_by_threshold_calc(
        coords,
        atc_mask,
        pure_radii_sum,
        distances_apart,
        vdW_overlap=vdW_overlap,
        )


def evaluate_vdw_clash_by_threshold_from_disk(
        name,
        pdb_data,
        atoms_to_consider,
        elements_to_consider,
        **kwargs,
        ):
    """
    Evaluate clashes in a structure.

    Created to evaluate conformers from disk.

    .. see-also::
        validate_vdw_clashes
    """
    s = Structure(pdb_data)
    s.build()
    da = s.data_array
    atom_names = da[:, col_name]
    atom_elements_pure = da[:, col_element]
    atom_elements_from_names = atom_names.astype('<U1')

    # elements_dont_exist = np.logical_not(atom_elements_pure.astype(np.bool))
    #
    # ValueError: invalid literal for int() with base 10:
    # really? okay, lets do a list comprehension...
    elements_dont_exist = np.array([not bool(i) for i in atom_elements_pure])

    # considers PDBs that lab element information
    best_observation_ele = np.zeros(atom_elements_pure.shape, dtype='<U2')
    best_observation_ele[:] = atom_elements_pure[:]
    best_observation_ele[elements_dont_exist] = \
        atom_elements_from_names[elements_dont_exist]

    res_numbers = da[:, col_resSeq].astype(np.int)
    coords = s.coords

    rows, cols, distances, radii_sum, overlap = vdw_clash_by_threshold(
        coords,
        atom_names,  # protein atoms
        best_observation_ele,  # protein elements
        atoms_to_consider,  # atoms to consider in the validation
        elements_to_consider,  # elements to consider in the validation
        res_numbers,
        **kwargs,
        )

    report = report_vdw_clash(da, rows, cols, distances, radii_sum, overlap)

    # rows.size is the number of clashes
    return name, rows.size, report


def validate_bb_bond_len(coords, tolerance=0.01):
    """
    Validate backbone bond lengths of `coords`.

    Considers `coords` are already sorted to (N, CA, C) per residue.
    Considers only (N, CA, C) atoms are present in `coords`.
    Evalutes against N-CA, CA-C and C-Np1 distances used in IDPConfGen.

    Parameters
    ----------
    tolerance : float
        A tolerance in the same units as `coords`.
        Conflicts under the tolerance are consider valid.

    Returns
    -------
    np.array, dtype=bool, shape (N-1,)
        True if bond length is invalid.
        False if bond length is valid.


    np.array, dtype=np.float, shape (N-1,)
        The computed bond distances.

    np.array, dtype=np.float, shape (N-1,)
        The expected bond lengths
    """
    expected_bond_length = np.tile(
        [
            average_distance_N_CA,
            average_distance_CA_C,
            average_distance_C_Np1
            ],
        coords.shape[0] // 3
        )

    # https://stackoverflow.com/questions/1401712
    # bond_distances = np.linalg.norm(coords[:-1] - coords[1:], axis=1)
    _coords = coords[:-1] - coords[1:]
    bond_distances = np.sqrt(np.einsum("ij,ij->i", _coords, _coords))

    # true when bond length is too different
    invalid = np.invert(np.isclose(
        bond_distances,
        expected_bond_length[:-1],
        atol=tolerance,
        ))

    return invalid, bond_distances, expected_bond_length[:-1]


def validate_bb_bonds_len_from_disk(
        name=None,
        pdb_data=None,
        tolerance=0.1,
        ):
    """
    Validate backbone bond lengths of a structure stored in disk.
    """
    # read structure
    pdb = Path(name) if name else pdb_data
    s = Structure(pdb)
    s.build()
    s.add_filter_backbone(minimal=True)
    fa = s.filtered_atoms
    coords = s.get_sorted_minimal_backbone_coords(filtered=False)

    invalid, bond_distances, expected_bond_length = \
        validate_bb_bond_len(coords, tolerance)

    labels = generate_backbone_pairs_labels(fa)

    report = report_sequential_bon_len(
        bond_distances,
        expected_bond_length,
        invalid,
        labels,
        )

    return name, np.sum(invalid), report


def eval_bb_bond_length_distribution(name, pdb_data):
    """."""
    s = Structure(pdb_data)
    s.build()
    coords = s.get_sorted_minimal_backbone_coords()
    _coords = coords[:-1] - coords[1:]
    return np.sqrt(np.einsum('ij,ij->i', _coords, _coords))


def report_vdw_clash(data_array, pair1, pair2, distances, radii_sum, overlap):
    """
    Prepare a report of the identified clashes.

    Parameters
    ----------
    data_array : Numpy ndarray, shape (N, M)
        A numpy data_array as given by
        :attr:`idpconfgen.libs.libstructure.Structure.data_array`.

    pair1, pair2 : ordered iterable of integers
        The row indexes where to retrieve atom information from `data_array´.
        pair1 and pair2 must be aligned in order for the report to make sense,
        that it, the first item of pair1 clashes with the first item of pair2.

    distances, threshold : indexable of length equal to pair1/2 length
        Contain distances and clash thresholds for the different
        identified clashes.

    Returns
    -------
    str
        The report.
    """
    report = []
    for i, (p1, p2) in enumerate(zip(pair1, pair2)):
        n1, r1, a1 = data_array[p1, cols_coords]
        n2, r2, a2 = data_array[p2, cols_coords]
        report.append(
            f"{n1:>5} {r1:>3} {a1:>5} - "
            f"{n2:>5} {r2:>3} {a2:>5} - "
            f'r: {radii_sum[i]:.3f} - '
            f'd: {distances[i]:.3f} - '
            f'o: {overlap[i]:.3f}'
            )

    return '\n'.join(report)


def report_sequential_bon_len(
        bond_distances,
        expected_bond_length,
        invalid_bool,
        labels,
        ):
    """
    Generate a report from bond distances of sequential bonds.

    Returns
    -------
    string
    """
    BD = bond_distances[invalid_bool]
    EB = expected_bond_length[invalid_bool]
    DIFF = np.abs(BD - EB)
    LABELS = labels[invalid_bool]

    rows = (
        f'{l} - f: {bd:.3f} e: {eb:.3f} - {d:.3f}'
        for bd, eb, d, l in zip(BD, EB, DIFF, LABELS)
        )

    return '\n'.join(rows)


def validate_conformer_for_builder(
        coords,
        atom_labels,
        residue_numbers,
        bb_mask,
        carbonyl_mask,
        LOGICAL_NOT=np.logical_not,
        ISNAN=np.isnan,
        ):
    """."""

    # no need to compute isnan in whole coordinates because coordinates
    # either are all nan or are all numbers, hence we use index 0
    coords_in_usagee = LOGICAL_NOT(ISNAN(coords[:, 0]))
    assert coords_in_usagee.shape == (coords.shape[0], )

    elements = atom_labels[coords_in_usagee].astype('<U1')
    # numeric elements refer to the different H atoms as described
    # in core/build_definitions.py
    elements[np.char.isnumeric(elements)] = 'H'

    rows, *_ = vdw_clash_by_threshold(
        coords[coords_in_usagee, :],
        atom_labels[coords_in_usagee],
        elements,
        False,
        False,
        residue_numbers[coords_in_usagee],
        residues_apart=2,
        )

    is_valid = all((
        not np.any(rows),
        ))

    if np.any(rows): # not valid
        return 1  # simulates positive energy
    else:
        return -1  # simulates negative energy
