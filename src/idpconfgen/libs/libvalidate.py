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

from idpconfgen.core.definitions import vdW_radii_dict
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resName,
    col_resSeq,
    )


def vdW_clash_common_preparation(
        vdW_elements,
        elements_to_consider,
        residue_numbers,
        residues_apart,
        vdW_radii='tsai1999',
        ):
    """
    Prepare masks for vdW clash calculation.

    .. see-also::
        vdW_clash_calc
    """
    atoms_to_consider = np.isin(vdW_elements, elements_to_consider)
    atc_mask = np.logical_and.outer(atoms_to_consider, atoms_to_consider)
    vdW_vector = np.zeros(vdW_elements.size)

    for atom, radius in vdW_radii_dict[vdW_radii].items():
        vdW_vector[vdW_elements == atom] = radius

    pure_radii_sum = np.add.outer(vdW_vector, vdW_vector)

    residue_distances = np.subtract.outer(residue_numbers, residue_numbers)

    return atc_mask, pure_radii_sum, residue_distances >= residues_apart


def vdW_clash_calc(
        coords,
        atc_mask,
        pure_radii_sum,
        distances_apart,
        ):
    #assert coords.size
    #assert atc_mask.size
    #assert pure_radii_sum.size
    #assert residue_distances.size
    #assert isinstance(residues_apart, int)

    distances = distance.cdist(coords, coords, 'euclidean')

    clashes_raw = distances <= pure_radii_sum

    rows_cols = np.logical_and(
        np.logical_and(clashes_raw, atc_mask),
        distances_apart,
        ).nonzero()

    rows, cols = rows_cols

    return rows, cols, distances[rows_cols], pure_radii_sum[rows_cols]


def vdW_clash(
        coords,
        vdW_elements,
        elements_to_consider,
        residue_numbers,
        residues_apart=3,
        vdW_radii='tsai1999',
        ):
    """
    Calculates vdW clashes from XYZ coordinates and identity masks.

    Parameters
    ----------
    coordinates : numpy array, dtype=float, shape (N, 3)
        The atom XYZ coordinates.

    #distances : numpy array, dtype=float, shape (N, N)
        #All-to-all distance array. The diagonal of this array is 0.

    vdW_elements : numpy array, dtype='<U1', shape (N,)
        The corresponding element character of the different atoms.

    elements_to_consider : list-like
        The elements in `vdW_elements` to consider in the vdW clash
        analysis.

    residue_number : numpy array, dtype=int, shape (N,)
        The residue number corresponding to each atom.

    residues_apart : int
        The minimum number of residues apart to consider for a clash.

    Returns
    -------
    What :func:`vdW_clash_calc` returns.
    """
    atc_mask, pure_radii_sum, distances_apart = \
        vdW_clash_common_preparation(
            vdW_elements,
            elements_to_consider,
            residue_numbers,
            residues_apart,
            vdW_radii=vdW_radii,
            )

    return vdW_clash_calc(
        coords,
        atc_mask,
        pure_radii_sum,
        distances_apart,
        )


def validate_conformer_from_disk(
        name,
        pdb_data,
        elements_to_consider,
        residues_apart,
        vdW_radii,
        ):
    """."""
    s = Structure(pdb_data)
    s.build()
    da = s.data_array
    atom_names = da[:, col_name]
    atom_elements = atom_names.astype('<U1')
    res_numbers = da[:, col_resSeq].astype(np.int)
    coords = s.coords

    rows, cols, distances, threshold = vdW_clash(
        coords,
        atom_elements,
        elements_to_consider,
        res_numbers,
        residues_apart,
        vdW_radii=vdW_radii,
        )

    report = clash_report(da, rows, cols, distances, threshold)

    # rows.size is the number of clashes
    return name, rows.size, report


def clash_report(data_array, pair1, pair2, distances, thresholds):
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
    items = [col_resSeq, col_resName, col_name]
    report = []
    for i, (p1, p2) in enumerate(zip(pair1, pair2)):
        n1, r1, a1 = data_array[p1, items]
        n2, r2, a2 = data_array[p2, items]
        report.append(
            f"{n1:>5} {r1:>3} {a1:>5} - "
            f"{n2:>5} {r2:>3} {a2:>5} - "
            f't: {thresholds[i]:.3f} - '
            f'd: {distances[i]:.3f}'
            )

    return '\n'.join(report)
