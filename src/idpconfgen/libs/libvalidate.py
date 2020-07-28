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
from scipy.spatial import KDTree
from numpy import array
import math



from idpconfgen.core.definitions import vdW_radii_dict
from idpconfgen.libs.libstructure import Structure



def vdW_clash_common_preparation(
        vdW_elements,
        elements_to_consider,
        residue_numbers,
        residues_apart,
        vdW_radii='tsai1999',
        ):
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

    rows, cols = np.logical_and(
        np.logical_and(clashes_raw, atc_mask),
        distances_apart,
        ).nonzero()

    #print(rows, cols)

    return rows, cols


def vdW_clash(
        coords,
        vdW_elements,
        elements_to_consider,
        residue_numbers,
        residues_apart=3,
        vdW_radii='tsai1999',
        ):
    """
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
    """
    atc_mask, pure_radii_sum, distances_apart = \
        vdW_clash_common_preparation(
            vdW_elements,
            elements_to_consider,
            residue_numbers,
            residues_apart,
            vdW_radii=vdW_radii,
            )

    results = vdW_clash_calc(
        coords,
        atc_mask,
        pure_radii_sum,
        distances_apart,
        )

    return results



def validate_conformer(name, data):
    """
    """
    structure = Structure(data)
    structure.build()
    structure.add_filter_backbone()
    bbcoords = structure.coords
    print(bbcoords)
    return
