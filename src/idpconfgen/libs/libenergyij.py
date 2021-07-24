"""
Functions to calculate energy from i-j pairs evaluation.

Contains only energy terms that evaluate energy for each atom-atom pair
separately.
"""
import numpy as np
from numba import njit

from idpconfgen import log
from idpconfgen.logger import S


post_calc_options = [
    'pairs',
    'whole',
    ]
"""Option keys to apply after calculating energy value for each pair."""

default_post_calc_option = post_calc_options[0]


def init_lennard_jones_calculator(
        acoeff,
        bcoeff,
        postf=default_post_calc_option,
        ):
    """
    Calculate Lennard-Jones full pontential.

    The LJ potential is calculated fully and no approximations to
    proximity of infinite distance are considered.

    Parameters
    ----------
    acoeff, bcoeff : np.ndarray, shape (N, 3), dtype=np.float
        The LJ coefficients prepared already for the ij-pairs upon which
        the resulting function is expected to operate.
        IMPORTANT: it is up to the user to define the coefficients such
        that resulting energy is np.nan for non-relevant ij-pairs, for
        example, covalently bonded pairs, or pairs 2 bonds apart.

    postf : str
        There are different implementations of the energy calculation.
        Current options are:
        "whole": applies np.nansum after calculating the energy per
        pair.
        "pairs": returns the energies per pair.

    Returns
    -------
    numba.njitted func
        Function closure with registered `acoeff`s and `bcoeff`s that
        expects an np.ndarray of distances with same shape as `acoeff`
        and `bcoeff`: (N,). The `func` return value depends on the
        `postf` options.
    """
    @njit
    def calculate(distances_ij):
        ar = acoeff / (distances_ij ** 12)
        br = bcoeff / (distances_ij ** 6)
        return ar - br

    @njit
    def calculate_nansum(distances_ij, NANSUM=np.nansum):
        energy_ij = calculate(distances_ij)
        return NANSUM(energy_ij)

    # has to be aligned with post_calc_options
    _options = [
        calculate,
        calculate_nansum,
        ]

    options = {k: v for k, v in zip(post_calc_options, _options)}

    log.info(S(f'Lennard-Jones type configured to: {postf!r}'))
    return options[postf]


def init_coulomb_calculator(charges_ij, postf=default_post_calc_option):
    """
    Calculate Coulomb portential.

    Parameters
    ----------
    charges_ij : np.ndarray, shape (N, 3), dtype=np.float
        The `charges_ij` prepared already for the ij-pairs upon which
        the resulting function is expected to operate.
        IMPORTANT: it is up to the user to define the charge such
        that resulting energy is np.nan for non-relevant ij-pairs, for
        example, covalently bonded pairs, or pairs 2 bonds apart.

    postf : str
        There are different implementations of the energy calculation.
        Current options are:
        "whole": applies np.nansum after calculating the energy per
        pair.
        "pairs": returns the energies per pair.

    Returns
    -------
    numba.njitted func
        Function closure with registered `charges_ij` that expects an
        np.ndarray of distances with same shape as `acoeff` and `bcoeff`:
        (N,). The `func` return value depends on the `postf` options.
    """
    @njit
    def calculate(distances_ij):
        return charges_ij / distances_ij

    @njit
    def calculate_nansum(distances_ij, NANSUM=np.nansum):
        return NANSUM(calculate(distances_ij))

    # has to be aligned with post_calc_options
    _options = [
        calculate,
        calculate_nansum,
        ]

    options = {k: v for k, v in zip(post_calc_options, _options)}

    log.info(S(f'Coulomb type configured to: {postf!r}'))
    return options[postf]


def energycalculator_ij(distf, efuncs):
    """
    Calculate the sum of energy terms.

    This function works as a closure.

    Accepts only energy terms that compute for non-redundant ij-pairs.

    Energy terms must have distances ij as unique positional parameter,
    and should return an integer.

    Example
    -------
    >>> ecalc = energycalculator_ij(calc_ij_pair_distances, [...])
    >>> total_energy = ecalc(coords)

    Where `[...]` is a list containing energy term functions.

    See Also
    --------
    init_lennard_jones_calculator
    init_coulomb_calculator

    Parameters
    ----------
    distf : func
        The function that will be used to calculate ij-pair distances
        on each call. If performance is a must, this function should be
        fast. `distf` function should receive `coords` as unique
        argument where `coords` is a np.ndarray of shape (N, 3), where N
        is the number of atoms, and 3 represents the XYZ coordinates.
        This function should return a np.ndarray of shape
        (N * (N - 1)) / 2,), dtype=np.float.

    efuncs : list
        A list containing the energy terms functions. Energy term
        functions are prepared closures that accept the output of
        `distf` function.

    Returns
    -------
    func
        A function that accepts coords in the form of (N, 3). The
        coordinates sent to the resulting function MUST be aligned with
        the labels used to prepare the `efuncs` closures.
    """
    def calculate(coords):
        dist_ij = distf(coords)
        energy = 0
        for func in efuncs:
            energy += func(dist_ij)
        return energy
    return calculate
