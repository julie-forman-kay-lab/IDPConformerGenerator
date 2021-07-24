"""Test libenergy ij."""
import numpy as np
import pytest

from idpconfgen.libs import libenergyij


def test_lennard_jones_calculator():
    c = libenergyij.init_lennard_jones_calculator(
        np.arange(10),
        np.arange(10),
        'pairs',
        )
    result = c(np.arange(1, 11))
    assert isinstance(result, np.ndarray)
    assert result.ndim == 1
    assert result.size == 10


def test_lennard_jones_calculator_whole():
    c = libenergyij.init_lennard_jones_calculator(
        np.arange(10),
        np.arange(10),
        'whole')
    result = c(np.arange(1, 11))
    assert isinstance(result, float)


def test_energy_calculator_ij():
    """."""
    def dist_calc(coords):
        return coords - 1

    def calc_dummy(arg1, arg2=5):
        def calc(distances):
            return sum(distances) + arg1 - arg2
        return calc

    def calc_dummy2(arg3, arg6=5):
        def calc(distances):
            return sum(distances + arg3 - arg6)
        return calc

    d1 = calc_dummy(2, 3)
    d3 = calc_dummy2(2, arg6=1)
    calc = libenergyij.energycalculator_ij(dist_calc, (d1, d3))
    coords = np.arange(10)
    energy = calc(coords)
    assert energy == 79


