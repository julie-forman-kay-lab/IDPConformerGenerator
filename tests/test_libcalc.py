"""Test libcalc."""
from math import radians

import numpy as np
import pytest

from idpconfgen.libs import libcalc


@pytest.mark.parametrize(
    'av,bv,cv,expected',
    [
        (
            np.array([0.000, 0.000, 0.000]),
            np.array([1.458, 0.000, 0.000]),
            np.array([2.009, 1.420, 0.000]),
            (
                np.array([0., 0., -1.]),
                np.array([-1., 0., 0.]),
                np.array([-0., -1., -0.]),
                ),
            ),
        (
            np.array([1.000, 1.000, 1.000]),
            np.array([1.000, 1.000, 1.000]),
            np.array([1.000, 1.000, 1.000]),
            (
                np.array([0., 0., 0.]),
                np.array([0., 0., 0.]),
                np.array([0., 0., 0.]),
                ),
            ),
        ],
    )
def test_make_axis_vectors_1(av, bv, cv, expected):
    """Test axis vectors from coordinates."""
    result = libcalc.make_axis_vectors(av, bv, cv)
    assert np.all(np.abs(np.subtract(result, expected)) < 0.0009)


@pytest.mark.parametrize(
    'av,bv,cv,expected',
    [
        (
            np.array([0.000, 0.000, 0.000]),
            np.array([1.458, 0.000, 0.000]),
            np.array([2.009, 1.420, 0.000]),
            (
                np.array([
                    [-1., 0., 0.],
                    [0., 1., 0.],
                    [0., 0., -1.],
                    ])
                ),
            ),
        ],
    )
def test_rotation_to_plane(av, bv, cv, expected):
    """Test rotation to plane."""
    result = libcalc.rotation_to_plane(av, bv, cv)
    assert np.all(np.abs(np.subtract(result, expected)) < 0.0009)


@pytest.mark.parametrize(
    'theta,phi,distance,expected',
    [
        (
            1,
            1,
            1,
            np.array([0.54030231, 0.45464871, 0.70807342])),
        ],
    )
def test_make_coord_from_angles(theta, phi, distance, expected):
    """Test making a coordinate in space from angles and radius."""
    result = libcalc.make_coord_from_angles(theta, phi, distance)
    assert np.all(np.abs(np.subtract(result, expected)) < 0.0009)


@pytest.mark.parametrize(
    'theta,phi,rad,parent,xaxis,yaxis,expected',
    [
        (
            1,
            1,
            1,
            np.array([1.000, 1.000, 1.000]),
            np.array([1.000, 1.000, 1.000]),
            np.array([1.000, 1.000, 1.000]),
            np.array([1., 1., 1.]),
            ),
        (
            1.1,
            0.3,
            4.5,
            np.array([.4, 1.000, .985]),
            np.array([1.001, 2.000, 2.54]),
            np.array([1.000, .234, .8]),
            np.array([-1.129, -2.281, -1.689]),
            ),
        (
            1,
            1,
            1,
            np.array([.002, 0.005, 0.070]),
            np.array([1123., 2432., 15232.12321]),
            np.array([6543.654, 1543., .0023428]),
            np.array([0.56525141, -0.66535531, -0.41308967]),
            ),
        ],
    )
def test_make_coord(theta, phi, rad, parent, xaxis, yaxis, expected):
    """Test make coord from angles, distance of bound and parent atoms."""
    res = libcalc.make_coord(theta, phi, rad, parent, xaxis, yaxis)
    res = np.round(res, decimals=3)
    assert np.all(np.abs(res - expected) < 0.0009)


@pytest.mark.parametrize(
    'coords, expected',
    [
        (np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]]), np.array([90], dtype=float)),  # noqa: E501
        (np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]), np.array([90, -90], dtype=float)),  # noqa: E501
        (
            np.array((
                (1, 0, 0.000),
                (0, 0, 0.000),
                (0, 1, 0.000),
                (-1, 1, 0.000),
                (-1, 1, 1),
                (0, 1, 1),
                (0, 0, 0)
                )),
            np.array([-180, 90, -0, -45], dtype=float),
            ),
        (
            np.array((
                (1, 0, 0.000),
                (0, 0, 0.000),
                (0, 1, 0.000),
                (-1, 1, 0.000),
                (-1, 1, 1),
                )),
            np.array([-180, 90], dtype=float),
            ),
        ]
    )
def test_calculate_torsions(coords, expected):
    """Tests torsion calculation."""
    result = np.degrees(libcalc.calc_torsion_angles(coords))
    assert np.all(np.equal(result, expected))


@pytest.mark.parametrize(
    'i,e',
    (
        [0, 0],
        [1, 0],
        [2, 0],
        [3, 1],
        [4, 1],
        [5, 1],
        [6, 2],
        [7, 2],
        [8, 2],
        [10, 3],
        )
    )
def test_residue_number_ste3(i, e):
    """
    Test obtain residue number by index.

    Considers negative indexes.

    i = input
    e = expected
    o = output
    """
    o = libcalc.calc_residue_num_from_index(i)
    assert e == o


@pytest.mark.parametrize(
    'i,e',
    (
        [0, 0],
        [1, 0],
        [2, 0],
        [3, 0],
        [4, 1],
        [5, 1],
        [6, 1],
        [7, 1],
        [8, 2],
        )
    )
def test_residue_number_step4(i, e):
    """
    Test obtain residue number by index.

    Considers negative indexes.

    i = input
    e = expected
    o = output
    """
    o = libcalc.calc_residue_num_from_index(i, step=4)
    assert e == o


@pytest.fixture(
    params=[
        [23, 20],
        [-34, -30],
        [180, 180],
        [-179, -180],
        [15, 20],
        [25, 20],
        [35, 40],
        [0, 0],
        [-1, 0],
        [8, 10],
        [5, 0],
        [-5, 0],
        ],
    )
def radians_to_degree_bins_data(request):
    return request.param


def test_rrd10(radians_to_degree_bins_data):
    """
    Test rounding radians to nearest degree.
    """
    x0, expected = radians_to_degree_bins_data
    result = libcalc.round_radian_to_degree_bin_10(radians(x0))
    assert result == expected


def test_rrd10_njit(radians_to_degree_bins_data):
    """
    Test rounding radians to nearest degree.
    """
    x0, expected = radians_to_degree_bins_data
    result = libcalc.rrd10_njit(radians(x0))
    assert result == expected


def test_lennard_jones_calculator():
    c = libcalc.init_lennard_jones_calculator(np.arange(10), np.arange(10))
    c(np.arange(1, 11))


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
    calc = libcalc.energycalculator_ij(dist_calc, (d1, d3))
    coords = np.arange(10)
    energy = calc(coords)
    assert energy == 79
