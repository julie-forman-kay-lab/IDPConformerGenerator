"""Test libcalc."""
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
    'labels,expected',
    [
        (['N', 'CA', 'C'], True),
        (np.array(['N', 'CA', 'C'] * 10), False),
        (np.array(['CA', 'N', 'C'] + ['N', 'CA', 'C'] * 9), True),
        (np.array(['N', 'CA', 'C', 'N'] * 10), True),
        (np.array(['N', 'CA', 'O'] * 10), True),
        ]
    )
def test_validate_backbone_labels_for_torsions(labels, expected):
    """Validate Backbone labels for torsions."""
    assert bool(libcalc.validate_backbone_labels_for_torsion(labels)) == expected  # noqa: E501


@pytest.mark.parametrize(
    'coords, expected',
    [
        (np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]]), np.array([90], dtype=float)),  # noqa: E501
        (np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]), np.array([90, -90], dtype=float)),  # noqa: E501
        ]
    )
def test_calculate_torsions(coords, expected):
    """Tests torsion calculation."""
    result = np.degrees(libcalc.calc_torsion_angles(coords))
    assert np.all(np.equal(result, expected))
