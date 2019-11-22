"""Test libcalc."""
import numpy as np
import pytest

from idpconfgen.libs import libcalc


backbone_seed = [
    {
        'N': np.array([0.000, 0.000, 0.000]),
        'CA': np.array([1.458, 0.000, 0.000]),
        'C': np.array([2.009, 1.420, 0.000]),
        },
    ]


@pytest.mark.parametrize(
    'av,bv,cv,expected',
    [
        (
            np.array([0.000, 0.000, 0.000]),
            np.array([1.458, 0.000, 0.000]),
            np.array([2.009, 1.420, 0.000]),
            (
                np.array([ 0.,  0., -1.]),
                np.array([-1.,  0.,  0.]),
                np.array([-0., -1., -0.]),
                ),
            ),
        (
            np.array([1.000, 1.000, 1.000]),
            np.array([1.000, 1.000, 1.000]),
            np.array([1.000, 1.000, 1.000]),
            (
                np.array([ 0.,  0., 0.]),
                np.array([ 0.,  0., 0.]),
                np.array([ 0.,  0., 0.]),
                ),
            ),
        ],
    )
def test_make_axis_vectors_1(av, bv, cv, expected):
    """Test axis vectors from coordinates."""
    result = libcalc.make_axis_vectors(av, bv, cv)
    assert np.all(np.equal(result, expected))


@pytest.mark.parametrize(
    'av,bv,cv,expected',
    [
        (
            np.array([0.000, 0.000, 0.000]),
            np.array([1.458, 0.000, 0.000]),
            np.array([2.009, 1.420, 0.000]),
            (
                np.array([
                    [-1.,  0., 0.],
                    [0.,  1.,  0.],
                    [0., 0., -1.],
                    ])
                ),
            ),
        ],
    )
def test_RT_to_plane(av, bv, cv, expected):
    """Test RT to plane."""
    result = libcalc.RT_to_plane(av, bv, cv)
    assert np.all(np.equal(result, expected))

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
    print(result)
    assert np.all(np.subtract(result, expected) < 0.000001)


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
        ],
    )
def test_make_coord(theta, phi, rad, parent, xaxis, yaxis, expected):
    """Test make coord from angles, distance of bound and parent atoms."""
    res = libcalc.make_coord(theta, phi, rad, parent, xaxis, yaxis)
    assert np.all(np.equal(res, expected))
