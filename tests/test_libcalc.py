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
                array([ 0.,  0., -1.]),
                array([-1.,  0.,  0.]),
                array([-0., -1., -0.]),
                ),
            ),
        ],
    )
def test_make_axis_vectors_1(av, bv, cv, expected):
    result = libcalc.make_axis_vectors(av, bv, cv)
    assert np.all(np.equal(result, expected))
