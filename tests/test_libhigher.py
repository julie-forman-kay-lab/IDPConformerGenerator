"""
Test lib higher routines.
"""
from collections import defaultdict

import pytest

from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libhigherlevel import *
from . import tcommons


def test_timer_torsion():
    """Test creation of trimer/torsion library."""
    results = defaultdict(dict)
    read_trimer_torsion_planar_angles(tcommons.EXPL_A, results)

    expected = read_dictionary_from_disk(tcommons.expl_a_bgeo)

    assert set(results.keys()) == set(expected.keys())

    for key in results.keys():
        for k2 in results[key].keys():
            r = results[key][k2]
            e = expected[key][k2]

            assert all(r[i] - r[i] < 0.0000000000001 for i in range(len(r)))

    return

