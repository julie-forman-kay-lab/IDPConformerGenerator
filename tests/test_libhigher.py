"""
Test lib higher routines.
"""
from collections import defaultdict

import pytest

from idpconfgen.libs.libio import read_dictionary_from_disk, save_dictionary
from idpconfgen.libs.libhigherlevel import *
from . import tcommons


def test_bgeo_database():
    """Test creation of trimer/torsion library."""
    results = defaultdict(dict)
    read_trimer_torsion_planar_angles(tcommons.EXPL_A, results)

    expected = read_dictionary_from_disk(tcommons.expl_a_bgeo)

    assert set(results.keys()) == set(expected.keys())

    for key in results.keys():
        for k2 in results[key].keys():
            r = results[key][k2]
            e = expected[key][k2]

            assert all(r[i] - e[i] < 0.00000000001 for i in range(len(r)))

    return


def test_bgeo_database_convert():
    """Test creation of trimer/torsion library."""
    inputdb = read_dictionary_from_disk(tcommons.expl_a_bgeo)
    results = convert_bond_geo_lib(inputdb)
    expected = read_dictionary_from_disk(tcommons.expl_a_bgeo_converted)

    assert set(results.keys()) == set(expected.keys())

    for btype in results.keys():
        for res in results[btype].keys():
            for pairs in results[btype][res].keys():
                for tor in results[btype][res][pairs].keys():
                    r = results[btype][res][pairs][tor]
                    e = expected[btype][res][pairs][tor]

                    assert all(
                        i - j < 0.00000000001
                        for i, j in zip(e, r)
                        )

    return


def test_bgeo_database_reduce():
    """Test reduce bgeo database to trimer and res."""
    inputdb = read_dictionary_from_disk(tcommons.expl_a_bgeo_converted)
    rtrimer, rres = bgeo_reduce(inputdb)
    etrimer = read_dictionary_from_disk(tcommons.expl_a_bgeo_trimer)
    eres = read_dictionary_from_disk(tcommons.expl_a_bgeo_res)

    assert list(eres.keys()) == list(rres.keys())
    for btype in eres.keys():
        assert list(eres[btype].keys()) == list(rres[btype].keys())
        for res in eres[btype].keys():
            r = rres[btype][res]
            e = eres[btype][res]

            assert all(
                i - j < 0.00000000001
                for i, j in zip(e, r)
                )


    assert list(etrimer.keys()) == list(rtrimer.keys())
    for btype in etrimer.keys():
        assert list(etrimer[btype].keys()) == list(rtrimer[btype].keys())

        for res in etrimer[btype].keys():
            assert list(etrimer[btype][res].keys()) == list(rtrimer[btype][res].keys())  # noqa: E501

            for pair in etrimer[btype][res].keys():

                r = rtrimer[btype][res][pair]
                e = etrimer[btype][res][pair]

                assert all(
                    i - j < 0.00000000001
                    for i, j in zip(e, r)
                    )


