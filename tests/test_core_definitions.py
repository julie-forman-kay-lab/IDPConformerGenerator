"""Test core definitions."""

from idpconfgen.core import definitions as DEFS


def test_aa_letters():
    assert len(DEFS.aa3to1) == len(DEFS.aa1to3)
    assert list(DEFS.aa3to1.keys()) == list(DEFS.aa1to3.values())
    assert list(DEFS.aa1to3.keys()) == list(DEFS.aa3to1.values())
