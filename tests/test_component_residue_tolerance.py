"""Test residue tolerance."""
import argparse

import pytest

from idpconfgen.components import residue_tolerance as rt



option_5 = {
    'A': 'A',
    'C': 'C',
    'D': 'D',
    'E': 'E',
    'F': 'FY',
    'G': 'G',
    'H': 'H',
    'I': 'I',
    'K': 'K',
    'M': 'M',
    'N': 'N',
    'P': 'P',
    'Q': 'Q',
    'R': 'R',
    'S': 'S',
    'T': 'T',
    'W': 'W',
    'Y': 'YF',
    }


option_53 = {
    'A': 'A',
    'C': 'CY',
    'D': 'D',
    'E': 'E',
    'F': 'FY',
    'G': 'G',
    'H': 'H',
    'I': 'IVM',
    'K': 'K',
    'M': 'MI',
    'N': 'N',
    'P': 'P',
    'Q': 'Q',
    'R': 'R',
    'S': 'S',
    'T': 'T',
    'W': 'W',
    'Y': 'YF',
    }


def test_argument_add_usubs():
    """Test adding argument."""
    ap = argparse.ArgumentParser()
    rt.add_substitution_groups(ap)
    cmd = ap.parse_args(['-usubs', "\"{'A': 'AC'}\"" ])
    result = vars(cmd)
    assert "residue_tolerance" in result
    assert result["residue_tolerance"] == {'A': 'AC'}


@pytest.mark.skip
def test_argument_add_edss50():
    """Test adding argument."""
    ap = argparse.ArgumentParser()
    rt.add_substitution_groups(ap)
    cmd = ap.parse_args(['-edss50'])
    assert "residue_tolerance" in vars(cmd)


@pytest.mark.skip
def test_argument_add_edss50_5():
    """Test adding argument 5."""
    ap = argparse.ArgumentParser()
    rt.add_substitution_groups(ap)
    cmd = ap.parse_args(['-edss50', '5'])
    assert vars(cmd.residue_tolerance) == option_5


@pytest.mark.skip
def test_argument_add_edss50_5_3():
    """Test adding argument 0."""
    ap = argparse.ArgumentParser()
    rt.add_substitution_groups(ap)
    cmd = ap.parse_args(['-edss50', '5', '3'])
    assert vars(cmd.residue_tolerance) == option_53
