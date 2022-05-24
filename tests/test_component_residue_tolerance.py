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


option_0 = {
    'A': 'ASG',
    'C': 'CSVH',
    'D': 'D',
    'E': 'E',
    'F': 'FVH',
    'G': 'GA',
    'H': 'HFRC',
    'I': 'IT',
    'K': 'KQ',
    'M': 'MT',
    'N': 'NSTD',
    'P': 'P',
    'Q': 'QKR',
    'R': 'RHQW',
    'S': 'STCNA',
    'T': 'TVSMIN',
    'W': 'WR',
    'Y': 'YWH',
    }


@pytest.mark.parametrize(
    'name',
    ['-urestol', '--user-residue-tolerance'],
    )
def test_argument_add_usubs(name):
    """Test adding argument."""
    ap = argparse.ArgumentParser()
    rt.add_res_tolerance_groups(ap)
    cmd = ap.parse_args([name, '{"A": "AC"}'])
    result = vars(cmd)
    assert "residue_tolerance" in result
    assert result["residue_tolerance"] == {'A': 'AC'}


@pytest.mark.parametrize(
    'name',
    ['-edss50', '--edss50-residue-tolerance'],
    )
def test_argument_add_edss50_5(name):
    """Test adding argument."""
    ap = argparse.ArgumentParser()
    rt.add_res_tolerance_groups(ap)
    cmd = ap.parse_args([name, '5'])
    result = vars(cmd)
    assert "residue_tolerance" in result
    assert result["residue_tolerance"] == option_5


@pytest.mark.parametrize(
    'name',
    ['-edss50', '--edss50-residue-tolerance'],
    )
def test_argument_add_edss50_5_3(name):
    """Test adding argument 0."""
    ap = argparse.ArgumentParser()
    rt.add_res_tolerance_groups(ap)
    cmd = ap.parse_args(['-edss50', '5', '3'])
    result = vars(cmd)
    assert result["residue_tolerance"] == option_53


@pytest.mark.parametrize(
    'name',
    ['-edss50', '--edss50-residue-tolerance'],
    )
def test_argument_add_edss50_0(name):
    """Test adding argument 0."""
    ap = argparse.ArgumentParser()
    rt.add_res_tolerance_groups(ap)
    cmd = ap.parse_args(['-edss50', '0'])
    result = vars(cmd)
    assert result["residue_tolerance"] == option_0
