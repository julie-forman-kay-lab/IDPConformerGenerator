"""Test xmer_probs component."""
import argparse
import os
from pathlib import Path
from math import isclose

import pytest

from idpconfgen.components.xmer_probs import *
from idpconfgen.libs.libcalc import make_seq_probabilities

from . import tcommons


xmer_probs = Path(tcommons.data_folder, 'xmer_prob.col')


def test_read_xmer_prob():
    """."""
    result = read_xmer_probs_from_file(xmer_probs)
    expected = make_xmerprobs([1, 2, 5], [2, 4, 6])
    assert result.size == expected.size
    for i, j in zip(result.probs, expected.probs):
        assert isclose(i, j)


def test_prepare_xmerprob():
    """."""
    result = prepare_xmer_probs(xmer_probs)
    expected = make_xmerprobs([1, 2, 5], [2, 4, 6])
    assert result.size == expected.size
    for i, j in zip(result.probs, expected.probs):
        assert isclose(i, j)


def test_prepare_xmerprob2():
    """."""
    result = read_xmer_probs_file_or_default(None)
    assert result is default_XmerProbs


def test_prepare_xmerprob3():
    """."""
    result = prepare_xmer_probs(None)
    assert result is default_XmerProbs


@pytest.mark.parametrize(
    'value,error',
    [
        ('', IsADirectoryError),
        ('filedoesnotexistnevereverbefore', FileNotFoundError),
        (1, (ValueError, TypeError)),
        (1.1, (ValueError, TypeError)),
        ({1: 2}, (ValueError, TypeError)),
        ({}, (ValueError, TypeError)),
        ([], (ValueError, TypeError)),
        ],
    )
def test_prepare_xmerprob_error(value, error):
    """."""
    with pytest.raises(error):
        result = prepare_xmer_probs(value)


def test_isXmer():
    """."""
    a = make_xmerprobs([1, 2], [3, 4])
    assert is_XmerProbs(a)


def test_Actions1():
    """."""
    ap = argparse.ArgumentParser()
    add_xmer_arg(ap)
    cmd = vars(ap.parse_args([]))
    assert cmd['xmer_probs'] is None


def test_Actions2():
    """."""
    ap = argparse.ArgumentParser()
    add_xmer_arg(ap)
    cmd = vars(ap.parse_args(f'-xp {os.fspath(xmer_probs)}'.split()))
    assert is_XmerProbs(cmd['xmer_probs'])
