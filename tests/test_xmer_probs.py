"""Test xmer_probs component."""
from pathlib import Path

import pytest

from idpconfgen.components.xmer_probs import *

from . import tcommons


xmer_probs = Path(tcommons.data_folder, 'xmer_prob.col')


def test_read_xmer_prob():
    """."""
    result = read_xmer_probs_from_file(tcommons.xmer_probs)
    expected = XmerProbs([1, 2, 5], [2, 4, 6])
    assert result.size == expected.size
    assert result.probs == expected.probs


def test_prepare_xmerprob():
    """."""
    result = prepare_xmer_probs(xmer_probs)
    expected = XmerProbs([1, 2, 5], [2, 4, 6])
    assert result.size == expected.size
    assert result.probs == expected.probs


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
        ],
    )
def test_prepare_xmerprob_error(value, error):
    """."""
    with pytest.raises(error):
        result = prepare_xmer_probs(value)


def test_isXmer():
    """."""
    a = XmerProbs([1, 2], [3, 4])
    assert is_XmerProbs(a)
