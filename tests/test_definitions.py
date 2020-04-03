"""
This module tests interfaces in definitions.

No enphasis is placed on actual definitions.
"""
import pytest

from idpconfgen.core.definitions import (
    aa3to1,
    aa1to3,
    )


@pytest.mark.parametrize(
    'aadict',
    [aa3to1, aa1to3],
    )
def test_aa_letter_dictionaries(aadict):
    assert isinstance(aadict, dict)
    assert len(aadict) == 20  # there are 20 natural aminoacids
