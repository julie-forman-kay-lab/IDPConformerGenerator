"""
Test raw definitions.

No enphasis is placed on actual definitions.
"""
import pytest

from idpconfgen.core.definitions import aa1to3, aa3to1


# TODO: complete this tests
@pytest.mark.parametrize(
    'aadict',
    [aa3to1, aa1to3],
    )
def test_aa_letter_dictionaries(aadict):
    """Test amino acids letter conversion dictionaries."""
    assert isinstance(aadict, dict)
    # there are 20 natural aminoacids + 3 His protonation states
    assert len(aadict) == 23
