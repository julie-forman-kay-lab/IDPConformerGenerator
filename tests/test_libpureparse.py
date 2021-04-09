"""Test libpure.libparse module"""

import pytest

from idpconfgen.libpure.libparse import *


@pytest.mark.parametrize(
    'in1,expected',
    [
        ['1-59', ((1, 59),)],
        ['1-59,80-102,150-170', ((1, 59), (80, 102), (150, 170))],
        ],
    )
def test_parse_number_ranges(in1, expected):
    """Test parse number ranges."""
    assert parse_number_ranges(in1) == expected

@pytest.mark.parametrize(
    'in1',
    [
        '1',
        '1-59,80,100-120',
        ]
    )
def test_parse_number_ranges_ValueError(in1):
    with pytest.raises(ValueError):
        parse_number_ranges(in1)

