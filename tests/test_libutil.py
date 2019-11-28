"""Test libutil."""
import pytest

from idpconfgen.libs import libutil as UTIL


def test_random_fragment_return():
    """Test return type is slice object."""
    result = UTIL.random_fragment([1, 2])
    assert isinstance(result, slice)


@pytest.mark.parametrize(
    'in1,fragsize,expected',
    [
        (list(range(1000)), 7, 7),
        (list(range(1000)), 0, 0),
        (list(range(1000)), None, 1000),
        ([], None, 0),
        ([], 0, 0),
        ('abcdefgh', 3, 3),
        ((1,2,3,4), 1, 1),
        ],
    )
def test_random_fragment(in1, fragsize, expected):
    """
    Test fragment has expected length.

    Parametrize
    -----------
    1: list with fragsize < len(list)
    2: list with fragsize == 0, should return empty list
    3: list with fragsize is None, should return whole list
    4: empty list with fragsize is None, should return empty list
    5: empty list with fragsize is None, should return empty list
    6: functionality for strings
    7: functionality for tuples
    """
    result = UTIL.random_fragment(in1, fragsize)
    assert len(in1[result]) == expected


@pytest.mark.parametrize(
    'in1,fragsize,error',
    [
        ([], 7, ValueError),
        (list(range(100)), 700, ValueError),
        (9, 2, TypeError),
        (9.0, 2, TypeError),
        ],
    )
def test_random_fragment_errors(in1, fragsize, error):
    """
    Test errors raised with input.

    Parametrize
    -----------
    1: Empty list with fragsize > 0
    2: list with fragsize > len(list)
    3: int
    4: float
    """
    with pytest.raises(error):
        UTIL.random_fragment(in1, fragsize)
