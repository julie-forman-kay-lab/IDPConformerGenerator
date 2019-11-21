"""Test libcheck."""
import pytest

from idpconfgen.libs import libcheck

@libcheck.argstype(int, float)
@libcheck.kwargstype(str)
def dummy_function(arg1, arg2, arg3='string'):
    print(arg1, arg2, arg3)


@pytest.mark.parametrize(
    'arg1,arg2,arg3',
    [
        ('a', 1.0, 's'),
        (1, 'a', 's'),
        (1, 1.0, 1),
        ([1], ['a'], 's'),
        (1, 1.0, ['s']),
        ]
    )
def test_argstype(arg1, arg2, arg3):
    with pytest.raises(TypeError):
        dummy_function(arg1, arg2, arg3=arg3)
