"""Test libcheck."""
import pytest

from idpconfgen.libs import libcheck

@libcheck.argstype(int, (float,int))
@libcheck.kwargstype(str)
def dummy_function(arg1, arg2, arg3='string'):
    print(arg1, arg2, arg3)


@pytest.mark.parametrize(
    'arg1,arg2,arg3',
    [
        (1.0, 1, 'a'),
        ([1], 1, 'a'),
        ('1', 1, 'a'),
        ((1,), 1, 'a'),
        (1, [1], 'a'),
        (1, '1', 'a'),
        (1, (1,), 'a'),
        (1, 1, 1),
        (1, 1, 1.0),
        (1, 1, ['a']),
        ]
    )
def test_type_error(arg1, arg2, arg3):
    """Test type check with TypeError."""
    with pytest.raises(TypeError):
        dummy_function(arg1, arg2, arg3=arg3)


@pytest.mark.parametrize(
    'arg1,arg2,arg3',
    [
        (1, 1.0, 'a'),
        (1, 1, 'a'),
        ]
    )
def test_types_no_erros(arg1, arg2, arg3):
    """Test args without errors."""
    dummy_function(arg1, arg2, arg3)
