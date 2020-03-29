"""Test custom Exceptions."""
import inspect

import hypothesis.strategies as st
import pytest
from hypothesis import given

#from idpconfgen import contactus as CONTACTUS
from idpconfgen.core import exceptions as EXCPTNS


EXCPT_classes = inspect.getmembers(EXCPTNS, predicate=inspect.isclass)
error_classes = [
    t[1] for t in EXCPT_classes
        if issubclass(t[1], EXCPTNS.IDPConfGenException)
            and t[0].endswith('Error')
    ]


@pytest.fixture(params=error_classes)
def ErrorClass(request):
    """Custom Error Classes in exception module."""
    return request.param


def test_all_errors_names_end_in_error():
    """Test whether all custom error classes end with 'Error'."""
    endswitherror = [t[1] for t in EXCPT_classes if t[0].endswith('Error')]
    subclss = [
        t[1] for t in EXCPT_classes
            if issubclass(t[1], EXCPTNS.IDPConfGenException)
        ]
    assert len(endswitherror) == len(subclss) - 1 # IDPConfGenException itself


def test_IDPConfGenException_type():
    """Test IDPConfGenException is Exception."""
    assert issubclass(EXCPTNS.IDPConfGenException, Exception)


def test_IDPConfGenExc_no_error_mg_0():
    """Test clean instation gives error msg."""
    errmsg = EXCPTNS.IDPConfGenException.errmsg
    assert str(EXCPTNS.IDPConfGenException()) == errmsg


@given(st.none())
def test_IDPConfGenExc_errmsg_None(errmsg):
    """Test init with errmsg=None."""
    err = EXCPTNS.IDPConfGenException(errmsg=errmsg)
    assert str(err) == EXCPTNS.IDPConfGenException.errmsg


@pytest.mark.parametrize(
    'args',
    [
        (['some error']),
        (['some error {}', 1]),
        (['some error {} {}', 1, 2]),
        (['some error {} {} {}', 1, 2, 'asd']),
        ]
    )
def test_IDPCalcException_with_args(args):
    """Test IDPCalcException to errmsg receive."""
    err = EXCPTNS.IDPConfGenException(*args)
    assert err.errmsg == args[0]
    assert str(err) == args[0].format(*args[1:])


@pytest.mark.parametrize(
    'args,errmsg',
    [
        (['this should be ignored {} {}', 1, 2], 'some error.'),
        ]
    )
def test_IDPCalcException_errmsg(args, errmsg):
    """Test IDPCalcException to errmsg without formatting args."""
    err = EXCPTNS.IDPConfGenException(*args, errmsg=errmsg)
    assert err.errmsg == errmsg


def test_ErrorClasses_are_IDPConfGenExc_subclasses(ErrorClass):
    """Is subclass of IDPCalcException."""
    assert issubclass(ErrorClass, EXCPTNS.IDPConfGenException)
