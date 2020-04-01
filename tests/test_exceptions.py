"""Test custom Exceptions."""
import inspect

import hypothesis.strategies as st
import pytest
from hypothesis import given

#from idpconfgen import contactus as CONTACTUS
from idpconfgen.core import exceptions as EXCPTNS


class MyExceptionNoMsg(EXCPTNS.IDPConfGenException):
    pass


class MyExceptionUnformattable(EXCPTNS.IDPConfGenException):
    errmsg = 'Unformattable error message.'


class MyExceptionFormattable(EXCPTNS.IDPConfGenException):
    # Formattable error message: {}.
    errmsg = 'fer: {}'


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
    """Test init with errmsg=None, should be ignored."""
    err = EXCPTNS.IDPConfGenException(errmsg=errmsg)
    assert str(err) == EXCPTNS.IDPConfGenException.errmsg


### new

@pytest.mark.parametrize(
    'args,expected',
    [
        ([None], 'fer: None'),
        (['some error'], 'fer: some error'),
        (['some error {}', 1], 'some error 1'),
        (['some error {} {}', 1, 2], 'some error 1 2'),
        (['some error {} {} {}', 1, 2, 'asd'], 'some error 1 2 asd'),
        ]
    )
def test_IDPExceptionFormattableError(args, expected):
    """Test IDPConfGenExceptions with Formattable Error."""
    err = MyExceptionFormattable(*args)
    assert str(err) == expected


@pytest.mark.parametrize(
    'args,expected',
    [
        ([], MyExceptionUnformattable.errmsg),
        ([None], 'None'),
        (['some error'], 'some error'),
        (['some error {}', 1], 'some error 1'),
        (['some error {} {}', 1, 2], 'some error 1 2'),
        (['some error {} {} {}', 1, 2, 'asd'], 'some error 1 2 asd'),
        ]
    )
def test_IDPExceptionUnformattableError(args, expected):
    """Test IDPConfGenExceptions with Unformattable Error."""
    err = MyExceptionUnformattable(*args)
    assert str(err) == expected


@pytest.mark.parametrize(
    'args,expected',
    [
        ([], MyExceptionNoMsg.errmsg),
        ([None], 'None'),
        (['some error'], 'some error'),
        (['some error {}', 1], 'some error 1'),
        (['some error {} {}', 1, 2], 'some error 1 2'),
        (['some error {} {} {}', 1, 2, 'asd'], 'some error 1 2 asd'),
        ]
    )
def test_IDPExceptionBase(args, expected):
    """Test IDPConfGenException Base has unformattable behaviour."""
    err = MyExceptionNoMsg(*args)
    assert str(err) == expected


# #old




#@pytest.mark.parametrize(
#    'args',
#    [
#        (['some error']),
#        (['some error {}', 1]),
#        (['some error {} {}', 1, 2]),
#        (['some error {} {} {}', 1, 2, 'asd']),
#        ]
#    )
#def test_IDPCalcException_with_args(args):
#    """Test IDPCalcException to errmsg receive."""
#    err = EXCPTNS.IDPConfGenException(*args)
#    assert err.errmsg == args[0]
#    assert str(err) == args[0].format(*args[1:])
#
#
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
