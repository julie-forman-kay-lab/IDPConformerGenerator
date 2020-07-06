"""Test custom Exceptions."""
import inspect

import pytest
from hypothesis import given
from hypothesis import strategies as st

from idpconfgen.core import count_string_formatters
from idpconfgen.core import exceptions as EXCPTNS
from idpconfgen.core import has_string_formatters

from .tcommons import random_type


EXCPT_classes = inspect.getmembers(EXCPTNS, predicate=inspect.isclass)
error_classes = [
    t[1] for t in EXCPT_classes
    if issubclass(t[1], EXCPTNS.IDPConfGenException)
    and t[0].endswith('Error')
    ]


@pytest.fixture(params=error_classes)
def ErrorClass(request):
    """Return custom Error Classes in exception module."""
    return request.param


exceptions_with_formattable_erromsg = list(filter(
    lambda x: has_string_formatters(x.errmsg),
    error_classes,
    ))


@pytest.fixture(params=exceptions_with_formattable_erromsg)
def ExcptsFormattable(request):
    """Return IDPConfGen Exceptions with formattable errmsg."""
    return request.param


def test_all_errors_names_end_in_error():
    """Test whether all custom error classes end with 'Error'."""
    endswitherror = [t[1] for t in EXCPT_classes if t[0].endswith('Error')]
    subclss = [
        t[1] for t in EXCPT_classes
        if issubclass(t[1], EXCPTNS.IDPConfGenException)
        ]
    assert len(endswitherror) == len(subclss) - 1  # IDPConfGenException itself


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


def test_IDPExceptionFormattableError(ExcptsFormattable):
    """
    Test Exceptions with formattable errmsgs.

    # have to use random_type() because of the incompatibility
    # betwee pytest.fixtures and hypotesis.given

    HypothesisDeprecationWarning: [...] uses the 'ExcptsFormattable' fixture,
    but function-scoped fixtures should not be used with @given(...) tests,
    because fixtures are not reset between generated examples!
    since="2020-02-29",

    -- Docs: https://docs.pytest.org/en/latest/warnings.html
    """
    num = count_string_formatters(ExcptsFormattable.errmsg)
    args = [random_type() for i in range(num)]
    str(ExcptsFormattable(*args))


@pytest.fixture(
    params=[
        ('some error {}', 1, 'some error 1'),
        ('some error {} {}', 1, 2, 'some error 1 2'),
        ('some error {} {} {}', 1, 2, 'asd', 'some error 1 2 asd'),
        ]
    )
def forcing_messages(request):
    """Formattable messages that override the default errmsg."""
    return request.param


def test_forcing_messages(ErrorClass, forcing_messages):
    """Test behaviour against when formattable string."""
    args, expected = forcing_messages[:-1], forcing_messages[-1]
    err = ErrorClass(*args)
    assert str(err) == expected


@pytest.mark.parametrize(
    'args,errmsg',
    [
        (['this should be ignored {} {}', 1, 2], 'some error.'),
        ]
    )
def test_IDPCalcException_errmsg(ErrorClass, args, errmsg):
    """Test IDPCalcException to errmsg without formatting args."""
    err = ErrorClass(*args, errmsg=errmsg)
    assert err.errmsg == errmsg


def test_ErrorClasses_are_IDPConfGenExc_subclasses(ErrorClass):
    """Is subclass of IDPCalcException."""
    assert issubclass(ErrorClass, EXCPTNS.IDPConfGenException)
