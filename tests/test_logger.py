"""Test logger."""
import pytest

from idpconfgen import Path, log
from idpconfgen.logger import S, T, init_files, report_on_break
from idpconfgen.core.exceptions import IDPConfGenException


def test_init_files():
    """Test init log files."""
    init_files(log, 'dummy')
    paths = [Path('dummy').with_suffix(p) for p in ['.log', '.error', '.debug']]
    assert all(p.exists() for p in paths)
    for p in paths:
        p.unlink()


def test_T():
    """Test T formatter."""
    logmsg = T('my title {}', 'IDP')
    assert str(logmsg) == '\n* My Title IDP ...'


@pytest.mark.parametrize(
    'msg,args,spacer,indent,expected',
    [
        (
            'a log message with param {}',
            'IDP',
            '+',
            8,
            '++++++++a log message with param IDP',
            ),
        ('string {}', (), '', 1, 'string ()'),
        ('string', (), ' ', 4, '    string'),
        ],
    )
def test_S(msg, args, spacer, indent, expected):
    """Test S formatter."""
    sobj = S(msg, args, spacer=spacer, indent=indent)
    assert str(sobj) == expected


def test_report_on_break():
    """Test record func on error to file."""
    def funca(a, b, c=1, d=2):
        raise TypeError

    ext = 'testing_report_on_break'
    #execute = report_on_break(TypeError, ext=ext)(funca)
    with pytest.raises(TypeError):
        report_on_break(
            funca,
            'idp', 'confgen', c=range(10), d=dict.fromkeys('qwerty'),
            exception=TypeError,
            ROB_ext=ext,
            )

    errfiles = list(Path.cwd().glob(f'*.{ext}'))
    assert len(errfiles) > 0
    for p in errfiles:
        p.unlink()
