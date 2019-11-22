"""Test logger."""
import pytest

from idpconfgen import Path, log
from idpconfgen.logger import S, T, init_files


def test_init_files():
    """Test init log files."""
    init_files(log, 'dummy')
   
    files_created = [
        Path(f).exists() for f in ['dummy.log', 'dummy.error', 'dummy.debug']
        ]
    
    assert all(files_created)


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
        ],
    )
def test_S(msg, args, spacer, indent, expected):
    """Test S formatter."""
    sobj = S(msg, args, spacer=spacer, indent=indent)
    assert str(sobj) == expected
