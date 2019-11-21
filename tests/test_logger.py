"""Test logger."""
import pytest

from idpconfgen import log, Path
from idpconfgen.logger import init_files, S, T


def test_init_files():
    """Test init log files."""
    init_files(log, 'dummy')
   
    files_created = [
        Path(f).exists() for f in ['dummy.log', 'dummy.error', 'dummy.debug']
        ]
    
    assert all(files_created)


def test_T():
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
def test_S(msg,args,spacer,indent,expected):
    sobj = S(msg, args, spacer=spacer, indent=indent)
    assert str(sobj) == expected
    
