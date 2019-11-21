"""Test logger."""

from idpconfgen import log, Path
from idpconfgen.logger import init_files, S, T


def test_init_files():
    """Test init log files."""
    init_files(log, 'dummy')
   
    files_created = [
        Path(f).exists() for f in ['dummy.log', 'dummy.error', 'dummy.debug']
        ]
    
    assert all(files_created)


