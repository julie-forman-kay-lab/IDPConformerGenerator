"""Test client for PDB Downloader."""
import shutil
from idpconfgen import Path
from idpconfgen.cli_pdbdownloader import main

from . import tcommons


def test_main():
    """Test main from cli_pdbdownloader."""
    folderout = Path('downloaded_test')
    folderout.mkdir(exist_ok=True)
    main(
        [tcommons.cull],
        destination=folderout,
        update=True,
        )
    shutil.rmtree(folderout)
