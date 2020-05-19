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

    pdbs = [
        '12E8_H.pdb',
        '16PK_A.pdb',
        '16VP_A.pdb',
        '1A04_A.pdb',
        '6XY7_AAA.pdb',
        ]
    length = [
        1682,
        3135,
        2457,
        1587,
        3704,
        ]
    for p, l in zip(pdbs, length):
        lines = [li for li in Path(folderout, p).read_text().split('\n') if li]
        assert len(lines) == l

    shutil.rmtree(folderout)
