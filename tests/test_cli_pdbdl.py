"""Test client for PDB Downloader."""
import os

import pytest

from idpconfgen import Path
from idpconfgen.cli_pdbdownloader import main

from . import tcommons


@pytest.mark.parametrize(
    'ncores',
    [1, os.cpu_count() - 1],
    )
def test_main(ncores):
    """Test main from cli_pdbdownloader."""
    with tcommons.TmpFolder('downloaded_test') as DF:
        main(
            [tcommons.cull],
            destination=DF,
            update=True,
            ncores=ncores,
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
            lines = [li for li in Path(DF, p).read_text().split('\n') if li]
            assert len(lines) == l
