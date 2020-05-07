"""Test cli_fastaext."""
from pathlib import Path

from idpconfgen.cli_fastaext import main

from . import tcommons


def test_main_1():
    """Test main from cli_fastaext."""
    fn = 'result.fasta'
    main(
        [
            tcommons.PFA1,
            tcommons.PFA2,
            tcommons.PFA3,
            ],
        output=fn,
        )
    fout = Path(fn)
    result = fout.read_text()
    fout.unlink()
    assert result == tcommons.pdbs_fasta.read_text().rstrip()
