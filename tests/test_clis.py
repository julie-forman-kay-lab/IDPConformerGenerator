"""Test clients."""
import shutil

from idpconfgen import Path
from idpconfgen.cli_pdbdownloader import main as pdbdl
from idpconfgen.cli_segext import main as segext
from idpconfgen.cli_ssext import main as ssext
from idpconfgen.cli_fastaext import main as fastaext

from . import tcommons


def test_pdbdl_main_1():
    """Test main from cli_pdbdownloader."""
    folderout = Path('downloaded_test')
    folderout.mkdir(exist_ok=True)
    pdbdl(
        [tcommons.cull],
        destination=folderout,
        update=True,
        )
    shutil.rmtree(folderout)


def test_pdbdl_main_2():
    """Test main from cli_pdbdownloader."""
    folderout = Path('downloaded_test')
    folderout.mkdir(exist_ok=True)
    pdbdl(
        [tcommons.cull],
        destination=folderout,
        update=True,
        ncores=2,
        )
    shutil.rmtree(folderout)


def test_pdbdl_main_3():
    """Test main from cli_pdbdownloader."""
    folderout = Path('downloaded_test')
    folderout.mkdir(exist_ok=True)
    pdbdl([tcommons.cull])
    folderout.rmdir()


def test_fasta_main_1():
    """Test main from cli_fastaext."""
    fn = 'result.fasta'
    fastaext(
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
