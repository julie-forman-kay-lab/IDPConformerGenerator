
from idpconfgen import Path
from idpconfgen.cli_pdbdownloader import main
from . import tcommons


class TestCliPDBdl:

    cull = Path(file_path, 'data', 'cull.list')
    fout = Path(file_path, 'data_out_scratch')

    def test_main(self):
        """Test main from cli_pdbdownloader."""

        main(
            [self.cull],
            destination=tcommons.folder_output,
            update=True,
            )
