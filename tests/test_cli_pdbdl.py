
from idpconfgen import Path
from idpconfgen.cli_pdbdownloader import main
from . import tcommons


class TestCliPDBdl:

    cull = Path(tcommons.data_folder, 'cull.list')

    def test_main(self):
        """Test main from cli_pdbdownloader."""

        main(
            [self.cull],
            destination=tcommons.folder_output,
            update=True,
            )
