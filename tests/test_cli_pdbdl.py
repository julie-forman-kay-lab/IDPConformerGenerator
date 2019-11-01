
from idpconfgen import Path
from idpconfgen.cli_pdbdownloader import main


file_path = Path(__file__).myparents()


class TestCliPDBdl:

    cull = Path(file_path, 'data', 'cull.list')
    fout = Path(file_path, 'data_out_scratch')

    def test_main(self):
        """Test main from cli_pdbdownloader."""

        main(
            [self.cull],
            destination=self.fout,
            update=True,
            )
