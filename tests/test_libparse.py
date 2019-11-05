
from idpconfgen import Path
from idpconfgen.libs import libparse


class TestDSSPParser:

    dssp_file = Path(
        Path(__file__).myparents(),
        'data',
        '1A12_A.dssp',
        )
    
    dsspdata = dssp_file.read_text()

    def test_dssp1(self):
        """Test initiation empty args.""" 
        libparse.DSSPParser()

    def test_dssp2(self):
        """Test initiation from data."""
        libparse.DSSPParser(data=self.dsspdata)

    def test_dssp3(self):
        """Test initiation from file."""
        libparse.DSSPParser(fin=self.dssp_file)

    
