"""Test client for secondary structure extract."""

from idpconfgen import Path
from idpconfgen.cli_ssext import main
from idpconfgen.libs import libfile

file_path = Path(__file__).myparents()


#class TestCliSSEXT:
#    
#    pdbds = libfile.glob_folder(Path('data', 'pdbs_ssext'), '.pdb')
#
#    def test_main():
#        """Test main cli_ssext."""
#
#        main(self.pdbs)
