"""
Contain different parsing strategies for different files.
"""
from idpconfgen import Path, log


class DSSPParser:
    """
    Provides an interface for `DSSP files`_.

    .. _DSSP files: https://github.com/cmbi/xssp
    """
    
    def __init__(self, *, fin=None, data=None, pdbid=None):
        
        if fin and Path(fin).exists():
            self.read_dssp_data(Path(fin).read_text())
        elif data:
            self.read_dssp_data(data)
        else:
            self.data = None

        self.pdbid = None

    def read_dssp_data(self, data):
        self.data = data
        

