"""
Contain different parsing strategies for different files.
"""
import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.core import definitions as DEFS


class DSSPParser:
    """
    Provides an interface for `DSSP files`_.

    .. _DSSP files: https://github.com/cmbi/xssp
    """
    
    def __init__(self, *, fin=None, data=None, pdbid=None):
        
        self.pdbid = pdbid
        
        if fin and Path(fin).exists():
            self.read_dssp_data(Path(fin).read_text())
        elif data:
            self.read_dssp_data(data)
        else:
            self.data = None

    def read_dssp_data(self, data):
        
        try:
            data = data.split('\n')
        except AttributeError:  # data already a list
            pass
        
        data = [i for i in data if i]  # removes empty strings
       
        try:
            data_header_index = self._finds_data_index(data)
        except IndexError as err:
            raise EXCPTS.DSSPParserError(self.pdbid) from err
            
        self.data = data[data_header_index  + 1:]  # data starts afterthe header

        #self.data = self._finds_data_index(data)
        self.read_sec_structure()

    def read_sec_structure(self):
        """
        Read secondary structure information from DSSP files.
        
        Assigns :attr:`ss`.
        """
        self.ss = list_index_to_array(self.data, index=16)
        self._confirm_ss_data(self.ss)

    @staticmethod
    def _confirm_ss_data(data):
        # confirms data makes sense
        if not all((i in vars(DEFS.dssp_ss_keys).values() for i in data)):
            raise EXCPTS.DSSPSecStructError(data)

    @staticmethod
    def _finds_data_index(data):
        """
        Find index for where item startswith '#'.

        Evaluates after left striping white spaces.
        """
        # brute force index finding
        i = 0
        line = data[i]
        # if not exausted, while breaks with IndexError
        while not line.lstrip().startswith('#'):
            i += 1
            line = data[i]
        else:
            return i


def list_index_to_array(
        list_,
        index=None,
        sObj=None,
        start=None,
        stop=None,
        step=None,
        ):
    """
    Extract slices of strings in lists to an array.
    """
     
    if any((i is not None for i in (start, stop, step))):
        index = slice(start, stop, step)
    
    elif sObj is not None:
        if isinstance(sObj, slice):
            index = sObj
        else:
            raise ValueError(
                'sObj must be slice type: {} given.'.format(type(sObj))
                )

    elif index is not None:
        index = index

    else:
        raise ValueError('At least one of the index args most be provided.')

    len_of_slice = len(list_[0][index])

    array = np.empty(
        len(list_),
        dtype='<U{}'.format(len_of_slice),
        )

    for i, line in enumerate(list_):
        array[i] = line[index]
    
    return array
