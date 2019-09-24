import logging
import numpy as np
from pathlib import Path as _Path
import string
import sys

class Path(type(_Path())):
    def str(self):
        return os.fspath(self)


LOGFILE = 'idpcalc_pdb_splitter.log'
ERRORLOG = 'idbpcalc_pdb_splitter.error'

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

FORMATTER = logging.Formatter('%(message)s')

ch.setFormatter(FORMATTER)

log.addHandler(ch)


def _initiate_logfiles():
    
    # initiates log file only if main is run
    logfile = logging.FileHandler(LOGFILE, mode='w')
    logfile.setLevel(logging.DEBUG)
    log.addHandler(logfile)
    logfile.setFormatter(FORMATTER)
    
    errorlog = logging.FileHandler(ERRORLOG, mode='w')
    errorlog.setLevel(logging.ERROR)
    log.addHandler(errorlog)
    errorlog.setFormatter(FORMATTER)

    return

class PDBParams:
    """
    PDB Format string slicing according to:
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
    """
    line_name = slice(0, 6)
    serial = slice(6, 11)
    atom_name = slice(12, 16)
    altloc = slice(16, 17)
    resname = slice(17, 20)
    chainid = slice(21, 22)
    resseq = slice(22, 26)
    icode = slice(26, 27)
    xcoord = slice(30, 38)
    ycoord = slice(38, 46)
    zcoord = slice(46, 54)
    occupancy = slice(54, 60)
    tempfactor = slice(60, 66)
    element = slice(76, 78)
    charge = slice(78, 80)


class PDBStructure:
    
    def __init__(self, data):
        self.rawdata = data
   
    @property
    def n_residues(self):
        try:
            return len(set(self.pdb_array_data[:, 6]))
        except AttributeError:
            log.error('Execute build() before accessing data.')
    
    @property
    def res_list(self):
        res_list = []
        prev_res = None
        for line in self.pdb_array_data:
            res = line[6]
            if res != prev_res:
                res_list.append(res)
            prev_res = res
        return res_list

    def build(self):
        self.data = list(filter(None, self.rawdata.split('\n')))
        self.read_pdb_data_to_array()
    
    def read_pdb_data_to_array(self):
        """
        Transforms PDB data into an array.
        """
        slicers = list(filter(
            lambda x: x[0].startswith(tuple(string.ascii_lowercase)),
            PDBParams.__dict__.items(),
            ))
        
        coordinate_data = list(filter(
            lambda x: x.startswith(('ATOM', 'HETATM', 'ANISOU')),
            self.data,
            ))
        
        self.pdb_array_data = np.empty(
            (len(coordinate_data), len(slicers)),
            dtype='<U8')
        
        for ii, line in enumerate(coordinate_data):
            for column, slicer_item in enumerate(slicers):
                self.pdb_array_data[ii, column] = line[slicer_item[1]]
    

def read_all_pdbs_in_folder_tree(folder):
    
    pdb_list = []
    try:
        for root, subdirs, files in os.walk(folder):
            
            only_pdbs = filter(
                lambda x: x.endswith('.pdb'),
                files,
                )
            
            for file_ in only_pdbs:
                pdb_list.append(Path(root, file_))
    except TypeError:
        pass
    return pdb_list


def load_args():
    
    ap = argparse.ArgumentParser(
        usage=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        ) 
    
    ap.add_argument(
        '-s',
        '--sourcedir',
        help='The folder to recursively search for PDB files.',
        type=Path,
        default=None,
        )
  
    ap.add_argument(
        '-f',
        '--files',
        help="Specific PDB files upon which to perform DSSP.",
        nargs='+',
        )
    
    ap.add_argument(
        '-o',
        '--output_folder',
        help="The folder where splitted PDBs will be saved",
        default='splittedPDBs',
        type=str,
        )
    
    cmd = ap.parse_args()

    return cmd


def main(
        sourcedir=None,
        files=None,
        output_folder=None,
        ):
    
    _initiate_logfiles()
   
    pdb_file_list = read_all_pdbs_in_folder_tree(sourcedir)
    
    try:
        pdb_file_list.extend(files)
    except (TypeError):
        pass
    
    if not pdb_file_list:
        sys.exit('No input provided, nothing to do')
   


    return


class TestPDBStructure:
    
    with open('12AS_A.pdb') as fh:
        pd = PDBStructure(fh.read())
    
    pd.build()

    def test1(self):
        assert self.pd.data[-1]
    def test2(self):    
        assert all(i for i in self.pd.data)
    def test3(self):
        assert len(self.pd.data) == 2559
    def test4(self):
        # columns are created based on PDBParams slicing
        assert self.pd.pdb_array_data.shape == (2559, 15)
    def test5(self):
        # starts at 4
        assert self.pd.n_residues == 327
    def test6(self):
        chunk_size = 30                                                         
        chunks = []    
        previ = 0
        i = 0
        #while i < len(self.pd):
        #    i += chunk_size
#            chunks.append(slice(previ, i))
 #           previ = i
    def test7(self):
        assert len(self.pd.res_list) == self.pd.n_residues
        
        

if __name__ == '__main__':

    cmd = load_args()
    main(**vars(cmd))
