"""
************************************
IDPCalculator: PDB to DSSP extractor
************************************

2019
Dr. Julie Forman-Kay Lab
http://abragam.med.utoronto.ca/~JFKlab/

version: 0.3

DESCRIPTION:

    Extracts DSSP information from a PDB file.

USAGE:
    
    Current implementations:
        - extracts secondary structure of multiple files to a single
            ID indexed text file of the format:

            1ABC_D|HHHEEELLLL(...)
            1DEF_H|EEELLLHHHH(...)

            where 1ABC_D is the PDBID_Chain identifier and H, E and L
            are the secondary structure information as extracted from
            the DSSP results and where every result different from H or E
            is reduced to L (loop).

REQUIREMENTS:

    Requires a third party DSSP program to perform the actual
    DSSP calculation. Path to DSSP is given by the '-d' option.

    Parsing the DSSP output is implemented considering the output
    as given by XSSP: https://github.com/cmbi/xssp

HELP:
    
    use the following command for help:
    
    $ python idpcalc_pdb2dssp.py -h

CONTRIBUTORS TO THIS FILE:

- Joao M.C. Teixeira (https://github.com/joaomcteixeira)
"""
# IDPCalculator PDB Downloader is free software:
# you can redistribute it and/or modify
# it under the terms of the LGPL - GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# LGPL - GNU Lesser General Public License for more details.
#
# You should have received a copy of the this license along
# with this library. If not, see <http://www.gnu.org/licenses/>.
import argparse
import contextlib
import logging
import multiprocessing
import os
from pathlib import Path as _Path
import subprocess
import sys
import traceback

version = '0.3'

LOGFILE = 'idpcalc_pdb2dssp.log'
ERRORLOG = 'idpcalc_pdb2dssp.error'

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

FORMATTER = logging.Formatter('%(message)s')

ch.setFormatter(FORMATTER)
log.addHandler(ch)


def _initiate_logfile():
    
    logfile = logging.FileHandler(LOGFILE, mode='w')
    logfile.setLevel(logging.DEBUG)
    logfile.setFormatter(FORMATTER)
    log.addHandler(logfile)
    
    errorlogfile = logging.FileHandler(ERRORLOG, mode='w')
    errorlogfile.setLevel(logging.ERROR)
    errorlogfile.setFormatter(FORMATTER)
    log.addHandler(errorlogfile)


class TitleLog:
    
    def __init__(self, msg):
        self.msg = msg.title()
    
    def __str__(self):
        return '\n* {} ...'.format(self.msg)


class SubLog:
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return '    {}'.format(self.msg)


T = TitleLog
S = SubLog


class DSSPParserError(Exception):
    def __init__(self, pdbid):
        self.pdbid = pdbid
    
    def __str__(self):
        return 'Error while parsing {}'.format(self.pdbid)


class Path(type(_Path())):
    def str(self):
        return os.fspath(self)


class Worker(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, timeout=10):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.timeout = timeout
    
    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                log.info(S(f'exiting... {proc_name}'))
                self.task_queue.task_done()
                break
            result = next_task()
            self.task_queue.task_done()
            self.result_queue.put(result)
        return


class DSSPTask:
    def __init__(self, dssp_exec, pdbpath):
        self.cmd_exec = dssp_exec
        self.pdbpath = pdbpath
        
        self.cmd = [
            self.cmd_exec,
            '-i',
            str(self.pdbpath),
            ]

    def __str__(self):
        return ' '.join(self.cmd)

    def __call__(self):
        log.info(S(f'running {self.cmd[-1]}'))
        result_dssp_string = subprocess.run(
            self.cmd,
            capture_output=True,
            )
        
        decoded = result_dssp_string.stdout.decode('utf8')
        return DSSPParse(Path(self.pdbpath).stem, decoded)


class DSSPParse:
    def __init__(self, pdbid, inputstring):
        self.pdbid = pdbid.upper()
        self.data = inputstring.split('\n')
    
    def get_secondary_structure(self):
        """
        Extracts secondary structure information from DSSP data.
        """
        index_where_data_starts = self._get_index()

        valid_lines = filter(
            lambda line: bool(line),
            self.data[index_where_data_starts:],
            )
        
        self.ssd = [self.parse_ss_char(line) for line in valid_lines]
            
        return
    
    def _get_index(self, token_char='#'):
        
        # index from where structural data starts
        # in DSSP output structural data goes until the end of the file
        data_strip = (l.strip() for l in self.data)
        
        found_index = [
            i for i, line in enumerate(data_strip)
                if line.startswith('#')
            ]
        try:
            return found_index[0] + 1
        except IndexError as e:
            raise DSSPParserError(self.pdbid) from e

    @staticmethod
    def parse_ss_char(line, ss_char_index=16):
        char = line[ss_char_index]
        
        if char not in ('H', 'E'):
            return 'L'
        else:
            return char


class DSSPExec:
    def __init__(
            self,
            command=None,
            pdb_file_list=None,
            ncores=1,
            ):

        self.command = command
        self.pdb_file_list = pdb_file_list
        self.ncores = 1
        return

    @property
    def pdb_file_list(self):
        return self._pdbfl

    @pdb_file_list.setter
    def pdb_file_list(self, pdbfl):
        if all(Path(f).suffix == '.pdb' for f in pdbfl):
            self._pdbfl = pdbfl
        else:
            raise TypeError(f'Not all files have .pdb extension')

    def run(self):
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        workers = [Worker(tasks, results) for i in range(self.ncores)]
        
        for w in workers:
            w.start()
        
        for pdbfile in self.pdb_file_list:
            tasks.put(DSSPTask(self.command, pdbfile))
        
        # Adds a poison pill
        for w in workers:
            tasks.put(None)

        tasks.join()
        
        numjobs = len(self.pdb_file_list)
        self.dssp_results = []
        while numjobs:
            self.dssp_results.append(results.get())
            numjobs -= 1


class DSSPDataBase:

    def __init__(self):
        self._ids = []
        self._ssd = []
    
    def add(self, pdbid_ssd_tuple):
        self._append_ids(pdbid_ssd_tuple[0])
        self._append_ssd(pdbid_ssd_tuple[1])

    def _append_ids(self, pdbid):
        self._ids.append(pdbid)

    def _append_ssd(self, ssd):
        try:
            self._ssd.append(''.join(ssd))
        except AttributeError:
            self._ssd.append(ssd)

    def write(self, destination=None):
        lines = self._prepare_data(self._ids, self._ssd)
        if destination is None:
            sys.stdout.write(''.join(lines))
        else:
            try:
                with open(destination, 'w') as fh:
                    fh.write('\n'.join(lines))
            except TypeError:
                destination.write('\n'.join(lines))
    
    @staticmethod
    def _prepare_data(ids, ssds):
        lines = []
        for pdbid, ssd in zip(ids, ssds):
            lines.append('{}|{}'.format(pdbid, ssd))
        return lines


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
        '-d',
        '--dssp_command',
        help='The path to the DSSP executable file.',
        default=None,
        )
    
    ap.add_argument(
        '-o',
        '--output',
        help=(
            "The output file containing the PDBID and "
            "respective secondary structure information. "
            "Defaults to sys.stdout (prints to console)."
            ),
        type=Path,
        default=None,
        )
    
    ap.add_argument(
        '-n',
        '--ncores',
        help='Number of cores to use',
        default=1,
        type=int,
        )

    cmd = ap.parse_args()
   
    if cmd.dssp_command is None:
        ap.error('DSSP_COMMAND option is mandatory')
    
    if not (cmd.files or cmd.sourcedir):
        ap.error('At least one input source is required')

    return cmd


def main(
        dssp_command=None,
        output=None,
        sourcedir=None,
        files=None,
        ncores=1,
        ):
    _initiate_logfile()
 
    pdb_file_list = read_all_pdbs_in_folder_tree(sourcedir)
    
    try:
        pdb_file_list.extend(files)
    except (TypeError):
        pass
    
    if not pdb_file_list:
        sys.exit('No input provided, nothing to do')
    
    dssp_exec = DSSPExec(
        command=dssp_command,
        pdb_file_list=pdb_file_list,
        ncores=1,
        )

    dssp_exec.run()
    
    ssdb = DSSPDataBase()
    
    for dssp_parser in dssp_exec.dssp_results:
        
        try:
            dssp_parser.get_secondary_structure()
        except Exception as e:
            log.error(traceback.format_exc())
            log.error(repr(e))
            continue

        ssdb.add((dssp_parser.pdbid, dssp_parser.ssd))

    ssdb.write(destination=output)

    return


class TestDSSPParser:

    example_string = \
"""==== Secondary Structure Definition by the program DSSP, CMBI version 3.0.8                          ==== DATE=2019-09-23      .
REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              .
                                                                                                                               .
COMPND                                                                                                                         .
SOURCE                                                                                                                         .
AUTHOR                                                                                                                         .
  327  1  0  0  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .
 14694.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .
  219 67.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .
    9  2.8   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
   50 15.3   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
    1  0.3   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5), SAME NUMBER PER 100 RESIDUES                              .
    2  0.6   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4), SAME NUMBER PER 100 RESIDUES                              .
    1  0.3   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-3), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-2), SAME NUMBER PER 100 RESIDUES                              .
    2  0.6   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-1), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+0), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+1), SAME NUMBER PER 100 RESIDUES                              .
   22  6.7   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+2), SAME NUMBER PER 100 RESIDUES                              .
   26  8.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+3), SAME NUMBER PER 100 RESIDUES                              .
   91 27.8   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4), SAME NUMBER PER 100 RESIDUES                              .
    6  1.8   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5), SAME NUMBER PER 100 RESIDUES                              .
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           .
  0  0  0  0  0  1  2  1  1  0  1  1  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0    RESIDUES PER ALPHA HELIX         .
  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    PARALLEL BRIDGES PER LADDER      .
  2  1  1  0  2  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    ANTIPARALLEL BRIDGES PER LADDER  .
  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    LADDERS PER SHEET                .
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN AUTHCHAIN
    1    4 A A     >        0   0   96      0, 0.0     4,-3.0     0, 0.0     5,-0.3   0.000 360.0 360.0 360.0 157.9   12.5   39.0   28.5                A         A
    2    5 A Y  H  >  +     0   0  158      2,-0.2     4,-2.8     1,-0.2     5,-0.2   0.959 360.0  48.5 -41.6 -56.8   15.6   39.4   26.3                A         A
    3    6 A I  H  > S+     0   0  117      2,-0.2     4,-1.4     1,-0.2    -1,-0.2   0.918 112.7  47.7 -55.5 -51.1   17.8   39.3   29.4                A         A
    4    7 A A  H >> S+     0   0   49      1,-0.2     4,-2.0     2,-0.2     3,-0.7   0.954 114.9  44.5 -54.4 -58.3   16.0   36.2   30.7                A         A
    5    8 A K  H 3X S+     0   0   23     -4,-3.0     4,-2.6     1,-0.2    -1,-0.2   0.832 108.7  58.3 -54.3 -39.3   16.1   34.3   27.4                A         A
    6    9 A Q  H 3X S+     0   0   57     -4,-2.8     4,-1.6    -5,-0.3    -1,-0.2   0.823 109.2  44.9 -61.4 -33.5   19.8   35.3   26.9                A         A
    7   10 A R  H <X S+     0   0  186     -4,-1.4     4,-1.6    -3,-0.7    -2,-0.2   0.794 109.2  55.7 -86.1 -32.3   20.7   33.7   30.1                A         A
    8   11 A Q  H  X S+     0   0   45     -4,-2.0     4,-1.8     2,-0.2     3,-0.3   0.958 108.3  48.7 -52.9 -57.3   18.6   30.6   29.3                A         A
    9   12 A I  H  X S+     0   0   11     -4,-2.6     4,-2.3     1,-0.2     5,-0.3   0.939 110.3  50.8 -53.8 -50.5   20.5   30.1   26.1                A         A
   10   13 A S  H  X S+     0   0   80     -4,-1.6     4,-1.3     1,-0.2    -1,-0.2   0.817 108.6  51.3 -55.7 -41.7   24.0   30.5   27.8                A         A
   11   14 A F  H  X S+     0   0   55     -4,-1.6     4,-2.8    -3,-0.3    -1,-0.2   0.880 111.9  46.9 -64.6 -44.9   23.2   27.9   30.5                A         A
   12   15 A V  H  X S+     0   0    0     -4,-1.8     4,-2.3    -3,-0.2     5,-0.3   0.964 113.2  45.7 -61.1 -54.2   22.0   25.3   28.0                A         A
   13   16 A K  H  X S+     0   0   49     -4,-2.3     4,-0.7     2,-0.2    -1,-0.2   0.787 118.5  45.6 -62.5 -29.9   25.0   25.6   25.6                A         A
   14   17 A S  H >X S+     0   0   80     -4,-1.3     4,-1.9    -5,-0.3     3,-0.6   0.979 111.4  47.9 -80.5 -53.6   27.3   25.6   28.5                A         A
   15   18 A H  H 3X S+     0   0   31     -4,-2.8     4,-1.8     1,-0.3     3,-0.2   0.879 115.0  46.2 -53.4 -50.3   25.9   22.6   30.5                A         A
   16   19 A F  H 3X S+     0   0    0     -4,-2.3     4,-1.7     2,-0.2    -1,-0.3   0.774 104.2  62.6 -67.3 -28.8   25.7   20.4   27.4                A         A
   17   20 A S  H <X S+     0   0   10     -4,-0.7     4,-1.2    -3,-0.6     3,-0.4   0.950 107.9  43.7 -57.8 -42.9   29.3   21.5   26.5                A         A
   18   21 A R  H  X S+     0   0  130     -4,-1.9     4,-1.9     1,-0.2     5,-0.3   0.848 107.1  60.6 -70.4 -32.9   30.3   19.9   29.8                A         A
   19   22 A Q  H  X S+     0   0   16     -4,-1.8     4,-1.6    -5,-0.2    -1,-0.2   0.854 105.3  49.2 -59.1 -34.3   28.0   16.9   28.9                A         A
   20   23 A L  H  X S+     0   0    3     -4,-1.7     4,-1.1    -3,-0.4     6,-0.6   0.824 110.9  48.3 -75.1 -35.7   30.1   16.3   25.9                A         A
   21   24 A E  H  X S+     0   0   77     -4,-1.2     4,-0.7     4,-0.2    -2,-0.2   0.818 118.0  37.8 -74.7 -37.6   33.5   16.4   27.7                A         A
   22   25 A E  H  < S+     0   0  125     -4,-1.9    -2,-0.2     2,-0.2    -3,-0.2   0.918 120.8  42.5 -87.6 -40.5   32.6   14.1   30.5                A         A
   23   26 A R  H  < S+     0   0  115     -4,-1.6    -3,-0.2    -5,-0.3    -2,-0.2   0.772 132.5  21.9 -75.4 -28.9   30.4   11.6   28.6                A         A
   24   27 A L  H  < S-     0   0   20     -4,-1.1    64,-0.2    -5,-0.2    -1,-0.2   0.445  98.2-122.3-123.5   0.3   32.7   11.4   25.6                A         A
   25   28 A G     <  +     0   0   40     -4,-0.7    63,-1.3     1,-0.2     2,-0.4   0.909  55.0 155.4  61.8  44.3   36.1   12.6   26.9                A         A
   26   29 A L  E     -a   88   0A   4     -6,-0.6     2,-0.3    61,-0.2    63,-0.3  -0.880  30.6-151.3-104.2 137.6   36.5   15.4   24.4                A         A
   27   30 A I  E     -a   89   0A 101     61,-2.1    63,-2.3    -2,-0.4     2,-0.3  -0.724  27.1-104.8-102.2 152.9   38.5   18.6   24.9                A         A
   28   31 A E  E     -a   90   0A  91     -2,-0.3     2,-0.3    61,-0.2    63,-0.2  -0.558  38.8-174.5 -76.6 134.1   37.8   21.9   23.4                A         A
   29   32 A V        -     0   0   34     61,-1.8     2,-0.4    -2,-0.3    63,-0.3  -0.883  33.2 -93.1-124.3 155.1   40.0   22.9   20.5                A         A
   30   33 A Q        -     0   0  160     -2,-0.3    61,-0.1    61,-0.1    -2,-0.0  -0.643  43.1-134.4 -73.8 129.1   40.2   26.2   18.5                A         A
   31   34 A A        -     0   0   41     -2,-0.4     2,-0.1    61,-0.2    -1,-0.0  -0.540   3.8-134.8 -87.0 137.6   38.0   26.3   15.4                A         A
   32   35 A P        -     0   0   38      0, 0.0    62,-0.4     0, 0.0     3,-0.1  -0.319  12.5-172.0 -63.3 164.0   38.8   27.5   11.9                A         A
   33   36 A I  S    S+     0   0   81      1,-0.2    34,-2.3    32,-0.1     2,-0.2   0.388  71.7  47.7-128.7 -27.9   36.4   29.7   10.1                A         A
   34   37 A L  E     -H   66   0B  65     32,-0.2     2,-0.3    30,-0.0    32,-0.2  -0.657  62.3-170.1-106.8 165.0   38.1   29.8    6.8                A         A
   35   38 A S  E     -H   65   0B  16     30,-1.9    30,-2.3    -2,-0.2     2,-0.3  -0.967  29.4-107.2-148.5 151.0   39.6   27.0    4.7                A         A
   36   39 A R  E >   -H   64   0B 132     -2,-0.3     3,-1.0    28,-0.2     5,-0.3  -0.617  33.1-118.1 -81.0 137.0   41.8   27.1    1.6                A         A
   37   40 A V  T 3  S+     0   0   50     26,-1.9    26,-0.3    -2,-0.3    -1,-0.1  -0.557  96.6  28.0 -70.7 132.4   40.2   26.2   -1.7                A         A
   38   41 A G  T 3  S+     0   0   62     -2,-0.2   236,-0.5     4,-0.0    -1,-0.2   0.215  90.8  97.2 107.3  -8.5   41.8   23.1   -3.2                A         A
   39   42 A D  S <  S-     0   0   43     -3,-1.0   236,-0.8   234,-0.2    -2,-0.1   0.210  84.2-123.3-106.1  14.2   43.1   21.3   -0.1                A         A
   40   43 A G  S    S+     0   0    7    234,-0.2    -3,-0.1     1,-0.1   233,-0.1   0.365  84.4 104.8  71.2  -8.7   40.4   18.7    0.5                A         A
   41   44 A T        +     0   0    7     -5,-0.3    32,-0.5    31,-0.1     2,-0.4   0.232  52.8 100.2 -95.3  14.3   39.6   19.9    4.0                A         A
   42   45 A Q  S    S-     0   0   11     30,-0.1     2,-1.1    -6,-0.1    28,-0.1  -0.833  75.8-124.1-101.8 132.4   36.4   21.7    3.3                A         A

"""

    def test1(self):
        a = DSSPParse('1abc', self.example_string)
        
        index = a._get_index()

        assert index == 28
        
        return

    def test2(self):
        s1 = '   42   45 A Q  S    S-     0   0   11     30,-0.1     2,-1.1    -6,-0.1    28,-0.1  -0.833  75.8-124.1-101.8 132.4   36.4   21.7    3.3                A         A'
        s2 = '   42   45 A Q  E    S-     0   0   11     30,-0.1     2,-1.1    -6,-0.1    28,-0.1  -0.833  75.8-124.1-101.8 132.4   36.4   21.7    3.3                A         A'
        s3 = '   42   45 A Q  H    S-     0   0   11     30,-0.1     2,-1.1    -6,-0.1    28,-0.1  -0.833  75.8-124.1-101.8 132.4   36.4   21.7    3.3                A         A'

        char = DSSPParse.parse_ss_char(s1)
        assert char == 'L'
        
        char = DSSPParse.parse_ss_char(s2)
        assert char == 'E'

        char = DSSPParse.parse_ss_char(s3)
        assert char == 'H' 

        return


if __name__ == '__main__':
    
    cmd = load_args()
    print(cmd)
    main(**vars(cmd))
