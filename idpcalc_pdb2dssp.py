"""
*****************************
IDPCalculator: PDB to DSSP extractor
*****************************

2019
Dr. Julie Forman-Kay Lab
http://abragam.med.utoronto.ca/~JFKlab/

version: 0.1

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

version = '0.1'

LOG_NAME = 'pdb_dssp_exec.log'
ERRORLOG = 'dssp_error.log'

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

FORMATTER = logging.Formatter('%(message)s')

ch.setFormatter(FORMATTER)
log.addHandler(ch)


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
            self.data[index_where_data_starts],
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


@contextlib.contextmanager
def register_if_error():
    try:
        yield
    except Exception as e:
        log.error(traceback.format_exc())
        log.error(repr(e))
        return


def _initiate_logfile():
    
    logfile = logging.FileHandler(LOG_NAME, mode='w')
    logfile.setLevel(logging.DEBUG)
    logfile.setFormatter(FORMATTER)
    log.addHandler(logfile)
    
    errorlogfile = logging.FileHandler(ERRORLOG, mode='w')
    errorlogfile.setLevel(logging.ERROR)
    errorlogfile.setFormatter(FORMATTER)
    log.addHandler(errorlogfile)


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
        with register_if_error():
            dssp_parser.get_secondary_structure()
            ssdb.add((dssp_parser.pdbid, dssp_parser.ssd))

    ssdb.write(destination=output)

    return


if __name__ == '__main__':
    
    cmd = load_args()
    print(cmd)
    main(**vars(cmd))
