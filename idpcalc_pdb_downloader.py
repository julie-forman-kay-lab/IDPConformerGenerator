"""
*****************************
IDPCalculator: PDB Downloader
*****************************

2019
Dr. Julie Forman-Kay Lab
http://abragam.med.utoronto.ca/~JFKlab/

version: 0.1

DESCRIPTION:

    Download structural information from rcsb.org (PDB or CIF files)
    saving it in separate chain specific files.

USAGE:
    
    A text file containing a list of PDB/mmCIF IDs and chains to download.
    The list of PDB chains to download has to have the following format:
    
    XXXX
    XXXXC
    
    where XXXX is the PDB identifier and C (optional) is the specific
    chain to download. If no chain (C) is provided, saves all the chains
    from the PDB to separate files, otherwise saves only the given C.
    To specific multiple chains, repeat the PDB entry as follows:
    
    1ABCD
    1ABCE
    1ABCF
    
    With this example chains D, E and F from PDB 1ABC will be saved to disk.

HELP:
    
    use the following command for help:
    
    $ python idpcalc_pdb_downloader.py -h

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
from abc import ABC, abstractmethod
import argparse
import contextlib
import functools
import glob
import logging
from multiprocessing.pool import ThreadPool
import numpy as np
import os
from pathlib import Path as _Path
import re
import string
import time
import urllib.request

version = '0.1'

PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]
NOT_DOWNLOADED = 'PDBs_not_downloaded.txt'

CIF_ATOM_KEYS = [
    '_atom_site.group_PDB',  # 1 atom HETATM
    '_atom_site.id',  # 2 atom number
    '_atom_site.label_atom_id',  # 3 ATOM NAME
    '_atom_site.label_alt_id',  # altloc
    '_atom_site.label_comp_id',  # 4 resname
    '_atom_site.label_asym_id',  # 5 chainid
    '_atom_site.label_seq_id',  # 6 resnum
    '_atom_site.label_entity_id',
    '_atom_site.auth_atom_id',  # 3 ATOM NAME
    '_atom_site.auth_alt_id',  # altloc
    '_atom_site.auth_comp_id',  # 4 resname
    '_atom_site.auth_asym_id',  # 5 chainid
    '_atom_site.auth_seq_id',  # 6 resnum
    '_atom_site.auth_entity_id',
    '_atom_site.pdbx_PDB_ins_code',  # iCode
    '_atom_site.Cartn_x',  # 7
    '_atom_site.Cartn_y',  # 8
    '_atom_site.Cartn_z',  # 9
    '_atom_site.occupancy',  # 10
    '_atom_site.B_iso_or_equiv',  # 11
    '_atom_site.type_symbol',  # 12
    '_atom_site.pdbx_formal_charge',
    '_atom_site.pdbx_PDB_model_num',
    '_atom_site.group_PDB',
    ]


class Path(type(_Path())):
    def str(self):
        return os.fspath(self)


LOG_NAME = 'idpcalc_pdb_downloader.log'

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

formatter = logging.Formatter('%(message)s')

ch.setFormatter(formatter)

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


def record_time(process_name='', *args, **kwargs):
    """
    Decorator to record time of function execution.
    """
    def decorator(func):
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            
            start = time.time()
            
            result = func(*args, **kwargs)
            
            log.info(S(f'elapsed time :{process_name}: {time.time() - start}'))
            
            return result
        return wrapper
    return decorator


def glob_folder(folder, ext):
    ext = f"*.{ext.lstrip('*').lstrip('.')}"
    files = sorted(glob.glob(Path(folder, ext).str()))
    log.debug(f'folder {folder} read {len(files)} files with extension {ext}')
    return files


class PDBIDFactoryError(Exception):
    pass


class DownloadFailedError(Exception):
    pass


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


@functools.total_ordering
class PDBID:
    """
    PDB object identifier.
    
    Identifies unique downloadable/stored units.
    
    In the current implmentation each unit is one PDB chain,
    which is identified by the PDBID and its chain identifier.
    
    Parameters
    ----------
    name : obj:`str`
        The PDBID, for example: 1ABC
    
    chain : obj:`str`
        The chain identifier.
        Defaults to None.
    
    Attributes
    ----------
    name:
        The four character PDB identifier.
    
    chain:
        The chain identifier.
    """
    def __init__(self, name, chain=None):
        self.name = name.upper()
        self.chain = chain
        
        # made manual to completely control order
        ids = {
            'chain': chain,
            }
        self.identifiers = {}
        for name, identifier in ids.items():
            if identifier:
                self.identifiers[name] = identifier
    
    def __repr__(self):
        
        iditems = self.identifiers.items()
        
        kwargs = ', '.join(f'{key}={val!r}' for key, val in iditems)
        
        if kwargs:
            kwargs = ', ' + kwargs
        
        return '{}(name={!r}{})'.format(
            self.__class__.__name__,
            self.name,
            kwargs,
            )
    
    def __lt__(self, other):
        return str(self) < str(other)
    
    def __hash__(self):
        return hash(str(self))
    
    def __eq__(self, other):
        return str(self) == str(other)
    
    def __str__(self):
        name = f'{self.name}'
        ids = '_'.join(self.identifiers.values())
        
        if ids:
            return f'{name}_' + ids
        else:
            return name


class PDBIDFactory:
    """
    Parses input for PDBID instatiation.
    
    Returns a PDBID object.
    
    (okay... this might not be exactly a factory,
    but I wanted to have a object instantiation interface type
    instead of function call)
    """
    cull_onlypdbid_regex_single = re.compile(r'^[0-9a-zA-Z]{4}$')
    cull_onlypdbid_regex = re.compile(r'^[0-9a-zA-Z]{4}[\s\n]')
    cull_with_chains_regex_single = re.compile(r'^[0-9a-zA-Z]{5,}$')
    cull_with_chains_regex = re.compile(r'^[0-9a-zA-Z]{5,}[\s\n]')
    pdb_filename_regex = re.compile(r'[0-9a-zA-Z]{4}_[0-9a-zA-Z]{1,}\.pdb$')
    
    def __new__(cls, name):
       
        if isinstance(name, PDBID):
            return name
        
        pdb_filename_regex = {
            cls.cull_onlypdbid_regex_single: cls._parse_cull_onlypdbid,
            cls.cull_onlypdbid_regex: cls._parse_cull_onlypdbid,
            cls.cull_with_chains_regex_single: cls._parse_cull_with_chains,
            cls.cull_with_chains_regex: cls._parse_cull_with_chains,
            cls.pdb_filename_regex: cls._parse_filename_with_chain,
            }
        
        for regex, parser in pdb_filename_regex.items():
            if regex.search(str(name)):  # in case Path obj
                break
        else:
            emsg = f"PDB code format not valid: {name}. No regex matches."
            raise PDBIDFactoryError(emsg)
        
        return PDBID(*parser(name))
    
    @staticmethod
    def _parse_cull_onlypdbid(pdbid):
        return pdbid[:4], None
    
    @staticmethod
    def _parse_cull_with_chains(pdbid):
        pdbinfo = pdbid.split()[0]
        return pdbinfo[:4], pdbinfo[4:]
    
    @staticmethod
    def _parse_filename_with_chain(pdbid):
        pdbid, chainid = Path(pdbid).stem.split('_')
        return pdbid, chainid


class PDBList:
    """
    List of PDBID objects.
    
    Parameters
    ----------
    pdb_names : obj:iterator
        An iterator containing the PDB names.
        PDB names can be in the form accepted by PDBIDFactory or
        PDBID objects.
    """
    def __new__(cls, pdb_names):
        
        try:
            if isinstance(pdb_names, cls):
                return pdb_names
            else:
                return super().__new__(cls)
        except IndexError:
            return super().__new__(cls)
        
    def __init__(self, pdb_names):
        self.set = set(PDBIDFactory(element) for element in pdb_names)
    
    def __repr__(self):
        return '{}(\n    {})\n'.format(
            self.__class__.__name__,
            ',\n    '.join(repr(x) for x in self),
            )
    
    def __str__(self):
        return '{} with {} elements'.format(
            self.__class__.__name__,
            len(self),
            )
    
    def __eq__(self, other):
        try:
            return self.set == other.set
        except AttributeError:
            return self.set == other
    
    def __iter__(self):
        return iter(self.to_tuple())
    
    def __getitem__(self, index):
        return self.to_tuple()[index]
    
    def __len__(self):
        return len(self.set)
    
    def to_tuple(self):
        """
        Converts PDBList to sorted tuple.
        """
        return tuple(sorted(self.set))
    
    def difference(self, other):
        """
        Returns a PDBList that's the set difference between
        self and other.
        """
        return PDBList(tuple(self.set.difference(other.set)))
    
    def write(self, filename='PDBIDs.list'):
        """
        Writes to filename the PDBIDs in the PDBList.
        """
        with open(filename, 'w') as fh:
            fh.write('\n'.join(str(pdbid) for pdbid in self.to_tuple()))
        
        log.info(S(f'PDBIDs written to {filename}'))
    

class PDBDownloader:
    def __init__(
            self,
            pdb_list,
            destination,
            report_file=NOT_DOWNLOADED,
            ncores=1,
            record_name=('ATOM',),
            **kwargs,
            ):
        """
        
        destination (folder)
        """
        
        self.pdb_list = pdb_list

        self.destination = Path(destination)
        self.destination.mkdir(parents=True, exist_ok=True)
        
        self.report_file = report_file
        
        self.ncores = ncores
        
        self.record_name = record_name
        
        self.kwargs = kwargs
        
        self.prepare_pdb_list()
        
        return
    
    def prepare_pdb_list(self):
        
        self.pdbs_to_download = {}
        for pdbid in self.pdb_list:
            pdbentry = self.pdbs_to_download.setdefault(pdbid.name, [])
            pdbentry.append(pdbid.chain)
    
    @record_time()
    def run(self):
        log.info(T('starting raw PDB download'))
        self._create_report_file()
        log.info(S(
            f'{len(self.pdbs_to_download)} PDBs will be downloaded '
            f'and at least {len(self.pdb_list)} chains will be saved'
            ))
        
        results = ThreadPool(self.ncores).imap_unordered(
            self._download_single_pdb,
            self.pdbs_to_download.items(),
            )
        
        for r in results:
            continue
    
    def _create_report_file(self):
        self.not_downloaded = Path(self.destination, self.report_file)
        with contextlib.suppress(FileNotFoundError):
            self.not_downloaded.unlink()
    
    def _download_single_pdb(self, pdbid_and_chains_tuple):
        
        pdbname = pdbid_and_chains_tuple[0]
        chains = pdbid_and_chains_tuple[1]
        
        possible_links = [l.format(pdbname) for l in POSSIBLELINKS]
        
        with self._attempt_download(pdbname):
            response = self._download_data(possible_links)
        
        try:
            downloaded_data = response.read()
        except (AttributeError, UnboundLocalError):  # response is None
            return
        
        pdbdata = PDBDataConstructor(downloaded_data)
        
        pdbdata.build()
        
        if chains[0] is None:
            chains = pdbdata.chain_set
        
        for chain in chains:
            
            pdbdata.add_filter_record_name(self.record_name)
            pdbdata.add_filter_chain(chain)
            pdbdata.write(Path(self.destination, f'{pdbname}_{chain}.pdb'))
            pdbdata.clear_filters()
        
    def _download_data(self, possible_links):
        for weblink in possible_links:
            try:
                response = urllib.request.urlopen(weblink)
            except urllib.error.HTTPError:
                log.info(S(f'failed from {weblink}'))
                continue
            else:
                log.info(S(f'completed from {weblink}'))
                return response
        else:
            raise DownloadFailedError
    
    @contextlib.contextmanager
    def _attempt_download(self, pdbname):
        try:
            yield
        except DownloadFailedError:
            with open(self.not_downloaded, 'a') as fh:
                fh.write(f'{pdbname}\n')
            log.info(S(f'FAILED {pdbname}, see {self.not_downloaded}'))
            return


class PDBDataConstructor:
    
    def __new__(cls, data):
        
        if isinstance(data, bytes):
            data = data.decode('utf_8')
        elif isinstance(data, str):
            pass
        else:
            raise NotImplementedError('str or bytes allowed')
        
        loop_ = re.compile('[lL][oO][oO][pP]_')
        if loop_.search(data):
            return DataFromCIF(data)
        else:
            return DataFromPDB(data)


class PDBData(ABC):
    
    def __init__(self):
        self.clear_filters()
    
    def clear_filters(self):
        """
        Clears up all the previously defined filters.
        """
        self.filters = []
    
    @abstractmethod
    def build(self):
        """
        Builds the required data and attributes from the raw input.
        """
        return
    
    @property
    @abstractmethod
    def chain_set(self):
        """
        Returns a set of the different chain identifiers contained
        in the PDB data.
        """
        return
    
    @property
    @abstractmethod
    def filtered(self):
        return
    
    @abstractmethod
    def add_filter_record_name(self, record_name):
        return
    
    @abstractmethod
    def add_filter_chain(self, chain):
        return
    
    def write(self, filename):
        """
        Writes filtered data to file in PDB format.
        """
        with open(filename, 'w') as fh:
            fh.write('\n'.join(self.filtered) + '\n')
        log.info(S(f'saved: {filename}'))


class DataFromPDB(PDBData):
    
    def __init__(self, data):
        self.rawdata = data
        super().__init__()
    
    def build(self):
        self.data = self.rawdata.split('\n')
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
    
    @property
    def chain_set(self):
        return set(self.pdb_array_data[:, 5])
        
    @property
    def filtered(self):
        filtered_data = self.data
        for f in self.filters:
            filtered_data = filter(f, filtered_data)
        return filtered_data
    
    def add_filter_record_name(self, record_name):
        self.filters.append(lambda x: x.startswith(record_name))
    
    def add_filter_chain(self, chain):
        self.filters.append(lambda x: self._filter_chain(x, chain))
    
    def _filter_chain(self, line, chain):
        try:
            return line[21] == chain
        except IndexError:
            return False


class DataFromCIF(PDBData):
    
    def __init__(self, data):
        
        self.rawdata = data
        self._atom_info_keys = CIF_ATOM_KEYS
        self.pdbdata = {}
        self._void_translation_table = str.maketrans('?.', '  ')
        super().__init__()
    
    def build(self):
        self.data = MMCIFDataParser(self.rawdata).pdb_dict
        self.read_PDB_data()
        self.all_true = np.repeat(
            True,
            len(self.pdbdata['_atom_site.group_PDB']),
            )
    
    def read_PDB_data(self):
        for key in self._atom_info_keys:
            try:
                self.pdbdata[key] = np.array(self.data[key])
            except KeyError:
                continue
    
    @property
    def chain_set(self):
        try:
            chainset = self.pdbdata['_atom_site.auth_asym_id']
        except KeyError:
            chainset = self.pdbdata['_atom_site.label_asym_id']
        
        return set(chainset)
    
    @property
    def filtered(self):
        return self.convert_to_pdb()
    
    def add_filter_chain(self, chain):
        try:
            chainids = self.pdbdata['_atom_site.auth_asym_id']
        except KeyError:
            chainids = self.pdbdata['_atom_site.label_asym_id']
        
        self.filters.append(chainids == chain)
    
    def add_filter_record_name(self, record_name):
        self.filters.append(
            np.isin(
                self.pdbdata['_atom_site.group_PDB'],
                np.array(record_name),
                ),
            )
    
    def convert_to_pdb(self):
        to_output = self.all_true
        for bool_array in self.filters:
            to_output = np.logical_and(to_output, bool_array)
        
        line_formatter = (
            "{:6s}"
            "{:5d} "
            "{:<4s}"
            "{:1s}"
            "{:3s} "
            "{:1s}"
            "{:4d}"
            "{:1s}   "
            "{:8.3f}"
            "{:8.3f}"
            "{:8.3f}"
            "{:6.2f}"
            "{:6.2f}      "
            "{:<4s}"
            "{:<2s}"
            "{:2s}"
            )
        
        pdb_lines = []
        serial_count = 0
        for i in range(len(self.all_true)):
            if not to_output[i]:
                continue
            
            serial_count += 1
            
            serial = serial_count
            chainid = "0"
            segid = chainid
            
            record,\
                atname,\
                altloc,\
                resname,\
                resnum,\
                icode,\
                x,\
                y,\
                z,\
                occ,\
                bfactor,\
                element,\
                charge = self._get_line_elements(i)
        
            line = line_formatter.format(
                record,
                serial,
                atname,
                altloc,
                resname,
                chainid,
                int(resnum),
                icode,
                float(x),
                float(y),
                float(z),
                float(occ),
                float(bfactor),
                segid,
                element,
                charge,
                )
            
            pdb_lines.append(line)
            
        return pdb_lines
    
    def _get_line_elements(self, i):
        
        record = self.pdbdata.get('_atom_site.group_PDB')[i]
        try:
            atname = self.pdbdata['_atom_site.auth_atom_id'][i]
        except KeyError:
            atname = self.pdbdata['_atom_site.label_atom_id'][i]
        
        try:
            altloc = self.pdbdata['_atom_site.auth_alt_id'][i]
        except KeyError:
            altloc = self.pdbdata['_atom_site.label_alt_id'][i]
        altloc = altloc.translate(self._void_translation_table)
        
        try:
            resname = self.pdbdata['_atom_site.auth_comp_id'][i]
        except KeyError:
            resname = self.pdbdata['_atom_site.label_comp_id'][i]
        
        try:
            resnum = self.pdbdata['_atom_site.auth_seq_id'][i]
        except KeyError:
            resnum = self.pdbdata['_atom_site.label_seq_id'][i]
        
        try:
            icode = self.pdbdata['_atom_site.pdbx_PDB_ins_code'][i]
        except KeyError:
            icode = " "
        icode = icode.translate(self._void_translation_table)
        
        x = self.pdbdata.get('_atom_site.Cartn_x')[i]
        y = self.pdbdata.get('_atom_site.Cartn_y')[i]
        z = self.pdbdata.get('_atom_site.Cartn_z')[i]
        
        occ = self.pdbdata.get('_atom_site.occupancy')[i]
        bfactor = self.pdbdata.get('_atom_site.B_iso_or_equiv')[i]
        
        element = self.pdbdata.get('_atom_site.type_symbol')[i]
        charge = self.pdbdata.get('_atom_site.pdbx_formal_charge')[i]
        
        return record, atname, altloc, resname, resnum, icode, x, y, z,\
            occ, bfactor, element, charge


class MMCIFDataParser:
    def __init__(self, data):
        self.data = data.split('\n')
        self.cifatomkeys = CIF_ATOM_KEYS
        self.read()
    
    def read(self):
        
        self.pdb_dict = {}
        found = False
        for ii, line in enumerate(self.data):
            if line.startswith('_atom_site.'):
                found = True
                self.pdb_dict.setdefault(line.strip(), [])
            elif found:
                atom_start_index = ii
                break
        
        for line in self.data[atom_start_index:]:
            
            if line.startswith('#'):
                break
            ls = line.split()
            
            for i, key in enumerate(self.pdb_dict.keys()):
                self.pdb_dict[key].append(ls[i])


def load_args():
    
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=__doc__,
        )
    
    ap.add_argument(
        'pdb_list',
        help="PDB curated list",
        type=Path,
        )
    
    ap.add_argument(
        '-d',
        '--destination',
        help=(
            "The folder to where PDBs will be downloader. "
            "If none provided uses the current working directory"
            ),
        type=Path,
        default=Path.cwd(),
        )
    
    ap.add_argument(
        '-u',
        '--update',
        help='Updates destination folder according to input PDB list',
        action='store_true',
        )
    
    ap.add_argument(
        '-rn',
        '--record_name',
        help="The coordinate lines record name",
        default=('ATOM',),
        type=tuple,
        )
    
    ap.add_argument(
        '-n',
        '--ncores',
        help="The number of cores to use during the download",
        type=int,
        default=1,
        )
    
    cmd = ap.parse_args()
    
    return cmd


def main(
        pdb_list,
        destination=None,
        update=False,
        record_name=('ATOM',),
        ncores=1,
        ):
    
    # initiates log file only if main is run
    logfile = logging.FileHandler(LOG_NAME, mode='w')
    logfile.setLevel(logging.DEBUG)
    log.addHandler(logfile)
    logfile.setFormatter(formatter)
    
    with pdb_list.open('r') as fh:
        pdblist = PDBList(fh.readlines())
    
    log.info(T('reading input PDB list'))
    log.info(S(f'from: {pdb_list}'))
    log.info(S(f'{str(pdblist)}'))
    log.info(S('done\n'))
    
    if destination:
        pdblist_destination = PDBList(glob_folder(destination, '*.pdb'))
        log.info(T('reading destination folder'))
        log.info(S(f'from: {destination}'))
        log.info(S(f'{str(pdblist_destination)}'))
        log.info(S('done\n'))
        
        pdblist_comparison = pdblist.difference(pdblist_destination)
        log.info(T('Comparison between input and destination'))
        log.info(S(f'{str(pdblist_comparison)}'))
        log.info(S('done\n'))
    
    if update:
        pdbdownloader = PDBDownloader(
            pdblist_comparison,
            destination,
            record_name=record_name,
            ncores=ncores,
            )
        pdbdownloader.prepare_pdb_list()
        pdbdownloader.run()
        pdblist_updated = PDBList(glob_folder(destination, '*.pdb'))
        log.info(T('Reading UPDATED destination'))
        log.info(S(f'{str(pdblist_updated)}'))
        log.info(S('done\n'))
    
    return


class TestPDBID:
    
    import pytest
    
    def test_PDBID1(self):
        pdbid = PDBID('1ABC', chain='D')
        assert pdbid.name == '1ABC'
    
    def test_PDBID_null_chain(self):
        pdbid = PDBID('1ABC')
        assert pdbid.name == '1ABC'
        assert pdbid.chain is None
    
    def test_PDBID_str(self):
        pdbid = PDBID('1ABC', chain='Z')
        assert str(pdbid) == '1ABC_Z'
        
    def test_PDBID_str2(self):
        pdbid = PDBID('1ABC')
        assert str(pdbid) == '1ABC'
    
    def test_repr(self):
        pdbid = PDBID('1ABC', chain='Z')
        assert repr(pdbid) == "PDBID(name='1ABC', chain='Z')"
    
    def test_repr2(self):
        pdbid = PDBID('1ABC')
        assert repr(pdbid) == "PDBID(name='1ABC')"
    
    def test_equality(self):
        pdbid1 = PDBID('1ABC', chain='Z')
        pdbid2 = PDBID('1ABC', chain='Z')
        assert pdbid1 == pdbid2
    
    def test_lower(self):
        pdbid1 = PDBID('1ABC', chain='X')
        pdbid2 = PDBID('1ABC', chain='Z')
        assert pdbid1 < pdbid2
    
    def test_higher(self):
        pdbid1 = PDBID('1ABC', chain='Z')
        pdbid2 = PDBID('1ABC', chain='X')
        assert pdbid1 > pdbid2


class TestPDBIDFactory:
    
    import pytest
    
    def test_parse1(self):
        pdbid0 = PDBIDFactory('1ABC')
        pdbid1 = PDBID('1ABC')
        assert pdbid0 == pdbid1
    
    def test_parse2(self):
        pdbid0 = PDBIDFactory('1ABC ')
        pdbid1 = PDBID('1ABC')
        assert pdbid0 == pdbid1
    
    def test_parse3(self):
        pdbid0 = PDBIDFactory('1ABCX')
        pdbid1 = PDBID('1ABC', chain='X')
        assert pdbid0 == pdbid1
    
    def test_parse4(self):
        pdbid0 = PDBIDFactory('1ABCX            something')
        pdbid1 = PDBID('1ABC', chain='X')
        assert pdbid0 == pdbid1
    
    def test_parse5(self):
        pdbid0 = PDBIDFactory('1ABCXYZ')
        pdbid1 = PDBID('1ABC', chain='XYZ')
        assert pdbid0 == pdbid1
    
    def test_parse6(self):
        pdbid0 = PDBIDFactory('1ABCXYZ\n')
        pdbid1 = PDBID('1ABC', chain='XYZ')
        assert pdbid0 == pdbid1
    
    def test_parse7(self):
        pdbid0 = PDBIDFactory('some/folder/PDBs/1ABC_D.pdb')
        pdbid1 = PDBID('1ABC', chain='D')
        assert pdbid0 == pdbid1

    def test_parse8(self):
        pdbid0 = PDBIDFactory(Path('some', 'fldr', 'PDBs', '1ABC_DEF.pdb'))
        pdbid1 = PDBID('1ABC', chain='DEF')
        assert pdbid0 == pdbid1


class TestPDBList:
    
    def test_1(self):
        pdblist = PDBList(('1ABC',))
        assert pdblist[0] == PDBID('1ABC')
    
    def test_2(self):
        
        names = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            )
        
        pdblist = PDBList(names)
        pdbids = [PDBIDFactory(name) for name in names]
        assert len(pdblist) == len(names)
        assert pdblist == set(pdbids)
    
    def test_difference(self):
        names1 = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            )
        names2 = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            )
        
        pdblist1 = PDBList(names1)
        pdblist2 = PDBList(names2)
        
        pdblist3 = PDBList(names1[-1:])
        
        assert pdblist1.difference(pdblist2) == pdblist3
        
    def test_to_tuple(self):
        names_sorted = (
            '1ABC',
            '2ABCZ',
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            )
        
        names_unsorted = (
            '2ABCZ',
            Path('some', 'path', '3RFC_ZZZ.pdb'),
            str(Path('some', 'path', '3RFC_YYY.pdb')),
            '1ABC',
            )
        
        pdblist1 = PDBList(names_sorted)
        pdblist2 = PDBList(names_unsorted)
        
        assert pdblist1 == pdblist2
        assert isinstance(pdblist2.to_tuple(), tuple)
        assert tuple(pdblist1) == pdblist2.to_tuple()


class TestPDBDownloader:
    
    def test1(self):
        destination = Path('here')
        PDBDownloader(PDBList(('1ABC',)), destination)
        assert destination.exists()
        destination.rmdir()


class TestPDBDataConstructor:
    
    def test1(self):
        pdbc = PDBDataConstructor('loop_')
        assert isinstance(pdbc, DataFromCIF)
    
    def test2(self):
        pdbc = PDBDataConstructor('some str with l o o p from CIF')
        assert isinstance(pdbc, DataFromPDB)


class TestDataFromPDB:
    
    s = """REMARK some remark
ATOM      1 N    VAL 0   4     -78.583  98.529 -20.107  1.00 70.49      0   N ? 
ATOM      2 CA   VAL 0   4     -78.412  98.561 -18.649  1.00 72.59      0   C ? 
ATOM      3 C    VAL 0   4     -76.913  98.491 -18.327  1.00 76.61      0   C ? 
ATOM      4 O    VAL 0   4     -76.112  98.005 -19.146  1.00 80.22      0   O ? 
ATOM      5 CB   VAL 0   4     -79.159  97.384 -17.947  1.00 73.51      0   C ? 
ATOM      6 CG1  VAL 0   4     -79.193  97.580 -16.405  1.00 74.65      0   C ? 
ATOM      7 CG2  VAL 0   4     -80.581  97.265 -18.475  1.00 72.07      0   C ? 
ATOM      8 N    SER 1   5     -76.530  98.969 -17.143  1.00 75.44      0   N ? 
ATOM      9 CA   SER 1   5     -75.112  99.085 -16.790  1.00 74.03      0   C ? 
ATOM     10 C    SER 1   5     -74.665  98.040 -15.785  1.00 74.64      0   C ? 
END"""

    def test1(self):
        pdb = DataFromPDB(self.s)
        pdb.build()
        assert pdb.pdb_array_data[:, 5].size == 10
        assert pdb.chain_set == {'0', '1'}
        
        pdb.add_filter_chain('1')
        
        assert len(list(pdb.filtered)) == 3
    

class TestDataFromCIF:
    
    s = """
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.type_symbol
_atom_site.pdbx_formal_charge
ATOM      1 N   . VAL 0   4     -78.583  98.529 -20.107  1.00 70.49       N ?
ATOM      2 CA  . VAL 0   4     -78.412  98.561 -18.649  1.00 72.59       C ?
ATOM      3 C   . VAL 0   4     -76.913  98.491 -18.327  1.00 76.61       C ?
ATOM      4 O   . VAL 0   4     -76.112  98.005 -19.146  1.00 80.22       O ?
ATOM      5 CB  . VAL 0   4     -79.159  97.384 -17.947  1.00 73.51       C ?
ATOM      6 CG1 . VAL 0   4     -79.193  97.580 -16.405  1.00 74.65       C ?
ATOM      7 CG1 . VAL 0   4     -79.193  97.580 -16.405  1.00 74.65       C ?
ATOM      8 N   . SER 1   5     -76.530  98.969 -17.143  1.00 75.44       N ?
ATOM      9 CA  . SER 1   5     -75.112  99.085 -16.790  1.00 74.03       C ?
ATOM     10 C   . SER 1   5     -74.665  98.040 -15.785  1.00 74.64       C ?
#"""

    def test1(self):
        pdb = DataFromCIF(self.s)
        pdb.build()
        assert pdb.chain_set == {'0', '1'}
        
        pdb.add_filter_chain('1')
        
        assert len(list(pdb.filtered)) == 3


if __name__ == '__main__':
    
    cmd = load_args()
    main(**vars(cmd))
