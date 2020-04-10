"""Contain  handlers of PDB information."""
# '\n'.join(','.join(map(str, bb)) for bb in b)
import contextlib
import functools
import re
import string
import traceback
import urllib.request
from abc import ABC, abstractmethod
from multiprocessing.pool import ThreadPool

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs import libtimer
from idpconfgen.logger import S, T


PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]


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


class _PDBParams:
    """
    Namespace for `PDB format v3`_.

    .. _old PDB format: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
    """  # noqa: E501

    # slicers of the ATOM and HETATM lines
    atom_line_name = slice(0, 6)
    atom_serial = slice(6, 11)
    atom_atom_name = slice(12, 16)
    atom_altloc = slice(16, 17)
    atom_resname = slice(17, 20)
    atom_chainid = slice(21, 22)
    atom_resseq = slice(22, 26)
    atom_icode = slice(26, 27)
    atom_xcoord = slice(30, 38)
    atom_ycoord = slice(38, 46)
    atom_zcoord = slice(46, 54)
    atom_occupancy = slice(54, 60)
    atom_tempfactor = slice(60, 66)
    atom_element = slice(76, 78)
    atom_charge = slice(78, 80)

    def __setattr__(self, key, value):
        raise NotImplementedError(f'Can not set attributes to {self.__class__}')

    @property
    def atom_slicers(self):
        """Ordered list of ATOM and HETATOM slicers."""
        try:
            return self._atom_slicers
        except AttributeError:
            atom_slicers = list(filter(
                lambda x: x[0].startswith('atom_'),
                self.__dict__.items(),
                ))
            self._atom_slicers = [s[1] for s in atom_slicers]

            # ensure
            assert all(
                isinstance(s, slice)
                for s in self._atom_slicer
                ))
            return self._atom_slicers

PDBParams = _PDBParams()

# this servers read_pdb_data_to_array mainly
# it is here for performance
_pdb_atom_line_headings = {
    'which': ('ATOM', 'HETATM'),
    'ATOM': 'ATOM',
    'HETATM': 'HETATM',
    }
















class PDBIDFactory:
    r"""
    Parse input for PDBID instatiation.
    
    Parameters
    ----------
    name : str or Path
        The code name or ID that identified the PDB.
        Possible formats:

            - XXXX
            - XXXXC\*
            - XXXX_C\*
            - \*.pdb

        where XXXX is the PDBID code, C is the chain ID and * means
        any number of characters. PDB and chaind ID codes are any digits,
        lower and upper case letters.

    Returns
    -------
    :class:`PDBID` object.
    
    """

    rgx_XXXX = re.compile(r'^[0-9a-zA-Z]{4}(\s|$)')
    rgx_XXXXC = re.compile(r'^[0-9a-zA-Z]{5,}(\s|$)')
    rgx_XXXX_C = re.compile(r'^[0-9a-zA-Z]{4}_[0-9a-zA-Z]+(\s|$)')
    
    def __new__(cls, name):
        """Construct class."""
        if isinstance(name, PDBID):
            return name

        namep = Path(name)
        if namep.suffix == '.pdb':
            name = namep.stem

        # where XXXX is the PDBID and C the chain ID
        pdb_filename_regex = {
            cls.rgx_XXXX: cls._parse_XXXX,
            cls.rgx_XXXXC: cls._parse_XXXXC,
            cls.rgx_XXXX_C: cls._parse_XXXX_C,
            }
        
        for regex, parser in pdb_filename_regex.items():
            if regex.search(str(name)):  # in case Path obj
                return PDBID(*parser(name))
        else:
            emsg = f"PDB code format not valid: {name}. No regex matches."
            raise EXCPTS.PDBIDFactoryError(emsg)
    
    @staticmethod
    def _parse_XXXX(pdbid):
        return pdbid[:4], None
    
    @staticmethod
    def _parse_XXXXC(pdbid):
        pdbinfo = pdbid.split()[0]
        return pdbinfo[:4], pdbinfo[4:]
    
    @staticmethod
    def _parse_XXXX_C(pdbid):
        pdbid, chainid = pdbid.split()[0].split('_')
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

    def __new__(cls, pdb_names):  # noqa: D102
        
        try:
            if isinstance(pdb_names, cls):
                return pdb_names
            else:
                return super().__new__(cls)
        except IndexError:
            return super().__new__(cls)
        
    def __init__(self, pdb_names):
        valid_pdb_names = filter(
            # str() because may receive Paths
            lambda x: not str(x).startswith('#'),
            pdb_names,
            )
        self.set = set(PDBIDFactory(element) for element in valid_pdb_names)
    
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
        """Convert PDBList to sorted tuple."""
        return tuple(sorted(self.set))
    
    def difference(self, other):
        """
        Difference between self and other.

        Returns
        -------
        PDBList
        """
        return PDBList(tuple(self.set.difference(other.set)))
    
    def write(self, filename='PDBIDs.list'):
        """
        Write to a file the PDBIDs in the PDBList.

        Parameters
        ----------
        filename : str, optional
            The output file name.
        """
        with open(filename, 'w') as fh:
            fh.write('\n'.join(str(pdbid) for pdbid in self.to_tuple()))
        
        log.info(S(f'PDBIDs written to {filename}'))
        

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


class PDBData(ABC):
    """
    Abstract base class for PDB data filters.

    PDB data filters do what the name implies, filter PDB data according
    to some criteria.

    Filters can be accumulated to add restrictions.
    """

    def __init__(self):
        self.clear_filters()

    @property
    def filters(self):
        return self._filters

    def clear_filters(self):
        """Clear up all the previously defined filters."""
        self._filters = []

    def build(self):
        """Build the required data and attributes from the raw input."""
        self._build()
        assert isinstance(self.data_array, np.ndarray)
        assert self.data_array.shape[1] == 15
        assert self.data_array.shape[0] > 0
        return

    @property
    @abstractmethod
    def chain_set(self):
        """
        Set of PDB chain IDs.

        Considers only Raw data.
        """
        return

    @property
    @abstractmethod
    def filtered(self):  # noqa: D102
        return

    def write(self, filename):
        """Write filtered data to file in PDB format."""
        data2write = self.join_filtered()

        with open(filename, 'w') as fh:
            fh.write(data2write + '\n')
        log.info(S(f'saved: {filename}'))

    def join_filtered(self, char='\n'):
        """
        Join filtered PDB data into a single string.

        Parameters
        ----------
        char : str
            The character with which perform string join.
            Defaults to `new_line`.

        Returns
        -------
        string
            A unique string resulting from the filetered data.

        Raises
        ------
        :class:`EmptyFilterError`
            If resulting string is an empty string.
        """
        datafiltered = '\n'.join(self.filtered)
        if datafiltered.strip():
            return datafiltered
        else:
            raise EXCPTS.EmptyFilterError


class PDBStructure(PDBData):
    """
    Hold structural data from PDB files.

    Parameters
    ----------
    data : str
        Raw structural data from PDB formatted files.
    """
    def __init__(self, data):
        self.rawdata = data
        super().__init__()

    def _build(self):
        """
        Read structure raw data in :attr:`rawdata`.

        After `.build()`, filters and data can be accessed.
        """
        self.data_array = read_pdb_data_to_array(
            self.rawdata.split('\n'),
            which='both',
            )

    @property
    def chain_set(self):
        """All chain IDs present in the raw dataset."""  # noqa: D401
        return set(self.pdb_array_data[:, 5])

    @property
    def filtered(self):
        """
        Parse the raw data applying the selected filters.

        Returns
        -------
        list
            The data in PDB format after filtering.
        """
        filtered_data = self.data
        for f in self.filters:
            filtered_data = filter(f, filtered_data)
        return filtered_data

    def add_filter_record_name(self, record_name):
        """Add filter for record names."""
        self.filters.append(lambda x: x.startswith(record_name))

    def add_filter_chain(self, chain):
        """Add filters for chain."""
        self.filters.append(lambda x: self._filter_chain(x, chain))

    def _filter_chain(self, line, chain):
        try:
            return line[21] == chain
        except IndexError:
            return False


class DataFromCIF(PDBData):
    """
    Reads structural data for proteins in CIF files.

    Parameters
    ----------
    data : str
        The CIf contained data.
    """

    def __init__(self, data):

        self.rawdata = data
        self._atom_info_keys = CIF_ATOM_KEYS
        self.pdbdata = {}
        self._void_translation_table = str.maketrans('?.', '  ')
        super().__init__()

    def build(self):
        """
        Build data.

        After data build, filters can be applied.
        """
        self.data = MMCIFDataParser(self.rawdata).pdb_dict
        self.read_PDB_data()
        self.all_true = np.repeat(
            True,
            len(self.pdbdata['_atom_site.group_PDB']),
            )
    
    def read_PDB_data(self):
        """Read PDB data."""
        for key in self._atom_info_keys:
            try:
                self.pdbdata[key] = np.array(self.data[key])
            except KeyError:
                continue
    
    @property
    def chain_set(self):
        """All the chain ids present in the data set."""
        try:
            chainset = self.pdbdata['_atom_site.auth_asym_id']
        except KeyError:
            chainset = self.pdbdata['_atom_site.label_asym_id']
        
        return set(chainset)
    
    @property
    def filtered(self):
        """Filtered data in PDB format."""  # noqa: D401
        return self.convert_to_pdb()
    
    def add_filter_chain(self, chain):
        """
        Add filter for chain identifier.

        Only atoms of the specified chain will be parsed.
        """
        try:
            chainids = self.pdbdata['_atom_site.auth_asym_id']
        except KeyError:
            chainids = self.pdbdata['_atom_site.label_asym_id']
        
        self.filters.append(chainids == chain)
    
    def add_filter_record_name(self, record_name):
        """
        Add a filter for the RECORD NAME.

        Usually ATOM or/and HETATM. `Read further`_
        
        .. _Read further: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
        """  # noqa: E501
        self.filters.append(
            np.isin(
                self.pdbdata['_atom_site.group_PDB'],
                np.array(record_name),
                ),
            )
    
    def convert_to_pdb(self):
        """
        Convert internal data to PDB format.

        Returns
        -------
        list
            A list of the formatted PDB lines.
        """
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
    """
    Quick mmCIF data parser for Proteins.
    
    Parameters
    ----------
    data : str
        The structural data contained in mmCIF format.
    """
    
    def __init__(self, data):
        self.data = data.split('\n')
        self.cifatomkeys = CIF_ATOM_KEYS
        self.read()
    
    def read(self):
        """
        Read data string to pdb_dict attribute.

        Performed at instantiation.
        """
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


class Structure:
    """
    Create a protein structure representation.

    Implements readers for:

    - PDB format version 3
    - mmCIF

    Returns the corresponding representation of Structure instance
    based on the data input.

    Parameters
    ----------
    data:
        PDB structural data as string or bytes.
    """

    def __new__(cls, data):
        """Factor PDB data."""
        if isinstance(data, Path):
            datastr = data.read_text()
        elif isinstance(data, bytes):
            datastr = data.decode('utf_8')
        elif isinstance(data, str):
            pass
        else:
            raise NotImplementedError('str or bytes allowed')

        loop_ = re.compile('[lL][oO][oO][pP]_')
        if loop_.search(datastr):
            return DataFromCIF(datastr)
        else:
            return PDBStructure(datastr)


class PDBDownloader:
    """
    Control PDB downloading operations.

    Given a list of :class:`PDBID` downloads those PDB files.

    Parameters
    ----------
    pdb_list : list
        List of PDBIDs.

    destination : str or Path
        Destination folder.

    record_name : tuple
        A tuple containing the atom record names to store from the
        downloaded PDBs. Record names are as described for the PDB
        format v3. Usually 'ATOM' and 'HETATM' can be used.
        Defaults to ('ATOM',)

    ncores : int, optional
        Number of cores to use during the download phase.
        All downloads operations are stored in a pool, each core
        performs a download operation per time grabing those operations
        from the pool. The more cores the more operations are performed
        simultaneoursly.
        Defaults to 1.
    """

    def __init__(
            self,
            pdb_list,
            destination,
            ncores=1,
            record_name=('ATOM',),
            **kwargs,
            ):

        self.pdb_list = pdb_list

        self.destination = Path(destination)
        self.destination.mkdir(parents=True, exist_ok=True)

        self.ncores = ncores

        self.record_name = record_name

        self.kwargs = kwargs

        self.prepare_pdb_list()

        return

    def prepare_pdb_list(self):
        """Prepare a list with the PDBs to download."""
        self.pdbs_to_download = {}
        for pdbid in self.pdb_list:
            pdbentry = self.pdbs_to_download.setdefault(pdbid.name, [])
            pdbentry.append(pdbid.chain)
    
    @libtimer.record_time()
    def run(self):
        """Run download operation."""
        log.info(T('starting raw PDB download'))
        log.info(S(
            f'{len(self.pdbs_to_download)} PDBs will be downloaded '
            f'and at least {len(self.pdb_list)} chains will be saved'
            ))
        
        results = ThreadPool(self.ncores).imap_unordered(
            self._download_single_pdb,
            self.pdbs_to_download.items(),
            )
        
        for _r in results:
            continue
    
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
        
        pdbdata = Structure(downloaded_data)
        
        pdbdata.build()
        
        if chains[0] is None:
            chains = pdbdata.chain_set
        
        for chain in chains:
            
            pdbdata.add_filter_record_name(self.record_name)
            pdbdata.add_filter_chain(chain)
            destination = Path(self.destination, f'{pdbname}_{chain}.pdb')
            try:
                pdbdata.write(destination)
            except EXCPTS.EmptyFilterError:
                log.error(T(f'Empty Filter for {destination}'))
                log.debug(traceback.format_exc())
                log.error(S(f'record_name: {self.record_name}'))
                log.error(S(f'chain filter: {chain}'))
            pdbdata.clear_filters()
        
    def _download_data(self, possible_links):
        for weblink in possible_links:
            try:
                response = urllib.request.urlopen(weblink)
            except urllib.error.HTTPError:
                log.error(S(f'failed from {weblink}'))
                continue
            else:
                log.info(S(f'completed from {weblink}'))
                return response
        else:
            raise EXCPTS.DownloadFailedError
    
    @contextlib.contextmanager
    def _attempt_download(self, pdbname):
        try:
            yield
        except EXCPTS.DownloadFailedError as e:
            log.error(S(f'{repr(e)}: FAILED {pdbname}'))
            return


def read_pdb_file_to_array(lines, which='both'):
    """
    Transform PDB data into an array.

    Parameters
    ----------
    lines : list of strs
        List of the PDB format v3 file lines.

    which : str
        Which lines to consider ['ATOM', 'HETATM', 'both'].
        Defaults to `'both'`, considers both 'ATOM' and 'HETATM'.

    Returns
    -------
    numpy.ndarray of (N, 15)
        Where N are the number of ATOM and/or HETATM lines,
        and 15 the number of fields in ATOM/HETATM lines according
        to the PDB format v3.
    """
    coords_headings = _pdb_atom_line_headings

    try:
        atom_hetatm_lines = filter(
            lambda x: x.startswith(coord_headings[which]),
            lines,
            )
    except KeyError:
        err = ValueError(f'`which` got an unexpected value \'{which}\''.)
        raise err from None

    pdb_data_array = np.empty(
        (len(atom_hetatm_lines), len(PDBParams.atom_slicers)),
        dtype='<U8',
        )

    for row, line in enumerate(atom_hetatm_lines):
        for column, slicer_item in enumerate(PDBParams.atom_slicers):
            pdb_data_array[row, column] = line[slicer_item]

    return pdb_data_array

