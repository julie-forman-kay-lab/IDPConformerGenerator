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
    atom_record = slice(0, 6)
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
    atom_segid = slice(72, 76)
    atom_element = slice(76, 78)
    atom_charge = slice(78, 80)

    _atom_attrs = filter(
        lambda x: x[0].startswith('atom_'),
        self.__dict__.items(),
        )

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

    def __setattr__(self, key, value):
        raise NotImplementedError(f'Can not set attributes to {self.__class__}')

    @property
    def atom_slicers(self):
        """Ordered list of ATOM and HETATOM slicers."""
        try:
            return self._atom_slicers
        except AttributeError:
            self._atom_slicers = [s[1] for s in self._atom_attrs]

            # ensure
            assert all(
                isinstance(s, slice)
                for s in self._atom_slicer
                ))
            return self._atom_slicers

    def acol(self):
        try:
            return self._acol
        except AttributeError:
            ac = namedtuple(
                'AtomCols',
                (s[0].lstrip("atom_") for s in self._atom_attrs),
                )
            self._acol(*range(len(names)))
            return self._acol

PDBParams = _PDBParams()

# this servers read_pdb_data_to_array mainly
# it is here for performance
_pdb_atom_line_headings = {
    'which': ('ATOM', 'HETATM'),
    'ATOM': 'ATOM',
    'HETATM': 'HETATM',
    }


def is_cif(datastr):
    """Detect if `datastr` is a CIF file."""
    assert isinstance(datastr, str)
    cif_loop = re.compile('[lL][oO][oO][pP]_')
    return bool(cif_loop.search(datastr)


def is_pdb(datastr):
    """Detect if `datastr` if a PDB format v3 file."""
    return bool(datastr.count('\nATOM ') > 0)


def parse_pdb_to_array(datastr, which='both'):
    """
    Transform PDB data into an array.

    Parameters
    ----------
    datastr : str
        String representing the PDB format v3 file.

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
    lines = datastr.split('\n')

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


def parse_cif_to_array(datastr):
    """
    """
    pdb_dict = {}
    data = datastr.split('\n')
    found = False
    for ii, line in enumerate(self.data):
        if line.startswith('_atom_site.'):
            found = True
            pdb_dict.setdefault(line.strip(), [])
        elif found:
            atom_start_index = ii
            break

    for line in data[atom_start_index:]:
        if line.startswith('#'):
            break
        ls = line.split()

        for i, key in enumerate(pdb_dict.keys()):
            pdb_dict[key].append(ls[i])
    else:
        number_of_atoms = len(ls)

    data_array = np.empty(
        (len(number_of_atoms), len(PDBParams.atom_slicers)),
        dtype='<U8',
        )

    for ii in range(number_of_atoms):
        data_array[ii, :] = cif_get_line_elements(ii)
    return data_array


def cif_get_line_elements(self, i):
    """
    """
    # http://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html
    record = self.pdbdata.get('_atom_site.group_PDB')[i]

    serial = pdbdata.get('_atom_site.Id')[i]

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
        chainids = self.pdbdata['_atom_site.auth_asym_id']
    except KeyError:
        chainids = self.pdbdata['_atom_site.label_asym_id']


    try:
        resseq = self.pdbdata['_atom_site.auth_seq_id'][i]
    except KeyError:
        resseq = self.pdbdata['_atom_site.label_seq_id'][i]


    try:
        icode = self.pdbdata['_atom_site.pdbx_PDB_ins_code'][i]
    except KeyError:
        icode = " "
    icode = icode.translate(self._void_translation_table)


    x = self.pdbdata.get('_atom_site.Cartn_x')[i]
    y = self.pdbdata.get('_atom_site.Cartn_y')[i]
    z = self.pdbdata.get('_atom_site.Cartn_z')[i]
    occ = self.pdbdata.get('_atom_site.occupancy')[i]
    tempfactor = self.pdbdata.get('_atom_site.B_iso_or_equiv')[i]
    element = self.pdbdata.get('_atom_site.type_symbol')[i]
    charge = self.pdbdata.get('_atom_site.pdbx_formal_charge')[i]

    return (
        record,
        serial,
        atname,
        altloc,
        resname,
        chainids,
        resseq,
        icode,
        x,
        y,
        z,
        occ,
        tempfactor,
        ' ',  # segid
        element,
        charge,
        )


# order matters
structure_parsers = [
    (is_cif, parse_cif_to_array),
    (is_pdb, parse_pdb_to_array),
    ]


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


class Structure:
    """
    Hold structural data from PDB files.

    Parameters
    ----------
    data : str
        Raw structural data from PDB formatted files.
    """
    # all should return a string
    _data_input_types = {
        type(Path()): lambda x: x.read_text(),
        bytes: lambda x: x.decode('utf_8'),
        str: lambda x: x,
        }


    def __init__(self, data):
        data_type = type(data)

        try:
            datastr = self._data_input_types[data_type](data)
        except KeyError as err:
            err2 = NotImplementedError('Struture data not of proper type')
            raise err2 from err
        assert isinstance(datastr, str)

        for condition, parser in structure_parsers.items():
            if condition(datastr):
                self._structure_parser = parser
                break

        self.datastr = datastr
        self.data_array = None
        self.clear_filters()
        assert isinstance(self.filters, list)

    @property
    def filters(self):
        return = self._filters

    @property
    def filtered_atom(self):
        """
        Filter data array by the selected filters.

        Returns
        -------
        list
            The data in PDB format after filtering.
        """
        filtered_data = self.data_array
        for f in self.filters:
            filtered_data = filter(f, filtered_data)
        return filtered_data

    @property
    def chain_set(self):
        """All chain IDs present in the raw dataset."""  # noqa: D401
        return set(self.pdb_array_data[:, PDBParams.acol.chainid])

    def clear_filters(self):
        self._filters = []

    def pop_last_filter(self):
        self._filters.pop()

    def build(self):
        """
        Read structure raw data in :attr:`rawdata`.

        After `.build()`, filters and data can be accessed.
        """
        self.data_array = self._structure_parser(self.datastr)

    def add_filter(self, function):
        """Adds a function as filter."""
        self.filters.append(function)

    def add_filter_record_name(self, record_name):
        """Add filter for record names."""
        self.filters.append(
            lambda x: x[PDBParams.acol.record].startswith(record_name)
            )

    def add_filter_chain(self, chain):
        """Add filters for chain."""
        self.filters.append(
            lambda x: x[PDBParams.acol.chainid].startswith(record_name)
            )

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


