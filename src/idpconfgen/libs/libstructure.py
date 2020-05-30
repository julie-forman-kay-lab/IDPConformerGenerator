"""
Store internal protein structure representation.

Classes
-------
Structure
    The main API that represents a protein structure in IDPConfGen.
"""
import warnings
from collections import defaultdict
from functools import reduce

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import definitions as DEFS
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs import libpdb
from idpconfgen.libs.libcif import CIFParser, is_cif
from idpconfgen.logger import S


# module variables are defined at the end.
class Structure:
    """
    Hold structural data from PDB/mmCIF files.

    Run the ``.buil()`` method to read the structure.

    Cases for PDB Files:
    * If there are several MODELS only the first model is considered.

    Parameters
    ----------
    data : str, bytes, Path
        Raw structural data from PDB/mmCIF formatted files.

    Examples
    --------
    Opens a PDB file, selects only chain 'A' and saves selection to a file.
    >>> s = Structure('1ABC.pdb')
    >>> s.build()
    >>> s.add_filter_chain('A')
    >>> s.write_PDB('out.pdb')

    Opens a mmCIF file, selects only residues above 50 and saves
    selection to a file.
    >>> s = Structure('1ABC.cif')
    >>> s.build()
    >>> s.add_filter(lambda x: int(x[col_resSeq]) > 50)
    >>> s.write_PDB('out.pdb')
    """

    __slots__ = [
        '_data_array',
        '_datastr',
        '_filters',
        '_structure_parser',
        'kwargs',
        ]

    def __init__(self, data, **kwargs):

        datastr = get_datastr(data)
        self._structure_parser = detect_structure_type(datastr)

        self._datastr = datastr
        self.kwargs = kwargs
        self.clear_filters()
        assert isinstance(self.filters, list)

    def build(self):
        """
        Read structure raw data in :attr:`rawdata`.

        After `.build()`, filters and data can be accessed.
        """
        self._data_array = self._structure_parser(self._datastr, **self.kwargs)
        del self._datastr

    def clear_filters(self):
        """Clear/Deletes registered filters."""
        self._filters = []

    @property
    def data_array(self):
        """Contain structure data in the form of a Numpy array."""
        try:
            return self._data_array
        except AttributeError as err:
            errmsg = (
                'Please `.build()` the Structure before attempting to access'
                'its data'
                )
            raise EXCPTS.NotBuiltError(errmsg=errmsg) from err

    @property
    def filters(self):
        """Filter functions registered ordered by registry record."""
        return self._filters

    @property
    def filtered_atoms(self):
        """
        Filter data array by the selected filters.

        Returns
        -------
        list
            The data in PDB format after filtering.
        """
        apply_filter = _APPLY_FILTER
        return np.array(list(reduce(
            apply_filter,
            self.filters,
            self.data_array,
            )))

    @property
    def chain_set(self):
        """All chain IDs present in the raw dataset."""  # noqa: D401
        return set(self.data_array[:, col_chainID])

    @property
    def fasta(self):
        """
        FASTA sequence of the :attr:`filtered_atoms` lines.

        HETATM residues with non-canonical codes are represented as X.
        """
        c, rs, rn = col_chainID, col_resSeq, col_resName

        chains = defaultdict(dict)
        for row in self.filtered_atoms:
            chains[row[c]].setdefault(row[rs], DEFS.aa3to1.get(row[rn], 'X'))

        return {
            chain: ''.join(residues.values())
            for chain, residues in chains.items()
            }

    #@property
    #def residue_segments(self):
    #    fa = self.filtered_atoms
    #    return np.split(fa, np.where(np.diff(fa[:, col_resSeq]) != 1)[0]+1)

    @property
    def residue_segments(self):
    # this is crude implementation compatible with multiprocessing.Pool
    # the numpy version is commented above
    # if you are reading this, rise an issue and let me know how to
    # handle the creation of new arrays with Pool.map/imap :-)
        return split_residue_segments(self.filtered_atoms)

    def pop_last_filter(self):
        """Pop last filter."""
        self._filters.pop()

    def add_filter(self, function):
        """Add a function as filter."""
        self.filters.append(function)

    def add_filter_record_name(self, record_name):
        """Add filter for record names."""
        self.filters.append(
            lambda x: x[col_record].startswith(record_name)
            )

    def add_filter_chain(self, chain):
        """Add filters for chain."""
        self.filters.append(lambda x: x[col_chainID] == chain)

    def add_filter_backbone(self, minimal=False):
        """Add filter to consider only backbone atoms."""
        ib = is_backbone
        self.filters.append(
            lambda x: ib(x[col_name], x[col_element], minimal=minimal)
            )

    def renumber_atoms(self):
        self.data_array[:, col_serial] = np.arange(1, self.data_array.shape[0] + 1).astype('<U8')

    def get_PDB(self, pdb_filter=None):
        """
        """
        def _(i, f):
            return f(i)

        fs = self.filtered_atoms
        fs[:, col_serial] = np.arange(1, fs.shape[0] + 1).astype('<U8')
        pdb_filter = pdb_filter or []

        lines = reduce(_, pdb_filter, structure_to_pdb(fs))
        return lines

    def write_PDB(self, filename, **kwargs):
        """Write Structure to PDB file."""
        lines = self.get_PDB(**kwargs)
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                write_PDB(lines, filename)
            except UserWarning:
                raise EXCPTS.EmptyFilterError(filename)


def parse_pdb_to_array(datastr, which='both', **kwargs):
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
    numpy.ndarray of (N, len(`libpdb.atom_slicers`))
        Where N are the number of ATOM and/or HETATM lines,
        and axis=1 the number of fields in ATOM/HETATM lines according
        to the PDB format v3.
    """
    # require
    assert isinstance(datastr, str), \
        f'`datastr` is not str: {type(datastr)} instead'

    _ = datastr[datastr.find('MODEL'):datastr.find('ENDMDL')].split('\n')[1:-1]
    if _:
        lines = _
    else:
        lines = datastr.split('\n')

    record_lines = filter_record_lines(lines, which=which)
    data_array = gen_empty_structure_data_array(len(record_lines))
    populate_structure_array_from_pdb(record_lines, data_array)
    return data_array


def parse_cif_to_array(datastr, **kwargs):
    """
    Parse mmCIF protein data to array.

    Array is as given by :func:`gen_empty_structure_data_array`.
    """
    cif = CIFParser(datastr)
    number_of_atoms = len(cif)
    data_array = gen_empty_structure_data_array(number_of_atoms)

    for ii in range(number_of_atoms):
        data_array[ii, :] = cif.get_line_elements_for_PDB(line=ii)

    return data_array


def gen_empty_structure_data_array(number_of_atoms):
    """
    Generate an array data structure to contain structure data.

    Parameters
    ----------
    number_of_atoms : int
        The number of atoms in the structure.
        Determines the size of the axis 0 of the structure array.

    Returns
    -------
    np.ndarray of (N, :attr:`libpdb.atom_slicers), dtype = '<U8'
        Where N is the ``number_of_atoms``.
    """
    # require
    assert isinstance(number_of_atoms, int), \
        f'`number_of_atoms` is not int, {type(number_of_atoms)} '
    assert number_of_atoms > 0, \
        f'or number is less than zero: {number_of_atoms}.'

    return np.empty(
        (number_of_atoms, len(libpdb.atom_slicers)),
        dtype='<U8',
        )


def populate_structure_array_from_pdb(record_lines, data_array):
    """
    Populate structure array from PDB lines.

    Parameters
    ----------
    record_lines : list-like
        The PDB record lines (ATOM or HETATM) to parse.

    data_array : np.ndarray
        The array to populate.

    Returns
    -------
    None
        Populates array in place.
    """
    for row, line in enumerate(record_lines):
        for column, slicer_item in enumerate(libpdb.atom_slicers):
            data_array[row, column] = line[slicer_item].strip()


def filter_record_lines(lines, which='both'):
    """Filter lines to get record lines only."""
    record_headings = record_line_headings
    try:
        # returns lines because needs len after
        return list(
            filter(
                lambda x: x.startswith(record_headings[which]),
                lines,
                ),
            )
    except KeyError as err:
        err2 = ValueError(f'`which` got an unexpected value \'{which}\'.')
        raise err2 from err


def get_datastr(data):
    """
    Get data in string format.

    Can parse data from several formats:

    * Path, reads file content
    * bytes, converst to str
    * str, returns the input

    Returns
    -------
    str
        That represents the data
    """
    t2s = type2string
    data_type = type(data)
    try:
        datastr = t2s[data_type](data)
    except KeyError as err:
        err2 = NotImplementedError('Struture data not of proper type')
        raise err2 from err
    assert isinstance(datastr, str)
    return datastr


def detect_structure_type(datastr):
    """
    Detect Structure data parser.

    Uses `structure_parsers`.

    Returns
    -------
    func or class
        That which can parse `datastr` to a :class:`Structure'.
    """
    sp = structure_parsers
    for condition, parser in sp:
        if condition(datastr):
            return parser
    raise EXCPTS.ParserNotFoundError


def write_PDB(lines, filename):
    """
    Write Structure data format to PDB.

    Parameters
    ----------
    lines : list or np.ndarray
        Lines contains PDB data as according to `parse_pdb_to_array`.

    filename : str or Path
        The name of the output PDB file.
    """
    # use join here because lines can be a generator
    concat_lines = '\n'.join(lines)
    if concat_lines:
        with open(filename, 'w') as fh:
            fh.write(concat_lines)
            fh.write('\n')
        #log.info(S(f'saved: {filename}'))
    else:
        warnings.warn('Empty lines, nothing to write, ignoring.', UserWarning)


def structure_to_pdb(atoms):
    """
    Convert table to PDB formatted lines.

    Parameters
    ----------
    atoms : np.ndarray, shape (N, 16) or similar data structure
        Where N is the number of atoms and 16 the number of cols.

    Yields
    ------
    Formatted PDB line according to `libpdb.atom_line_formatter`.
    """
    for line in atoms:
        values = [func(i) for i, func in zip(line, libpdb.atom_format_funcs)]
        values[col_name] = libpdb.format_atom_name(
            values[col_name],
            values[col_element],
            )
        yield libpdb.atom_line_formatter.format(*values)


col_record = 0
col_serial = 1
col_name = 2
col_altLoc = 3
col_resName = 4
col_chainID = 5
col_resSeq = 6
col_iCode = 7
col_x = 8
col_y = 9
col_z = 10
col_occ = 11
col_temp = 12
col_segid = 13
col_element = 14
col_model = 15


cols_coords = [col_x, col_y, col_z]


# this servers read_pdb_data_to_array mainly
# it is here for performance
record_line_headings = {
    'both': ('ATOM', 'HETATM'),
    'ATOM': 'ATOM',
    'HETATM': 'HETATM',
    }

# possible data inputs in __init__
# all should return a string
# used for control flow
type2string = {
    type(Path()): lambda x: x.read_text(),
    bytes: lambda x: x.decode('utf_8'),
    str: lambda x: x,
    }

# order matters
structure_parsers = [
    (is_cif, parse_cif_to_array),
    (libpdb.is_pdb, parse_pdb_to_array),
    ]


def _APPLY_FILTER(it, func):
    return filter(func, it)



def is_backbone(atom, element, minimal=False):
    """
    Whether `atom` is a protein backbone atom or not.

    Parameters
    ----------
    atom : str
        The atom name.

    element : str
        The element name.

    minimal : bool
        If `True` considers only `C` and `N` elements.
        `False`, considers also `O`.
    """
    e = element.strip()
    a = atom.strip()
    pure_atoms = {
        True: ('N', 'C'),
        False: ('N', 'C', 'O'),
        }
    return e in pure_atoms[minimal] and a in ('N', 'CA', 'O', 'C')


def split_residue_segments(atom_lines):
    """
    Splits a sequence of atoms into consecutive residue segments.

    Parameters
    ----------
    atom_lines : np.ndarray of shape (N, col*)
        A :class:`Structure.data_array` like.
    """
    segments = [[]]
    prev = int(atom_lines[0, col_resSeq])
    for line in atom_lines:
        curr = int(line[col_resSeq])
        if 0 <= curr - prev < 2:
            segments[-1].append(line)
        else:
            segments.append([])
            segments[-1].append(line)
        prev = curr
    return segments


