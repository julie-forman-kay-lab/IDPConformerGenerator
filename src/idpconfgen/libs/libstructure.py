"""
Store internal protein structure representation.

Classes
-------
Structure
    The main API that represents a protein structure in IDPConfGen.
"""
import traceback
import warnings
from collections import defaultdict
from functools import reduce

import numpy as np

from idpconfgen import log
from idpconfgen.core.definitions import aa3to1, blocked_ids, pdb_ligand_codes
from idpconfgen.core.definitions import residue_elements as _allowed_elements
from idpconfgen.core.exceptions import (
    EmptyFilterError,
    NotBuiltError,
    ParserNotFoundError,
    PDBFormatError,
    )
from idpconfgen.libs import libpdb
from idpconfgen.libs.libcif import CIFParser, is_cif
from idpconfgen.libs.libparse import group_runs, sample_case, type2string
from idpconfgen.libs.libpdb import RE_ENDMDL, RE_MODEL, delete_insertions
from idpconfgen.logger import S


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

    def __len__(self):
        return self.data_array.shape[0]

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
            raise NotBuiltError(errmsg=errmsg) from err

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
    def coords(self):
        """
        Coordinates of the filtered atoms.

        As float.
        """
        return self.filtered_atoms[:, cols_coords].astype(np.float32)




    @property
    def consecutive_residues(self):
        """Consecutive residue groups from filtered atoms."""
        # the structure.consecutive_residues
        # this will ignore the iCode
        # please use pdb-tool pdb_delinsertions before this step
        # PDBs downloaded with IDPConfGen already correct for these
        # see libs.libparse.delete_insertions
        return group_runs(self.residues)

    @property
    def fasta(self):
        """
        FASTA sequence of the :attr:`filtered_atoms` lines.

        HETATM residues with non-canonical codes are represented as X.
        """
        c, rs, rn = col_chainID, col_resSeq, col_resName

        chains = defaultdict(dict)
        for row in self.filtered_atoms:
            chains[row[c]].setdefault(row[rs], aa3to1.get(row[rn], 'X'))

        return {
            chain: ''.join(residues.values())
            for chain, residues in chains.items()
            }

    @property
    def filtered_residues(self):
        """Filter residues according to :attr:`filters`."""
        FA = self.filtered_atoms
        return [int(i) for i in dict.fromkeys(FA[:, col_resSeq])]

    @property
    def residues(self):
        """
        Residues of the structure.

        Without filtering, without chain separation.
        """
        return [int(i) for i in dict.fromkeys(self.data_array[:, col_resSeq])]

    @property
    def sorted_minimal_backbone_coords(self, filtered=False):
        """
        Generate a copy of the backbone coords sorted.

        Sorting according N, CA, C.

        This method was created because some PDBs may not have the
        backbone atoms sorted properly.

        Parameters
        ----------
        filtered : bool, optional
            Whether consider current filters or raw data.
        """
        atoms = self.filtered_atoms if filtered else self.data_array

        N_coords = coords[atoms[:, col_name] == 'N']
        CA_coords = coords[atoms[:, col_name] == 'CA']
        C_coords = coords[atoms[:, col_name] == 'C']

        N_num = N_coords.shape[0]
        CA_num = CA_coords.shape[0]
        C_num = C_coords.shape[0]
        num_backbone_atoms = sum([N_num, CA_num, C_num])
        assert num_backbone_atoms / 3 == N_num

        minimal_backbone = np.zeros((num_backbone_atoms, 3), dtype=np.float32)
        minimal_backbone[0:-2:3] = N_coords
        minimal_backbone[1:-1:3] = CA_coords
        minimal_backbone[2::3] = C_coords

        return minimal_backbone

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

    def get_PDB(self, pdb_filters=None, renumber=True):
        """
        Convert Structure to PDB format.

        Considers only filtered lines.

        Returns
        -------
        generator
        """
        def _(i, f):
            return f(i)

        fs = self.filtered_atoms

        # renumber atoms
        if renumber:
            try:
                fs[:, col_serial] = np.arange(1, fs.shape[0] + 1).astype('<U8')
            except IndexError as err:
                errmsg = (
                    'Could not renumber atoms, most likely, because '
                    'there are no lines in selection.'
                    )
                err2 = EmptyFilterError(errmsg)
                raise err2 from err

        pdb_filters = pdb_filters or []

        lines = list(reduce(_, pdb_filters, structure_to_pdb(fs)))

        return lines

    def write_PDB(self, filename, **kwargs):
        """Write Structure to PDB file."""
        lines = self.get_PDB(**kwargs)
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                write_PDB(lines, filename)
            except UserWarning:
                raise EmptyFilterError(filename)


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
    numpy.ndarray of (N, len(`libpdb.atom_slicers`))
        Where N are the number of ATOM and/or HETATM lines,
        and axis=1 the number of fields in ATOM/HETATM lines according
        to the PDB format v3.
    """
    # require
    assert isinstance(datastr, str), \
        f'`datastr` is not str: {type(datastr)} instead'

    model_idx = RE_MODEL.search(datastr)
    endmdl_idx = RE_ENDMDL.search(datastr)

    if bool(model_idx) + bool(endmdl_idx) == 1:
        # only one is True
        raise PDBFormatError('Found MODEL and not ENDMDL, or vice-versa')

    if model_idx and endmdl_idx:
        start = model_idx.span()[1]
        end = endmdl_idx.span()[0] - 1
        lines = datastr[start: end].split('\n')

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
    cif = CIFParser(datastr, **kwargs)
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


def generate_residue_labels(*residue_labels, fmt=None, delimiter=' - '):
    """
    Generate residue labels column.

    Concatenate labels in `residue_labels` using
        `concatenate_residue_labels`.

    Parameters
    ----------
    fmt : str, optional
        The string formatter by default we consider backbone atoms
        of a protein with less than 1000 residues.
        Defaults to `None`, uses '{:<8}', 8 or multiple of 8 according
        to length of residue_labels.
    """
    if not fmt:
        # 11 because 8 + len(' - ')
        fmt = '{:<' + str(len(residue_labels) * (8 + len(delimiter))) + '}'

    concat = (concatenate_residue_labels(l) for l in residue_labels)
    return [fmt.format(delimiter.join(clabels)) for clabels in zip(*concat)]


def concatenate_residue_labels(labels):
    """
    Concatenate residue labels.

    This function is a generator.

    Parameters
    ----------
    labels : numpy array of shape (N, M)
        Where N is the number of rows, and M the number of columns
        with the labels to be concatenated.
    """
    empty_join = ''.join
    return (empty_join(res_label) for res_label in labels)


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
    RH = record_line_headings
    try:
        return [line for line in lines if line.startswith(RH[which])]
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
    raise ParserNotFoundError


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
col_name = 2  # atom name
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
    elements = {
        True: ('N', 'C'),
        False: ('N', 'C', 'O'),
        }
    # elements is needed because of atoms in HETATM entries
    # for example 'CA' is calcium
    return a in ('N', 'CA', 'C', 'O') and e in elements[minimal]


def save_structure_by_chains(
        pdb_data,
        pdbname,
        altlocs=('A', '', ' ', '1'),   # CIF: 6uwi chain D has altloc 1
        chains=None,
        record_name=('ATOM', 'HETATM'),
        renumber=True,
        **kwargs,
        ):
    """
    Save PDBs/mmCIF in separated chains (PDB format).

    Logic to parse PDBs from RCSB.
    """
    # local assignments for speed boost :D
    _DR = pdb_ligand_codes  # discarded residues
    _AE = _allowed_elements
    _S = Structure
    _BI = blocked_ids
    _DI = [delete_insertions]

    pdbdata = _S(pdb_data)
    pdbdata.build()

    chain_set = pdbdata.chain_set

    chains = chains or chain_set

    add_filter = pdbdata.add_filter
    pdbdata.add_filter_record_name(record_name)
    add_filter(lambda x: x[col_resName] not in _DR)
    add_filter(lambda x: x[col_element] in _AE)
    add_filter(lambda x: x[col_altLoc] in altlocs)

    for chain in chains:

        # writes chains always in upper case because chain IDs given by
        # Dunbrack lab are always in upper case letters
        chaincode = f'{pdbname}_{chain}'

        # this operation can't be performed before because
        # until here there is not way to assure if the chain being
        # downloaded is actualy in the blocked_ids.
        # because at the CLI level the user can add only the PDBID
        # to indicate download all chains, while some may be restricted
        if chaincode in _BI:
            log.info(S(
                f'Ignored code {chaincode} because '
                'is listed in blocked ids.'
                ))
            continue

        # upper and lower case combinations:
        possible_cases = sample_case(chain)
        # cases that exist in the structure
        cases_that_actually_exist = chain_set.intersection(possible_cases)
        # this trick places `chain` first in the for loop because
        # it has the highest probability to be the one required
        cases_that_actually_exist.discard(chain)
        probe_cases = [chain] + list(cases_that_actually_exist)

        for chain_case in probe_cases:

            pdbdata.add_filter_chain(chain_case)
            fout = f'{chaincode}.pdb'

            try:
                pdb_lines = pdbdata.get_PDB(pdb_filters=_DI)
            except EmptyFilterError as err:
                err2 = \
                    EmptyFilterError(f'for chain {pdbname}_{chain_case}')
                errlog = (
                    f'{repr(err)}\n'
                    f'{repr(err2)}\n'
                    f'{traceback.format_exc()}\n'
                    'continuing to new chain\n'
                    )
                log.debug(errlog)
                continue
            else:
                if all(line.startswith('HETATM') for line in pdb_lines):
                    log.debug(
                        f'Found only HETATM for {chain_case}, '
                        'continuing with next chain.'
                        )
                    continue
                yield fout, '\n'.join(pdb_lines)
                break
            finally:
                pdbdata.pop_last_filter()
        else:
            log.debug(f'Failed to download {chaincode}')
