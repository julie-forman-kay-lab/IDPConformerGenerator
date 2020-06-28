"""Handle CIF data."""
import re

from idpconfgen.core import exceptions as EXCPTS


class CIFParser:
    """
    mmCIFParser for structural data ONLY.

    That is, fields from the m_site.` group.

    Parameters
    ----------
    datastr : str
        The contect of the mmCIF file in string format.
    """

    __slots__ = [
        '_auth_label',
        '_len',
        '_line',
        'cif_dict',
        ]

    def __init__(self, datastr, label_priority='auth'):

        _label_priority = {
            'auth': self._auth_label_priority,
            'label': self._label_auth_priority,
            }

        self._auth_label = _label_priority[label_priority]
        self.cif_dict = {}
        self.read_cif(datastr)
        self.line = 0

    def read_cif(self, datastr):
        """Read 'atom_site.' entries to dictionary."""
        lines = datastr.split('\n')
        atom_start_index = find_cif_atom_site_headers(lines, self.cif_dict)
        self._len = populate_cif_dictionary(
            lines,
            atom_start_index,
            self.cif_dict,
            )

    def __len__(self):
        """The length is defined by the number of atoms."""  # noqa: D401
        return self._len

    @property
    def line(self):
        """The current line number."""  # noqa: D401
        return self._line

    @line.setter
    def line(self, ii):
        """Set the current line number."""
        assert isinstance(ii, int)
        assert 0 <= ii <= len(self)
        self._line = ii

    def get_line_elements_for_PDB(self, line=None):
        """
        Retrieve `line` elements.

        Parameters
        ----------
        line : int, optional
            If given retrieve values for that line and sets that as
            current line. Else, retrieves values for current :attr:`line`.

        Return
        ------
        list
            The values of the line.
        """
        # http://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html
        if line:
            self.line = line
        return [
            self.record,
            self.serial,
            self.atname,
            self.altloc,
            self.resname,
            self.chainid,
            self.resseq,
            self.icode,
            self.xcoord,
            self.ycoord,
            self.zcoord,
            self.occ,
            self.tempfactor,
            ' ',  # segid
            self.element,
            self.charge,
            ]

    def get_value(self, label):
        """Retrieve the label value for the current :attr:`line`."""
        return self.cif_dict[label][self.line]

    def _auth_label_priority(self, label):
        try:
            return self.get_value(f'_atom_site.auth_{label}')
        except KeyError:
            return self.get_value(f'_atom_site.label_{label}')

    def _label_auth_priority(self, label):
        try:
            return self.get_value(f'_atom_site.label_{label}')
        except KeyError:
            return self.get_value(f'_atom_site.auth_{label}')

    def _translate(self, item):
        return item.translate(_void_translation_table)

    @property
    def record(self):
        """The record field of the current :attr:`line`."""  # noqa: D401
        return self.get_value('_atom_site.group_PDB')

    @property
    def serial(self):
        """The serial field of the current :attr:`line`."""  # noqa: D401
        try:
            return self.get_value('_atom_site.Id')
        except KeyError:
            return self.get_value('_atom_site.id')

    @property
    def atname(self):
        """The atom name field of the current :attr:`line`."""  # noqa: D401
        return self._auth_label('atom_id')

    @property
    def altloc(self):
        """The altloc field of the current :attr:`line`."""  # noqa: D401
        altloc = self._auth_label('alt_id')
        return self._translate(altloc)

    @property
    def resname(self):
        """The resname field of the current :attr:`line`."""  # noqa: D401
        return self._auth_label('comp_id')

    @property
    def chainid(self):
        """The chainID field of the current :attr:`line`."""  # noqa: D401
        value = self._auth_label('asym_id')
        if value in ('?', '.'):
            value = self.get_value('_atom_site.auth_asym_id')
        return value

    @property
    def resseq(self):
        """The resSeq field of the current :attr:`line`."""  # noqa: D401
        return self._translate(self._auth_label('seq_id'))

    @property
    def icode(self):
        """The icode field of the current :attr:`line`."""  # noqa: D401
        try:
            icode = self.get_value('_atom_site.pdbx_PDB_ins_code')
            return self._translate(icode)
        except KeyError:
            return " "

    @property
    def xcoord(self):
        """The Xcoord field of the current :attr:`line`."""  # noqa: D401
        return self.get_value('_atom_site.Cartn_x')

    @property
    def ycoord(self):
        """The Ycoord field of the current :attr:`line`."""  # noqa: D401
        return self.get_value('_atom_site.Cartn_y')

    @property
    def zcoord(self):
        """The Zcoord field of the current :attr:`line`."""  # noqa: D401
        return self.get_value('_atom_site.Cartn_z')

    @property
    def occ(self):
        """The occupancy field of the current :attr:`line`."""  # noqa: D401
        return self.get_value('_atom_site.occupancy')

    @property
    def tempfactor(self):
        """The tempfactor field of the current :attr:`line`."""  # noqa: D401
        return self.get_value('_atom_site.B_iso_or_equiv')

    @property
    def element(self):
        """The element field of the current :attr:`line`."""  # noqa: D401
        return self.get_value('_atom_site.type_symbol')

    @property
    def charge(self):
        """The charge field of the current :attr:`line`."""  # noqa: D401
        charge = self.get_value('_atom_site.pdbx_formal_charge')
        return self._translate(charge)


def is_cif(datastr):
    """Detect if `datastr` is a CIF file."""
    assert isinstance(datastr, str), \
        f'`datastr` is not str: {type(datastr)} instead'
    cif_loop = re.compile('[lL][oO][oO][pP]_')
    return bool(cif_loop.search(datastr))


def find_cif_atom_site_headers(lines, cif_dict):
    """
    Find `_atom_site.` headers.

    Given a list of mmCIF lines, finds headers of `_atom_site.`

    Parameters
    ----------
    lines : list
        The lines of the mmCIF files are read by `file.readlines()`.

    cif_dict : dict
        A dictionary where to initiate the `atom_site.` keys.
        Keys are assigned empty list as value.

    Returns
    -------
    int
        The index of ``lines`` where the `_atom_site.` structural
        information starts.

    Raises
    ------
    CIFFileError
        If any `_atom_site.` keys are found in ``lines``.

    See Also
    --------
    populate_cif_dictionary
    """
    # require
    assert isinstance(lines, list)
    assert isinstance(cif_dict, dict)

    found = False
    for ii, line in enumerate(lines):
        if line.startswith('_atom_site.'):
            found = True
            cif_dict.setdefault(line.strip(), [])
        elif found:
            return ii
    errmsg = 'Could not find `_atom_site.` entries for CIF file'
    raise EXCPTS.CIFFileError(errmsg)


def populate_cif_dictionary(lines, start_index, cif_dict):
    """
    Populate mmCIF dictionary.

    Parameters
    ----------
    lines : list
        The mmCIF lines. Normally the whole mmCIF file content.

    start_index : int
        The line index, that is, line number 0-indexed, where the
        actual `atom_site.` structural information starts.

    cif_dict : dict
        The mmCIF dictionary to populate. The dict keys should be
        previously prepared and should match the data in the `lines`.

    Returns
    -------
    None
        Edits dictionary in place.

    Raises
    ------
    InvalidCIFLineError
        If the number of fields in parsed line do not match expected
        `_atom_site.` fields in the dictionary.

    See Also
    --------
    find_cif_atom_site_headers
    """
    assert len(lines) > start_index

    valid_len = len(cif_dict)
    for counter, line in enumerate(lines[start_index:]):
        if line.startswith('#'):
            return counter

        ls = parse_cif_line(line)

        if len(ls) != valid_len:
            errmsg = (
                "Fields in line do not match number of excepted `_atom_site.*` "
                "fields.\n"
                f"Line read: {line!r}\n"
                f"At: start index {start_index}, counter {counter}."
                )
            raise EXCPTS.CIFFileError(errmsg=errmsg)

        for i, key in enumerate(cif_dict.keys()):
            cif_dict[key].append(ls[i])


def parse_cif_line(line):
    """
    Parse CIF line according to `CIF_LINE_REGEX`.

    Parameters
    ----------
    line : str
        The line to parse.

    Returns
    -------
    list
        A list with the parsed line elements.
    """
    CLR = cif_line_regex
    return [''.join(t).strip() for t in CLR.findall(line)]


cif_atom_keys = [
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

# thanks to @JoaoRodrigues (GitHub) for this cif regex
cif_line_regex = re.compile(
    r'''
    '(.*?)' | # single quoted substring OR
    "(.*?)" | # double quoted substring OR
    (\S+)     # any consec. non-whitespace character(s)
    ''',
    re.VERBOSE,
    )

_void_translation_table = str.maketrans('?.', '  ')
