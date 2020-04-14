"""Handle CIF data."""
import re

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS

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


# thanks to @JoaoRodrigues (GitHub) for this cif regex
CIF_LINE_REGEX = re.compile(
    r'''
    '(.*?)' | # single quoted substring OR
    "(.*?)" | # double quoted substring OR
    (\S+)     # any consec. non-whitespace character(s)
    ''',
    re.VERBOSE,
    )

_void_translation_table = str.maketrans('?.', '  ')


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
    CIFFileInvalidError
        If any `_atom_site.` keys are found in ``lines``.
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
    raise EXCPTS.CIFFileInvalidError(errmsg)


def populate_cif_dictionary(lines, start_index, cif_dict):
    """

    Parameters
    ----------

    Returns
    -------

    Raises
    ------
    """
    CLR = CIF_LINE_REGEX
    counter = 0
    for line in lines[start_index:]:
        if line.startswith('#'):
            return counter

        ls = [''.join(t) for t in CLR.findall(line)]

        for i, key in enumerate(cif_dict.keys()):
            try:
                cif_dict[key].append(ls[i])
            except IndexError as err:
                errmsg = f'CIF line did not split properly: {ls}'
                EXCPTS.CIFFileInvalidError(errmsg)

        counter += 1

    errmsg = 'Could not find the \'#\' to end CIF reading.'
    raise EXCPTS.CIFFileInvalidError(errmsg)


class CIFParser:
    """
    """
    def __init__(self, datastr):
        """
        """
        self.cif_dict = {}
        self.number_of_atoms = None
        self.read_cif(datastr)

    def read_cif(self, datastr):
        """Read 'atom_site' entries to dictionary."""
        lines = datastr.split('\n')
        atom_start_index = find_cif_atom_site_headers(lines, self.cif_dict)
        self._len = populate_cif_dictionary(
            lines,
            atom_start_index,
            self.cif_dict,
            )

    def __len__(self):
        return self._len

    @property
    def line(self):
        return self._line

    @line.setter
    def line(self, ii):
        assert isinstance(ii, int)
        assert 0 <= ii <= len(self)
        self._line = ii

    def get_line_elements_for_PDB(self, i=None):
        """
        """
        # http://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html
        if i:
            self.line = i
        return [
            self.record,
            self.serial,
            self.atname,
            self.altloc,
            self.resname,
            self.chainid,
            self.resseq,
            self.icode,
            self.x,
            self.y,
            self.z,
            self.occ,
            self.tempfactor,
            '    ',  # segid
            self.element,
            self.charge,
            ]

    def get_value(self, label):
        return self.cif_dict[label][self.line]

    def _auth_label(self, label):
        try:
            return self.get_value(f'_atom_site.label_{label}')
        except KeyError:
            return self.get_value(f'_atom_site.auth_{label}')

    def _translate(self, item):
        return item.translate(_void_translation_table)

    @property
    def record(self):
        return self.get_value('_atom_site.group_PDB')

    @property
    def serial(self):
        try:
            return self.get_value('_atom_site.Id')
        except KeyError:
            return self.get_value('_atom_site.id')

    @property
    def atname(self):
        return self._auth_label('atom_id')

    @property
    def altloc(self):
        altloc = self._auth_label('alt_id')
        return self._translate(altloc)

    @property
    def resname(self):
        return self._auth_label('comp_id')

    @property
    def chainid(self):
        value = self._auth_label('asym_id')
        if value in ('?', '.'):
            value = self.get_value(f'_atom_site.auth_asym_id')
        return value

    @property
    def resseq(self):
        return self._auth_label('seq_id')

    @property
    def icode(self):
        try:
            icode = self.get_value('_atom_site.pdbx_PDB_ins_code')
            return self._translate(icode)
        except KeyError:
            return " "

    @property
    def xcoord(self):
        return self.get_value('_atom_site.Cartn_x')

    @property
    def ycoord(self):
        return self.get_value('_atom_site.Cartn_y')

    @property
    def zcoord(self):
        return self.get_value('_atom_site.Cartn_z')

    @property
    def tempfactor(self):
        return self.get_value('_atom_site.B_iso_or_equiv')

    @property
    def element(self):
        return self.get_value('_atom_site.type_symbol')

    @property
    def charge(self):
        charge = self.get_value('_atom_site.pdbx_formal_charge')
        return self._translate(charge)

