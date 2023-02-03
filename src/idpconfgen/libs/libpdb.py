"""Contain  handlers of PDB information."""
import functools
import re
from collections import defaultdict

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.core.definitions import aa3to1
from idpconfgen.logger import S


RE_MODEL = re.compile('(\n|^)MODEL\s*\d+\s*\n')  # noqa: W605
RE_ENDMDL = re.compile('\nENDMDL\s*\n')  # noqa: W605


# string formats for atom name
_3 = ' {:<3s}'
_4 = '{:<4s}'
# len of element, atom formatting string
# anything outside this is an error
_atom_format_dict = {
    # len of element
    1: {
        # len of atom name
        1: _3,
        2: _3,
        3: _3,
        4: _4,
        },
    2: {
        1: _4,
        2: _4,
        3: _4,
        4: _4,
        },
    }


def delete_insertions(lines):
    """
    Delete insertions.

    Adapted from pdbtools and optimized for this context.

    Visit pdb-tools at:
    https://github.com/haddocking/pdb-tools/blob/master/pdbtools/pdb_delinsertion.py
    """  # noqa: E501
    # Keep track of residue numbering
    # Keep track of residues read (chain, resname, resid)
    offset = 0
    prev_resi = None
    seen_ids = set()
    clean_icode = False
    for line in lines:
        res_uid = line[17:27]  # resname, chain, resid
        id_res = line[21] + line[22:26].strip()  # A99, B12
        icode = line[26]

        # unfortunately, this is messy but not all PDB files follow a nice
        # order of ' ', 'A', 'B', ... when it comes to insertion codes..
        if prev_resi != res_uid:  # new residue
            # Have we seen this chain + resid combination
            # catch insertions WITHOUT icode ('A' ... ' ' ... 'B')
            if id_res in seen_ids:
                # Should we do something about it?
                clean_icode = True
                offset += 1
            # Do we have an explicit icode?
            elif icode != ' ':
                if id_res in seen_ids:  # never saw this, do not offset!
                    offset += 1
                clean_icode = True
            else:
                clean_icode = False

            prev_resi = res_uid

        if clean_icode:
            line = f'{line[:26]} {line[27:]}'

        resid = int(line[22:26]) + offset
        line = f'{line[:22]}{str(resid).rjust(4)}{line[26:]}'
        seen_ids.add(id_res)
        yield line


def format_atom_name(atom, element, AFD=_atom_format_dict):
    """
    Format PDB Record line Atom name.

    Further Reading:

    * https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html

    Parameters
    ----------
    atom : str
        The atom name.

    element : str
        The atom element code.

    Returns
    -------
    str
        Formatted atom name.
    """
    atm = atom.strip()
    len_atm = len(atm)
    len_ele = len(element.strip())

    try:
        return AFD[len_ele][len_atm].format(atm)
    except KeyError as err:
        _ = f'Could not format this atom:type -> {atom}:{element}'
        # raising KeyError assures that no context in IDPConfGen
        # will handle it. @joaomcteixeira never handles pure Python
        # exceptions, those are treated as bugs.
        raise KeyError(_) from err


def format_chainid(chain):
    """
    Format chain identifier to one letter.

    This is required to receive chain IDs from mmCIF files, which
    may have more than one letter.
    """
    return chain.strip()[0]


atom_record = slice(0, 6)
atom_serial = slice(6, 11)
atom_name = slice(12, 16)
atom_altLoc = slice(16, 17)
atom_resName = slice(17, 20)
atom_chainID = slice(21, 22)
atom_resSeq = slice(22, 26)
atom_iCode = slice(26, 27)
atom_x = slice(30, 38)
atom_y = slice(38, 46)
atom_z = slice(46, 54)
atom_occ = slice(54, 60)
atom_temp = slice(60, 66)
atom_segid = slice(72, 76)
atom_element = slice(76, 78)
atom_model = slice(78, 80)

# order matters
atom_slicers = [
    atom_record,
    atom_serial,
    atom_name,
    atom_altLoc,
    atom_resName,
    atom_chainID,
    atom_resSeq,
    atom_iCode,
    atom_x,
    atom_y,
    atom_z,
    atom_occ,
    atom_temp,
    atom_segid,
    atom_element,
    atom_model,
    ]


atom_line_formatter = (
    "{:6s}"
    "{:5d} "
    "{}"
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
    "{:>2s}"
    "{:2s}"
    )


# functions to apply to each format field in atom line string
def _nothing(x):
    return x


atom_format_funcs = [
    _nothing, int, _nothing, _nothing, _nothing,
    format_chainid, int, _nothing, float, float,
    float, float, float, _nothing, _nothing,
    _nothing,
    ]


# USED OKAY
def get_fasta_from_PDB(pdbid):
    """Extract FASTA from PDB."""
    try:
        lines = pdbid[1].decode('utf_8').split('\n')
    except AttributeError:
        lines = pdbid[1].split('\n')
    rn = {line[atom_resSeq].strip(): line[atom_resName] for line in lines}
    fasta = (aa3to1.get(f, 'X') for f in rn.values())
    return Path(pdbid[0]).stem, ''.join(fasta)


def is_pdb(datastr):
    """Detect if `datastr` if a PDB format v3 file."""
    assert isinstance(datastr, str), \
        f'`datastr` is not str: {type(datastr)} instead'
    return bool(datastr.count('\nATOM ') > 0)


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
    rgx_XXXX_C_segS = re.compile(r'^[0-9a-zA-Z]{4}_[0-9a-zA-Z]+_seg\d+(\s|$)')
    rgx_any = re.compile(r'\w+(\.pdb)?')

    def __new__(cls, name):
        """Construct class."""
        if isinstance(name, PDBID):
            return name

        namep = Path(name)
        if namep.suffix in ('.pdb', '.cif'):
            name = namep.stem

        # where XXXX is the PDBID and C the chain ID
        pdb_filename_regex = {
            cls.rgx_XXXX: cls._parse_XXXX,
            cls.rgx_XXXXC: cls._parse_XXXXC,
            cls.rgx_XXXX_C: cls._parse_XXXX_C,
            cls.rgx_XXXX_C_segS: cls._parse_XXXX_C_segS,
            cls.rgx_any: cls._parse_any_name,
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
        pdbid, chainid, *_ = pdbid.split()[0].split('_')
        return pdbid, chainid

    @staticmethod
    def _parse_XXXX_C_segS(pdbid):
        pdbid, chainid, segID = pdbid.split()[0].split('_')
        return pdbid, chainid, segID.lstrip('seg')

    @staticmethod
    def _parse_any_name(pdbid):
        return [Path(pdbid).stem]


# USED OKAY
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

        if isinstance(pdb_names, cls):
            return pdb_names
        else:
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
        return '{} with {} element(s)'.format(
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

    @property
    def pdbids(self):
        """Generate the PDBID names."""
        return (pdbid.name for pdbid in self)

    @property
    def name_chains_dict(self):
        """Export PDBIDs: Chains dictionary map."""
        name_chains = defaultdict(list)
        for pdbid in self:
            pc = pdbid.chain
            if pc:
                name_chains[pdbid.name].append(pc)
            else:
                name_chains[pdbid.name]
        return name_chains

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


# USED OKAY
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

    # segment are fragments of the same chain
    def __init__(self, name, chain=None, segment=None):

        self.name = name.upper()
        self.chain = chain

        # made manual to completely control order
        ids = {
            'chain': chain,
            'seg': segment,
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
