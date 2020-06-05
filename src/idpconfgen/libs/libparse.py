"""Contain different parsing strategies for different files."""
import string
import subprocess
import sys
import traceback
from contextlib import contextmanager
from os import fspath

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.core.definitions import blocked_ids
from idpconfgen.core.definitions import dssp_trans, pdb_ligand_codes
from idpconfgen.libs import libcheck
from idpconfgen.libs.libpdb import PDBIDFactory
from idpconfgen.libs.libstructure import (
    Structure,
    col_altLoc,
    col_element,
    col_iCode,
    col_name,
    col_resName,
    col_resSeq,
    type2string,
    )
from idpconfgen.logger import S


_ascii_lower_set = set(string.ascii_lowercase)
_ascii_upper_set = set(string.ascii_uppercase)
_minimal_bb_atoms = ['N', 'CA', 'C']  # ordered!


def delete_insertions(lines):
    """
    Delete insertions.

    Adapted from pdbtools and optimized for this context.
    """
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


def parse_dssp(data, reduced=False):
    """Parse DSSP file data."""
    DT = dssp_trans

    try:
        data_ = type2string[type(data)](data).split('\n')
    except KeyError:  # maybe is list
        data_ = data

    # RM means removed empty
    RM1 = (i for i in data_ if i)

    # exausts generator until
    for line in RM1:
        if line.strip().startswith('#'):
            break
    else:
        # if the whole generator is exhausted
        raise IndexError

    structure_data = (i for i in RM1 if i[13] != '!')

    dssp = []
    fasta = []
    residues = []
    for line in structure_data:
        dssp.append(line[16])
        fasta.append(line[13])
        residues.append(line[6:10].strip())

    if reduced:
        dssp = list(map(lambda x: x.translate(DT), dssp))

    assert sum((len(dssp), len(fasta), len(residues))) / 3 == len(dssp)

    return ''.join(dssp), ''.join(fasta), ','.join(residues)


def find_dssp_data_index(data):
    """
    Find index for where item startswith '#'.

    Evaluates after left striping white spaces.
    """
    for line in data:
        if line.strip().startswith('#'):
            return
    else:
        raise IndexError


def list_index_to_array(
        list_,
        start=None,
        stop=None,
        step=None,
        sObj=None,
        index=None,
        ):
    """
    Extract slices of strings in lists to an array.

    At least one the named parameters must be provided.

    In case more than one parameter set is provided:
        the `start`, `stop` or `step` trio have the highest priority,
        followed by `sObj` and last `index`.

    Raises
    ------
    ValueError
        If no named parameter is provided.

    ValueError
        If sObj is provided and is not a ``slice`` object.
    """
    if any((i is not None for i in (start, stop, step))):
        index = slice(start, stop, step)

    elif sObj is not None:
        if isinstance(sObj, slice):
            index = sObj
        else:
            raise ValueError(
                'sObj must be slice type: {} given.'.format(type(sObj))
                )

    elif index is not None:
        index = index

    else:
        raise ValueError('At least one of the index args most be provided.')

    len_of_slice = len(list_[0][index])

    array = np.empty(
        len(list_),
        dtype='<U{}'.format(len_of_slice),
        )

    for i, line in enumerate(list_):
        array[i] = line[index]

    return array


@libcheck.kwargstype((str, Path, type(None)))
def export_ss_from_DSSP(*dssp_task_results, output='dssp.database'):
    """
    Export secondary structure information from :class:`DSSPParser` results.

    Parameters
    ----------
    dssp_task_results : multiple :class:`DSSPParser`
    """
    output_data = _concatenate_ss_from_dsspparsers(dssp_task_results)

    output_data = '\n'.join(output_data) + '\n'

    if output:
        opath = Path(output)
        opath.myparents().mkdir(parents=True, exist_ok=True)
        opath.write_text(output_data)
    else:
        sys.stdout.write(output_data)


def _concatenate_ss_from_dsspparsers(dsspparsers):
    output = []
    for dsspparser in dsspparsers:
        output.append(
            '{}|{}'.format(
                dsspparser.pdbid,
                ''.join(dsspparser.ss),
                )
            )

    output.sort(key=lambda x: x.split('|')[0])

    return output


# DEPECRATED
def read_pipe_file(text):
    """
    Read IDPConfGen pipe file.

    Pipe files have the following structure:

    CODE|DATA
    CODE|DATA
    (...)

    Parameters
    ----------
    text : str
        The crude text to read.

    Returns
    -------
    dict {str: str}
        Of codes and data.
    """
    lines = text.split('\n')
    elines = filter(bool, lines)
    return {c: s for line in elines for c, s in (line.split('|'),)}


def group_by(data):
    """
    Group data by indexes.

    Parameters
    ----------
    data : iterable
        The data to group by.

    Returns
    -------
    list : [[type, slice],]

    Examples
    --------
    >>> group_by('LLLLLSSSSSSEEEEE')
    [['L', slice(0, 5)], ['S', slice(5, 11)], ['E', slice(11,16)]]
    """
    assert len(data) > 0

    prev = data[0]
    start = 0
    current = [prev]
    groups = []
    for i, datum in enumerate(data):
        if datum != prev:
            current.append(slice(start, i))
            groups.append(current)
            start = i
            current = [datum]
        prev = datum
    else:
        current.append(slice(start, i + 1))
        groups.append(current)

    assert isinstance(groups[0][0], str)
    assert isinstance(groups[0][1], slice)
    return groups


# DEPRECATED
def group_consecutive_ints(sequence):
    """Groupd consecutive ints."""
    prev = sequence[0]
    start = 0
    slices = []
    for i, integer in enumerate(sequence[1:], start=1):
        if integer - prev not in (0, 1):
            slices.append(slice(start, i))
            start = i
        prev = integer
    else:
        slices.append(slice(start, i + 1))
    return slices


# from https://stackoverflow.com/questions/21142231
def group_runs(li, tolerance=1):
    """Group consecutive numbers given a tolerance."""
    out = []
    last = li[0]
    for x in li:
        if x - last > tolerance:
            yield out
            out = []
        out.append(x)
        last = x
    yield out


_allowed_elements = ('C', 'O', 'N', 'H', 'S', 'Se', 'D')


@contextmanager
def try_to_write(data, fout):
    """Context to download."""
    try:
        yield
    except Exception as err:
        log.debug(traceback.format_exc())
        log.error(S('error found for {}: {}', fout, repr(err)))
        erp = Path(fout.myparents(), 'errors')
        erp.mkdir(parents=True, exist_ok=True)
        p = Path(erp, fout.stem).with_suffix('.structure')
        if p.exists():
            return
        else:
            try:
                p.write_bytes(data)
            except TypeError:
                p.write_text(data)


def filter_structure(pdb_path, **kwargs):
    """
    Download a PDB/CIF structure chains.

    Parameters
    ----------
    pdbid : tuple of 2-elements
        0-indexed, the structure ID at RCSB.org;
        1-indexed, a list of the chains to download.

    **kwargs : as for :func:`save_structure_chains_and_segments`.
    """
    pdbid = PDBIDFactory(pdb_path)
    pdbname = pdbid.name
    chains = pdbid.chain

    save_structure_chains_and_segments(
        pdb_path,
        pdbname,
        chains=chains,
        **kwargs,
        )


def eval_chain_case(chain, chain_set):
    """
    Evalute chain case.

    Tests whether chain exists in chain_set in both UPPER and LOWER cases.

    Returns
    -------
    str
        The found chain.
    """
    if chain in chain_set:
        return chain
    else:
        cl = chain.lower()
        if cl in chain_set:
            return cl
    raise ValueError(f'Neither {chain}/{cl} are in set {chain_set}')


def save_structure_chains_and_segments(
        pdb_data,
        pdbname,
        chains=None,
        record_name=('ATOM', 'HETATM'),
        altlocs=('A', '', ' '),
        renumber=True,
        mdict=None,
        ):
    """
    Prase structure.

    Logic to parse PDBs from RCSB.
    """
    _DR = pdb_ligand_codes  # discarded residues
    _AE = _allowed_elements

    pdbdata = Structure(pdb_data)
    pdbdata.build()

    chain_set = pdbdata.chain_set

    chains = chains or chain_set

    pdbdata.add_filter_record_name(record_name)
    pdbdata.add_filter(lambda x: x[col_resName] not in _DR)
    pdbdata.add_filter(lambda x: x[col_element] in _AE)
    pdbdata.add_filter(lambda x: x[col_altLoc] in altlocs)

    for chain in chains:

        # writes chains always in upper case because chain IDs given by
        # Dunbrack lab are always in upper case letters
        # eval_chain_case evaluates for the need for lower case,
        # however if the lower case is kept in the final file
        # it may create incompatibilities
        # 03/Jun/2020
        chaincode = f'{pdbname}_{chain}'

        # this operation can't be performed before because
        # until here there is not way to assure if the chain being
        # downloaded is actualy in the blocked_ids.
        # because at the CLI level the user can add only the PDBID
        # to indicate download all chains, while some may be restricted
        if chaincode in blocked_ids:
            log.info(S(
                f'Ignored code {chaincode} because '
                'is listed in blocked ids.'
                ))
            continue

        fout = f'{chaincode}.pdb'

        try:
            chain = eval_chain_case(chain, chain_set)
        except ValueError as err:
            log.error(repr(err))
            log.error(f'Skiping chain {chain} for {pdbname}')
            continue

        pdbdata.add_filter_chain(chain)

        mdict[fout] = list(pdbdata.get_PDB(pdb_filter=[delete_insertions]))

        pdbdata.pop_last_filter()

    return


def get_segments_based_on_backbone_continuity(atoms):
    """
    Split backbones in chunks of continuity.

    Have to explain big thing here.
    """
    ref = _minimal_bb_atoms
    start = 0
    idx = 0
    slices = []
    max_size = atoms.shape[0]

    bb = atoms[:, col_name]
    resis = np.core.defchararray.add(atoms[:, col_resSeq], atoms[:, col_iCode])

    while idx < max_size:
        a = list(bb[idx: idx + 3])
        bb_continuity_lost = a != ref or len(set(resis[idx: idx + 3])) > 1
        if bb_continuity_lost:
            slices.append(slice(start, idx, None))
            idx += 1
            start = idx
            while idx <= max_size:
                a = list(bb[idx: idx + 3])
                bb_continuity_restored = \
                    a == ref and len(set(resis[idx: idx + 3])) == 1
                if bb_continuity_restored:
                    start = idx
                    idx += 3
                    break  # back to the main while
                idx += 1
        else:
            idx += 3
    else:
        slices.append(slice(start, idx + 1, None))

    return [list(dict.fromkeys(atoms[seg, col_resSeq])) for seg in slices]


def identify_backbone_gaps(atoms):
    """Atoms is expected already only minimal backbone."""
    # this is a tricky implementation, PRs welcomed :-)
    if not set(atoms[:, col_name]).issubset({'N', 'CA', 'C'}):
        raise EXCPTS.PDBFormatError(errmsg='Back bone is not subset')

    resSeq, slices = zip(*group_by(atoms[:, col_resSeq]))
    print(slices)
    if not all(isinstance(i, slice) for i in slices):
        raise TypeError('Expected slices found something else')

    # calculates the length of each slice
    # each slice corresponds to a residue
    slc_len = np.array([slc.stop - slc.start for slc in slices])
    print(slc_len)

    # identifies indexes where len slice differ from three.
    # considering that atoms contains only minimal backbone atoms
    # (N, CA, C) this would mean that some atoms are missing
    aw = np.argwhere(slc_len != 3)

    # splits the slices according to the indexes where != 3
    s = np.split(slices, aw[:, 0])

    # np.split keeps the index where it was split, contrarily to str.split
    # here I filter out those unwanted indexes
    # parsed_split should contain a list of list with only the slices
    # for residues with 3 backbone atoms
    # the first index is a special case because of np.split
    parsed_split = \
        [list(i) for i in s[0:1] if i[1:].size > 0] \
        + [list(i[1:]) for i in s[1:] if i[1:].size > 0]

    # concatenate residue slices in segments
    segment_slices = [slice(i[0].start, i[-1].stop) for i in parsed_split]

    # returns a list of array views
    # return [atoms[seg, :] for seg in segment_slices]
    return [
        list(dict.fromkeys(atoms[seg, col_resSeq]))
        for seg in segment_slices
        ]


def mkdssp(pdb, ss_cmd, mdict=None, reduced=False):
    """
    Execute `mkdssp` from DSSP.

    https://github.com/cmbi/dssp
    """
    cmd = [ss_cmd, '-i', fspath(pdb.resolve())]
    result = subprocess.run(cmd, capture_output=True)

    dssp, fasta, residues = parse_dssp(
        result.stdout.decode('utf-8'),
        reduced=reduced,
        )
    mdict[str(PDBIDFactory(pdb))] = {
        'dssp': dssp,
        'fasta': fasta,
        'resids': residues,
        }

    return


def get_segsplitter(destination):
    """
    Select the appropriate function to write the segments
    """
    return



# DEPECRATED
# class DSSPParser:
#    """
#    Provides an interface for `DSSP files`_.
#
#    .. _DSSP files: https://github.com/cmbi/xssp
#
#    If neither `fin` nor `data` parameters are given, initiates
#    an data empty parser.
#
#    Parameters
#    ----------
#    fin : str or Path object, optional
#        The path to the dssp text file to read.
#        Defaults to ``None``.
#
#    data : str or list, optional
#        A string containing the dssp file data, or a list of lines
#        of such file.
#        Defaults to ``None``.
#
#    pdbid : any, optional
#        An identification for the DSSP file being parsed.
#        Deafults to None.
#
#    reduced : bool
#        Whether or not to reduce secondary structure representation
#        only three types: loops 'L', helices 'H', and sheets 'E'.
#
#    Attributes
#    ----------
#    ss : array
#        If `data` or `fin` are given :attr:`ss` stores the secondary
#        structure information DSSP Keys of the protein.
#
#    """
#
#    def __init__(
#            self,
#            *,
#            fin=None,
#            data=None,
#            pdbid=None,
#            reduced=False,
#            ):
#
#        self.pdbid = pdbid
#        self.reduced = reduced
#
#        if fin:
#            self.read_dssp_data(Path(fin).read_text())
#        elif data:
#            self.read_dssp_data(data)
#        else:
#            self.data = None
#
#
#    def __eq__(self, other):
#        is_equal = [
#            self.data == other.data,
#            self.pdbid == other.pdbid,
#            ]
#
#        return all(is_equal)
#
#    def read_dssp_data(self, data):
#        """
#        Read DSSP data into the parser object.
#
#        Parameters
#        ----------
#        data : str or list of strings
#            Has the same value as Class parameter `data`.
#        """
#        try:
#            data = data.split('\n')
#        except AttributeError:  # data already a list
#            pass
#
#        data = [i for i in data if i]  # removes empty strings
#
#        try:
#            data_header_index = self._finds_data_index(data)
#        except IndexError as err:
#            raise EXCPTS.DSSPParserError(self.pdbid) from err
#
#        # data starts afterthe header
#        #
#        # '!' appears in gap regions, so that line needs to be ignored.
#        self.data = [d for d in data[data_header_index + 1:] if d[13] != '!']
#
#        #self.read_sec_structure()
#        #self.read_fasta()
#        #self.read_resnum()
#
#    def read_sec_structure(self):
#        """
#        Read secondary structure information from DSSP files.
#
#        Assigns :attr:`ss`.
#        """
#        ss_tmp = list_index_to_array(self.data, index=16)
#
#        if self.reduced:
#            ss_tmp = np.char.translate(ss_tmp, DEFS.dssp_trans)
#
#        self._confirm_ss_data(ss_tmp)
#        self.ss = ss_tmp
#
#    def read_fasta(self):
#        """
#        Read FASTA (primary sequence) information from DSSP data.
#
#        Assigns :attr:`fasta`.
#        """
#        self.fasta = list_index_to_array(self.data, index=13)
#
#    def read_resnum(self):
#        """Read residue numbers, consider only pdb residue numbers."""
#        self.resnums = list(
#            map(
#                str.strip,
#                list_index_to_array(self.data, start=6, stop=10),
#                )
#            )
#
#    @classmethod
#    def from_data_id_tuple(cls, subcmd_tuple, **kwargs):
#        """
#        Initiate from a tuple.
#
#        Tuple of type (PDBID, DSSP DATA).
#        """
#        return cls(
#            pdbid=PDBIDFactory(subcmd_tuple[0]),
#            data=subcmd_tuple[1],
#            **kwargs,
#            )
#
#    @staticmethod
#    def _confirm_ss_data(data):
#        """Confirm secondary structure data characters are valid."""
#        if not all((i in DEFS.dssp_ss_keys.valid for i in data)):
#            raise EXCPTS.DSSPSecStructError()
#
#    @staticmethod
#    def _finds_data_index(data):
#        """
#        Find index for where item startswith '#'.
#
#        Evaluates after left striping white spaces.
#        """
#        # brute force index finding
#        i = 0
#        line = data[i]
#        # if not exausted, while breaks with IndexError
#        while not line.lstrip().startswith('#'):
#            i += 1
#            line = data[i]
#        else:
#            return i
#
#
