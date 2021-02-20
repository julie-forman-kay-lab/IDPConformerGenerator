"""
Parsing routines for different data structure.

All functions in this module receive a certain Python native datastructure,
parse the information inside and return/yield the parsed information.
"""
import subprocess
from itertools import product
from pathlib import Path as Path_

from numba import njit

from idpconfgen import Path, log
from idpconfgen.core.definitions import dssp_trans_bytes, jsonparameters
from idpconfgen.core.exceptions import DSSPParserError
from idpconfgen.libs.libpdb import PDBIDFactory, atom_resSeq
from idpconfgen.logger import S


# _ascii_lower_set = set(string.ascii_lowercase)
# _ascii_upper_set = set(string.ascii_uppercase)


# USED OKAY
# all should return a string
# used for control flow
type2string = {
    type(Path_()): lambda x: x.read_text(),
    type(Path()): lambda x: x.read_text(),
    bytes: lambda x: x.decode('utf_8'),
    str: lambda x: x,
    }


# USED OKAY
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

    datum_hold = prev = data[0]
    start = 0
    groups = []
    GA = groups.append
    for i, datum in enumerate(data):
        if datum != prev:
            GA([datum_hold, slice(start, i)])
            datum_hold, start = datum, i
        prev = datum
    else:
        GA([datum_hold, slice(start, i + 1)])

    assert isinstance(groups[0][0], str)
    assert isinstance(groups[0][1], slice)
    return groups


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


def mkdssp_w_split(pdb, cmd, **kwargs):
    """
    Execute `mkdssp` from DSSP.

    Saves the data splitted accoring to backbone continuity
    as identified by `mkdssp`. Splits the input PDB into bb continuity
    segments.

    https://github.com/cmbi/dssp

    Parameters
    ----------
    pdb : Path
        The path to the pdb file.

    cmd : str
        The command to execute the external DSSP program.

    Yields
    ------
    from `split_pdb_by_dssp`
    """
    pdbfile = pdb.resolve()
    _cmd = [cmd, '-i', str(pdbfile)]
    result = subprocess.run(_cmd, capture_output=True)
    # if subpress.run fails, the parser will raise an error,
    # no need to assert result.returncode
    yield from split_pdb_by_dssp(pdbfile, result.stdout, **kwargs)


def split_pdb_by_dssp(pdbfile, dssp_text, minimum=2, reduced=False):
    """
    Split PDB file based on DSSP raw data.

    Parameters
    ----------
    minimum : int
        The minimum length allowed for a segment.

    reduce : bool
        Whether to reduce the DSSP nomemclature to H/E/L.
    """
    # local scope
    _ss = jsonparameters.ss
    _fasta = jsonparameters.fasta
    _resids = jsonparameters.resids

    # did not use the pathlib interface read_bytes on purpose.
    with open(pdbfile, 'rb') as fin:
        pdb_bytes = fin.readlines()

    segs = 0
    pdbname = str(PDBIDFactory(pdbfile.stem))
    # converted because end product is to be saved in json or used as string
    for dssp, fasta, residues in parse_dssp(dssp_text, reduced=reduced):
        if len(residues) < minimum:
            continue

        fname = f'{pdbname}_seg{segs}'

        residues_set = set(residues)

        lines2write = [
            line for line in pdb_bytes
            if line[atom_resSeq].strip() in residues_set
            ]

        if all(line.startswith(b'HETATM') for line in lines2write):
            continue

        # it is much faster to decode three small strings here
        # than a whole DSSP file at the beginning
        __data = {
            _ss: dssp.decode('utf-8'),
            _fasta: fasta.decode('utf-8'),
            _resids: b','.join(residues).decode('utf-8'),
            }

        yield fname, __data, b''.join(lines2write)
        segs += 1


# USED OKAY
def sample_case(input_string):
    """
    Sample all possible cases combinations from `string`.

    Examples
    --------
    >>> sample_case('A')
    {'A', 'a'}

    >>> sample_case('Aa')
    {'AA', 'Aa', 'aA', 'aa'}
    """
    possible_cases = set(
        map(
            ''.join,
            product(*zip(input_string.upper(), input_string.lower())),
            )
        )
    return possible_cases


# USED OKAY
def parse_dssp(data, reduced=False):
    """
    Parse DSSP file data.

    JSON doesn't accept bytes
    That is why `data` is expected as str.
    """
    DT = dssp_trans_bytes

    data_ = data.split(b'\n')

    # RM means removed empty
    RM1 = (i for i in data_ if i)

    # exausts generator until
    for line in RM1:
        if line.strip().startswith(b'#'):
            break
    else:
        # if the whole generator is exhausted
        raise DSSPParserError("File exhausted without finding '#'")

    dssp = []
    fasta = []
    residues = []
    _d_append = dssp.append
    _f_append = fasta.append
    _r_append = residues.append
    bjoin = b''.join
    for line in RM1:

        if line[13:14] != b'!':
            _d_append(line[16:17])
            _f_append(line[13:14])
            _r_append(line[6:10].strip())

        else:
            dssp_ = bjoin(dssp).translate(DT) if reduced else bjoin(dssp)
            yield dssp_, bjoin(fasta), residues

            dssp.clear()
            fasta.clear()
            residues.clear()  # Attention with this clear!

    else:
        dssp_ = bjoin(dssp).translate(DT) if reduced else bjoin(dssp)
        yield dssp_, bjoin(fasta), residues


def pop_difference_with_log(dict1, dict2):
    """
    Pop keys in `dict1` that are not present in `dict2`.

    Reports pop'ed keys to log INFO.

    Operates `dict1` in place.

    Parameters
    ----------
    dict1, dict2 : dict

    Returns
    -------
    None
    """
    d1k = set(dict1.keys())
    d2k = set(dict2.keys())

    diff = d1k.difference(d2k)
    if diff:
        log.info(S(
            'The following keys will be removed from the dictionary: {}\n',
            diff
            ))

        for key in diff:
            dict1.pop(key)


def remap_sequence(seq, target='A', group=('P', 'G')):
    """
    Remap sequence.

    Parameters
    ----------
    seq : Protein primary sequence in FASTA format.

    target : str (1-char)
        The residue to which all other residues will be converted
        to.

    group : tuple
        The list of residues that excape map/conversion.

    Return
    ------
    str
        The remaped string.

    Examples
    --------

    >>> remap_sequence('AGTKLPHNG')
    'AGAAAPAAG'
    """
    return ''.join(target if res not in group else res for res in seq)


# njit
def get_trimer_seq(seq, idx):
    pre = seq[idx - 1] if idx > 0 else 'G'
    curr_res = seq[idx]
    try:
        pos = seq[idx + 1]
    except:  #IndexError in plain python
        pos = 'G'

    return curr_res, pre + pos


get_trimer_seq_njit = njit(get_trimer_seq)
