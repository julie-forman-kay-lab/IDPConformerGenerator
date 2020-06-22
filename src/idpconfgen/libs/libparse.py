"""
Parsing routines for different data structure.

All functions in this module receive a certain Python native datastructure,
parse the information inside and return/yield the parsed information.
"""
import string
import subprocess
import sys
import traceback
from contextlib import contextmanager
from os import fspath

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.core.definitions import dssp_trans
from idpconfgen.libs import libcheck
from idpconfgen.libs.libpdb import PDBIDFactory, delete_insertions, atom_resSeq
from idpconfgen.logger import S
from idpconfgen.libs.libio import save_pairs_to_disk


_ascii_lower_set = set(string.ascii_lowercase)
_ascii_upper_set = set(string.ascii_uppercase)


# USED OKAY
# all should return a string
# used for control flow
type2string = {
    type(Path()): lambda x: x.read_text(),
    bytes: lambda x: x.decode('utf_8'),
    str: lambda x: x,
    }


# USED OKAY
def parse_dssp(data, reduced=False):
    """
    Parse DSSP file data.

    JSON doesn't accept bytes
    That is why `data` is expected as str.
    """
    DT = dssp_trans

    data_ = data.split('\n')

    # RM means removed empty
    RM1 = (i for i in data_ if i)

    # exausts generator until
    for line in RM1:
        if line.strip().startswith('#'):
            break
    else:
        # if the whole generator is exhausted
        raise IndexError

    dssp = []
    fasta = []
    residues = []
    for line in RM1:
        if line[13:14] != '!':
            dssp.append(line[16:17])
            fasta.append(line[13:14])
            residues.append(line[6:10].strip())
        else:
            dssp_ = ''.join(dssp)
            if reduced:
                dssp_ = dssp_.translate(DT)
            yield dssp_, ''.join(fasta), residues
            dssp = []
            fasta = []
            residues = []
    else:
        dssp_ = ''.join(dssp)
        if reduced:
            dssp_ = dssp_.translate(DT)
        yield dssp_, ''.join(fasta), residues


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
        for i in range(len(chain) + 1):
            cl = f'{chain[:i]}{chain[i:].lower()}'
            if cl in chain_set:
                return cl

    raise ValueError(f'Neither {chain}/{cl} are in set {chain_set}')



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


# NOT USED
def mkdssp_w_split_w_save(*args, destination=None, **kwargs):
     for fname, split_dssp, split_pdb_bytes in mkdssp_w_split(*args, **kwargs):
        save_pairs_to_disk(((f'{fname}.pdb', split_pdb_bytes),), destination=destination)
        yield fname, split_dssp


def mkdssp_w_split(
        pdb,
        cmd,
        minimum=2,
        reduced=False,
        ):
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

    minimum : int
        The minimum length allowed for a segment.

    reduce : bool
        Whether to reduce the DSSP nomemclature to H/E/L.
    """
    pdbfile = fspath(pdb.resolve())
    _cmd = [cmd, '-i', pdbfile]
    result = subprocess.run(_cmd, capture_output=True)

    # did not use the pathlib interface read_bytes on purpose.
    with open(pdbfile, 'rb') as fin:
        pdb_bytes = fin.readlines()

    segs = 0
    dssp_text = result.stdout.decode('utf-8')
    pairs = []
    for dssp, fasta, residues in parse_dssp(dssp_text, reduced=reduced):
        if len(residues) < minimum:
            continue

        fname =f'{str(PDBIDFactory(pdb))}_seg{segs}'

        residues_bytes = set(c.encode() for c in residues)

        lines2write = [
                line for line in pdb_bytes
                if line[atom_resSeq].strip() in residues_bytes
                ]

        if all(l.startswith(b'HETATM') for l in lines2write):
            continue

        __data = {
            'dssp': dssp,
            'fasta': fasta,
            'resids': ','.join(residues),
            }

        yield fname, __data , b''.join(lines2write)
        segs += 1
