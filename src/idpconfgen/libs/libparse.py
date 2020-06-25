"""
Parsing routines for different data structure.

All functions in this module receive a certain Python native datastructure,
parse the information inside and return/yield the parsed information.
"""
import subprocess

from idpconfgen import Path
from idpconfgen.core.definitions import dssp_trans, jsonparameters
from idpconfgen.libs.libpdb import PDBIDFactory, atom_resSeq


# _ascii_lower_set = set(string.ascii_lowercase)
# _ascii_upper_set = set(string.ascii_uppercase)


# USED OKAY
# all should return a string
# used for control flow
type2string = {
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

    """
    # begin
    pdbfile = pdb.resolve()
    _cmd = [cmd, '-i', str(pdbfile)]
    result = subprocess.run(_cmd, capture_output=True)
    yield from split_pdb_by_dssp(
        pdbfile,
        result.stdout.decode('utf-8'),
        **kwargs)


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

        residues_bytes = set(c.encode() for c in residues)

        lines2write = [
            line for line in pdb_bytes
            if line[atom_resSeq].strip() in residues_bytes
            ]

        if all(line.startswith(b'HETATM') for line in lines2write):
            continue

        __data = {
            _ss: dssp,
            _fasta: fasta,
            _resids: ','.join(residues),
            }

        yield fname, __data, b''.join(lines2write)
        segs += 1


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
            print(line)
            break
    else:
        # if the whole generator is exhausted
        raise IndexError

    dssp = []
    fasta = []
    residues = []
    _d_append = dssp.append
    _f_append = fasta.append
    _r_append = residues.append
    for line in RM1:

        if line[13:14] != '!':
            _d_append(line[16:17])
            _f_append(line[13:14])
            _r_append(line[6:10].strip())

        else:
            dssp_ = ''.join(dssp).translate(DT) if reduced else ''.join(dssp)

            yield dssp_, ''.join(fasta), residues

            dssp.clear()
            fasta.clear()
            residues.clear()  # Attention with this clear!

    else:
        dssp_ = ''.join(dssp).translate(DT) if reduced else ''.join(dssp)
        yield dssp_, ''.join(fasta), residues
