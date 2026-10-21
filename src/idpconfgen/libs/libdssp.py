"""
Run the external DSSP program and handle the differences between versions.

IDPConformerGenerator reads the classic DSSP output format and supports:

* DSSP v2 and v3, from https://github.com/cmbi/dssp
* DSSP v4, from https://github.com/PDB-REDO/dssp

The residue columns that IDPConformerGenerator reads are the same in the
classic output of all these versions. The differences handled here are:

* How DSSP is called. DSSP v2 and v3 take the input file with ``-i``.
  Most DSSP v4 releases do not accept ``-i``, and DSSP v4 writes mmCIF
  to the standard output unless it is called with ``--output-format dssp``.
* How PDB files are read. DSSP v4 releases up to v4.4.5 (for example,
  v4.2.2 in Ubuntu 24.04) read a file as PDB only if it starts with a
  HEADER record; otherwise they read it as mmCIF and fail. DSSP v4 builds
  using libcifpp older than 5.1.1 (also v4.2.2 in Ubuntu 24.04) never
  finish reading a PDB file whose last line does not end with a newline,
  and use memory until they are killed. DSSP v4 is therefore given a copy
  of such files with a HEADER line at the top and a newline at the end.
* Polyproline II (PPII) helices. DSSP v4 assigns them itself (code ``P``).
  For DSSP v2 and v3, IDPConformerGenerator assigns them with DSSP-PPII
  (see :mod:`idpconfgen.cli_dssppii`).
"""
import os
import re
import shlex
import shutil
import signal
import subprocess
import tempfile
from contextlib import contextmanager

from idpconfgen.core.exceptions import DSSPExecutionError, DSSPParserError


# Header of the residue table, the same in DSSP v2, v3 and v4.
DSSP_RESIDUE_HEADER = '  #  RESIDUE AA STRUCTURE'

# Column of the secondary structure code in the lines of the residue table.
DSSP_SS_COLUMN = 16

# Blank HEADER record, 80 columns wide, for PDB files given to DSSP v4.
PDB_HEADER_LINE = b'HEADER'.ljust(80) + b'\n'

_DSSP_VERSION_REGEX = re.compile(r'\bversion\s+(\d+(?:\.\d+)+)')


def parse_dssp_version(text):
    """
    Read the DSSP version from the output of ``mkdssp --version``.

    Parameters
    ----------
    text : str
        The output of ``mkdssp --version``.

    Returns
    -------
    tuple of int, or None
        The version, for example ``(4, 6, 1)``. ``None`` if no version
        is found in `text`.
    """
    match = _DSSP_VERSION_REGEX.search(text)
    if match is None:
        return None
    return tuple(int(number) for number in match.group(1).split('.'))


def dssp_version_str(dssp_version):
    """
    Represent a DSSP version as a string.

    Parameters
    ----------
    dssp_version : tuple of int, or None
        The DSSP version, as returned by :func:`get_dssp_version`.

    Returns
    -------
    str
        For example ``'4.6.1'``, or ``'unknown'`` if `dssp_version`
        is ``None``.
    """
    if dssp_version is None:
        return 'unknown'
    return '.'.join(str(number) for number in dssp_version)


def is_dssp4(dssp_version):
    """
    Check whether a DSSP version is DSSP v4 or later.

    Parameters
    ----------
    dssp_version : tuple of int, or None
        The DSSP version, as returned by :func:`get_dssp_version`.

    Returns
    -------
    bool
        ``True`` for DSSP v4 or later. ``False`` for older versions and
        when the version is not known (``None``).
    """
    return dssp_version is not None and dssp_version[0] >= 4


def split_dssp_cmd(dssp_cmd):
    """
    Split a DSSP command into a list of arguments.

    Parameters
    ----------
    dssp_cmd : str or Path
        The DSSP executable, or a command that runs it. For example,
        ``'mkdssp'``, ``'/opt/dssp/bin/mkdssp'`` or
        ``'singularity exec dssp.sif mkdssp'``.

    Returns
    -------
    list of str
        ``[dssp_cmd]`` if `dssp_cmd` is an executable, so that paths with
        spaces work. Otherwise, `dssp_cmd` split as a shell splits it.

    Raises
    ------
    DSSPExecutionError
        If `dssp_cmd` is empty.
    """
    dssp_cmd = os.fspath(dssp_cmd)
    if shutil.which(dssp_cmd):
        return [dssp_cmd]

    args = shlex.split(dssp_cmd)
    if not args:
        raise DSSPExecutionError('the DSSP command is empty.')
    return args


def _run(args):
    """Run a command and capture its output."""
    try:
        return subprocess.run(
            args,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            check=False,
            )
    except OSError as err:
        raise DSSPExecutionError(
            f'could not run `{shlex.join(args)}` ({err}). Make sure DSSP '
            'is installed, or give the DSSP executable with `-cmd`.'
            ) from err


def _exit_status(returncode):
    """Describe the exit status of a command."""
    if returncode >= 0:
        return f'exited with code {returncode}'

    try:
        name = signal.Signals(-returncode).name
    except ValueError:
        name = f'signal {-returncode}'

    status = f'was killed by {name}'
    if name == 'SIGKILL':
        status += ', usually because the system ran out of memory'
    return status


def get_dssp_version(dssp_cmd='mkdssp'):
    """
    Get the version of DSSP.

    Runs ``<dssp_cmd> --version``.

    Parameters
    ----------
    dssp_cmd : str or Path, optional
        The DSSP executable, or a command that runs it.
        See :func:`split_dssp_cmd`. Defaults to ``'mkdssp'``.

    Returns
    -------
    tuple of int, or None
        The DSSP version, for example ``(3, 1, 4)`` or ``(4, 6, 1)``.
        ``None`` if the version cannot be read from the output.

    Raises
    ------
    DSSPExecutionError
        If `dssp_cmd` cannot be run.
    """
    result = _run(split_dssp_cmd(dssp_cmd) + ['--version'])
    output = result.stdout + b'\n' + result.stderr
    return parse_dssp_version(output.decode('utf-8', errors='replace'))


def make_dssp_args(pdb_file, dssp_cmd='mkdssp', dssp_version=None):
    """
    Make the command line arguments to run DSSP on a PDB file.

    DSSP v4 is called with ``--output-format dssp`` so that it writes the
    classic DSSP format. Older and unknown versions are called with
    ``-i``, as DSSP v2 and v3 expect.

    Parameters
    ----------
    pdb_file : str or Path
        The PDB file.

    dssp_cmd : str or Path, optional
        The DSSP executable, or a command that runs it.
        See :func:`split_dssp_cmd`. Defaults to ``'mkdssp'``.

    dssp_version : tuple of int, optional
        The DSSP version, as returned by :func:`get_dssp_version`.

    Returns
    -------
    list of str
    """
    args = split_dssp_cmd(dssp_cmd)
    if is_dssp4(dssp_version):
        args.extend(('--output-format', 'dssp'))
    else:
        args.append('-i')
    args.append(os.fspath(pdb_file))
    return args


@contextmanager
def pdb_for_dssp4(pdb_file):
    """
    Provide a PDB file that all DSSP v4 releases can read.

    DSSP v4 releases up to v4.4.5 read a file as PDB only if it starts
    with a HEADER record. DSSP v4 builds using libcifpp older than 5.1.1
    never finish reading a PDB file whose last line does not end with a
    newline, which is the case for the chain PDB files written by
    ``idpconfgen pdbdl``.

    Parameters
    ----------
    pdb_file : str or Path
        The PDB file.

    Yields
    ------
    str or Path
        `pdb_file` if it starts with a HEADER record and ends with a
        newline. Otherwise, the path to a temporary copy of `pdb_file`
        with a blank HEADER line added at the top and a newline added at
        the end where missing. The copy is deleted when the context exits.
    """
    with open(pdb_file, 'rb') as fin:
        pdb_data = fin.read()

    has_header = pdb_data.startswith(b'HEADER')
    ends_with_newline = pdb_data.endswith(b'\n')
    if has_header and ends_with_newline:
        yield pdb_file
        return

    if not has_header:
        pdb_data = PDB_HEADER_LINE + pdb_data
    if not ends_with_newline:
        pdb_data += b'\n'

    fd, tmp_pdb = tempfile.mkstemp(suffix='.pdb', prefix='idpconfgen_dssp_')
    try:
        with os.fdopen(fd, 'wb') as fout:
            fout.write(pdb_data)
        yield tmp_pdb
    finally:
        os.remove(tmp_pdb)


def check_dssp_output(dssp_output):
    """
    Check that DSSP output is in the classic DSSP format.

    Parameters
    ----------
    dssp_output : bytes
        The output of DSSP.

    Raises
    ------
    DSSPParserError
        If the residue table header is not in `dssp_output`.
    """
    if DSSP_RESIDUE_HEADER.encode() in dssp_output:
        return

    if dssp_output.lstrip().startswith(b'data_'):
        reason = (
            'DSSP wrote mmCIF instead of the classic DSSP format. DSSP v4 '
            'writes the classic format only with `--output-format dssp`.'
            )
    else:
        reason = (
            f'the residue table header ({DSSP_RESIDUE_HEADER.strip()!r}) '
            'is not in the DSSP output.'
            )
    raise DSSPParserError(reason)


def run_dssp(pdb_file, dssp_cmd='mkdssp', dssp_version=None):
    """
    Run DSSP on a PDB file.

    DSSP v4 may be given a copy of `pdb_file` with a HEADER line and a
    final newline added (see :func:`pdb_for_dssp4`).

    Parameters
    ----------
    pdb_file : str or Path
        The PDB file.

    dssp_cmd : str or Path, optional
        The DSSP executable, or a command that runs it.
        See :func:`split_dssp_cmd`. Defaults to ``'mkdssp'``.

    dssp_version : tuple of int, optional
        The DSSP version, as returned by :func:`get_dssp_version`.
        If not given, it is read by running ``<dssp_cmd> --version``.
        Give it when running DSSP on many files.

    Returns
    -------
    bytes
        The output of DSSP in the classic DSSP format.

    Raises
    ------
    DSSPExecutionError
        If DSSP cannot be run or exits with an error.

    DSSPParserError
        If the output of DSSP is not in the classic DSSP format.
    """
    if dssp_version is None:
        dssp_version = get_dssp_version(dssp_cmd)

    if is_dssp4(dssp_version):
        with pdb_for_dssp4(pdb_file) as dssp_input:
            args = make_dssp_args(dssp_input, dssp_cmd, dssp_version)
            result = _run(args)
    else:
        args = make_dssp_args(pdb_file, dssp_cmd, dssp_version)
        result = _run(args)

    if result.returncode != 0:
        # DSSP v4 can write warnings before the error, keep the last lines
        stderr = result.stderr.decode('utf-8', errors='replace').strip()
        stderr = '\n'.join(stderr.splitlines()[-10:])
        details = f':\n{stderr}' if stderr else '.'
        raise DSSPExecutionError(
            f'{os.fspath(pdb_file)}: `{shlex.join(args)}` '
            f'{_exit_status(result.returncode)}{details}'
            )

    check_dssp_output(result.stdout)
    return result.stdout


def remove_native_ppii(dssp_lines):
    """
    Reset the polyproline II helices assigned by DSSP v4 to loops.

    DSSP v4 assigns polyproline II (``P``) only to loop residues, after
    all other secondary structures. Resetting ``P`` to loop (a blank)
    gives the secondary structure that DSSP v2 and v3 assign, to which
    PPII can then be added with DSSP-PPII.

    Parameters
    ----------
    dssp_lines : list of str
        The lines of DSSP output in the classic DSSP format.

    Returns
    -------
    list of str
        The lines with ``P`` replaced by a blank in the secondary
        structure column of the residue table. Other lines and columns
        are not changed.
    """
    col = DSSP_SS_COLUMN
    new_lines = []
    in_residue_table = False
    for line in dssp_lines:
        if in_residue_table and line[col:col + 1] == 'P':
            line = line[:col] + ' ' + line[col + 1:]
        elif line.startswith(DSSP_RESIDUE_HEADER):
            in_residue_table = True
        new_lines.append(line)
    return new_lines
