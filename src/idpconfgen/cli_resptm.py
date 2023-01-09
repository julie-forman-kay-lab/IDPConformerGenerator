"""
Client for renaming 3 letter residues in a PDB file.

USAGE:
    $ idpconfgen resptm <PDB-FILES> -pt <PATTERN>
    $ idpconfgen resptm <PDB-FILES> -pt <PATTERN> -o <OUTPUT> -n <CORES>
"""
import argparse
import shutil
import re
from functools import partial
from pathlib import Path

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.idpconfgen_resptm'
TMPDIR = '__tmpresptm__'

_name = 'resptm'
_help = 'Renames residues of interest within a PDB file for post-translational modifications.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-pt',
    '--pattern',
    help=(
        "Entire pattern in single quotes. Pattern is denoted by colons "
        "for residue number and desired new residue name. Additional residues "
        "are delimited by commas and pattern must end at a comma. "
        "For e.g. -pt '12:TPO,14:SEP,'"
        ),
    nargs='?',
    )

ap.add_argument(
    '-o',
    '--output',
    help=(
        "A path to the folder where renamed residue PDBs are stored. "
        "Defaults to overwrite working files."
        ),
    type=Path,
    )

libcli.add_argument_ncores(ap)

ap.add_argument(
    '-tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )

def main(
        pdb_files,
        pattern,
        output,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    init_files(log, LOGFILESNAME)
    REGEX_PATTERN = re.compile("[0-9]*:[a-z]{3},")
    RESNUM = []
    RESNAME = []
    
    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))
    
    log.info(T('reading input pattern'))
    if REGEX_PATTERN.match(pattern):
        residues = pattern.split(',')
        residues.pop()  # last element should be empty
        for res in residues:
            tmp = res.split(':')
            RESNUM.append(tmp[0])
            RESNAME.append(tmp[1])
    else:
        log.info(S('Incorrect pattern input.'))
        log.info(S('Pattern is as follows: 1:ABC,...'))
        return
    log.info(S('done'))
    
    if _istarfile:
        shutil.rmtree(tmpdir)


if __name__ == '__main__':
    libcli.maincli(ap, main)
