"""
Client for renaming 3 letter residues in a PDB file.

USAGE:
    $ idpconfgen resptm <PDB-FILES> -pt <PATTERN>
    $ idpconfgen resptm <PDB-FILES> -pt <PATTERN> -o <OUTPUT> -n <CORES>
"""
import argparse
import json
import os
import shutil
from functools import partial
from pathlib import Path

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.idpconfgen_resptm'
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
        "Pattern separated by colons for residue number and desired new "
        "residue name. Additional residues are delimited by commas, "
        "e.g. -pt '12:TPO,14:SEP'"
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

def main(
        pdb_files,
        pattern,
        output,
        ncores=1,
        **kwargs,
        ):
    pass


if __name__ == '__main__':
    libcli.maincli(ap, main)
