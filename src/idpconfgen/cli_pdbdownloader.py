"""
PDBDOWNLOADER

Downloads structures from RCSB Databank.

USAGE:
    >>> icgpdbdl CULL

"""
import argparse
from pathlib import Path
import sys

from idpconfgen.libs import libcli
from idpconfgen.libs import libpdb

LOGFILESNAME = 'pdbdownloader'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = argparse.ArgumentParser(
    prog=_prog,
    description=_des,
    usage=_us,
    )
# https://stackoverflow.com/questions/24180527 

ap.add_argument(
    'pdblist',
    help='A list of PDBID:CHAIN to download.',
    type=Path,
    )

ap.add_argument(
    '-d',
    '--destination',
    help=(
        'Destination folder where PDB files will be stored.'
        'Defaults to current working directory.'
        )
    type=Path,
    default=Path.cwd(),
    )

ap.add_argument(
    '-u',
    '--update',
    help='Updates destination folder according to input PDB list.',
    action='store_true',
    )

ap.add_argument(
    '-rn',
    '--record_name',
    help='The coordinate PDB record name.',
    detault=('ATOM',),
    type=tuple,
    )

ap.add_argument(
    '-n',
    '--ncores',
    help='Number of cores to use.',
    type=int,
    default=1,
    )


def load_args():
    if len(sys.argv) == 1:
        ap.print_help()
        ap.exit()
    cmd = ap.parse_args()
    return cmd


def main(
        pdblist,
        destination=None,
        update=False,
        record_name=('ATOM',),
        ncores=1,
        ):
    
    logger.init_files(log, LOGFILESNAME)
    
    with pdblist.open('r') as fh:
        pdblist = libpdb.PDBList(fh.readlines())

    print('I am the main of pdbdownloader')
    print(args)


def maincli():
    cmd = load_args()
    main(cmd)


if __name__ == '__main__':
    maincli()
