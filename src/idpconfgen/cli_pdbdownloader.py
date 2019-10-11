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

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = argparse.ArgumentParser(
    prog=_prog,
    description=_des,
    usage=_us,
    )
# https://stackoverflow.com/questions/24180527 

ap.add_argument(
    'cull',
    help='The CULL PDB file list.',
    type=Path,
    )


def load_args():
    if len(sys.argv) == 1:
        ap.print_help()
        ap.exit()
    cmd = ap.parse_args()
    return cmd


def main(args):
    print('I am the main of pdbdownloader')
    print(args)


def maincli():
    cmd = load_args()
    main(cmd)


if __name__ == '__main__':
    maincli()
