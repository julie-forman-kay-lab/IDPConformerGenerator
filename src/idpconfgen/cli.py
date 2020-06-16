"""
IDP CONFORMER GENERATOR.

Generates IDP Conformers by combining observed backbone angles
    and random distribution.

USAGE:
    For help:
    >>> idpconfgen -h

    Using PDB Downloader:
    >>> idpconfgen pdbdl

"""
import argparse
import sys
import tracemalloc

from idpconfgen import cli_pdbdownloader as pdbdl
from idpconfgen import cli_fetch
from idpconfgen import cli_filter as cli_filter
from idpconfgen import cli_bbsplit as cli_bbsplit
from idpconfgen import cli_segext as segext
from idpconfgen import cli_sscalc as sscalc
from idpconfgen import cli_fastaext as fastaext
from idpconfgen import cli_RCSB_dssp as cli_RCSB_dssp
from idpconfgen.libs import libcli


# https://stackoverflow.com/questions/14988397
# https://stackoverflow.com/questions/4042452
#def _load_args():

prog_, description_, usage_ = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog='idpconfgen',  # prog_,
    description=libcli.detailed.format(description_),
    usage=usage_,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

subparsers = ap.add_subparsers(
    title='SUBROUTINES',
    description='DESCRIPTION',
    help='IDP Conf Gen subroutines:',
    )

libcli.add_version(ap)
libcli.add_subparser(subparsers, cli_fetch)
libcli.add_subparser(subparsers, pdbdl)
libcli.add_subparser(subparsers, cli_filter)
libcli.add_subparser(subparsers, fastaext)
libcli.add_subparser(subparsers, sscalc)
libcli.add_subparser(subparsers, segext)
libcli.add_subparser(subparsers, cli_bbsplit)
libcli.add_subparser(subparsers, cli_RCSB_dssp)


def load_args():
    """Load user input arguments."""
    return ap.parse_args()


def maincli():
    """
    Execute subroutine.

    Arguments are read from user command line input.
    """
    # prints help if not arguments are passed
    # if >2 prints help of subroutines.
    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    cmd = load_args()

    tracemalloc.start()
    try:
        cmd.func(**vars(cmd))
    except KeyboardInterrupt:
        pass
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')
    print('top 10')
    for stat in top_stats[:10]:
        print(stat)


if __name__ == '__main__':
    maincli()
