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

from idpconfgen import cli_pdbdownloader as pdbdl
from idpconfgen import cli_pdb_raw_downloader as raw_pdbdl
from idpconfgen import cli_segext as segext
from idpconfgen import cli_ssext as ssext
from idpconfgen import cli_fastaext as fastaext
from idpconfgen.libs import libcli


# https://stackoverflow.com/questions/14988397
# https://stackoverflow.com/questions/4042452
def _load_args():

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

    ap_fastaext = subparsers.add_parser(
        'fastaext',
        help='Extract FASTAS from PDBs',
        parents=[fastaext.ap],
        add_help=False,
        description=libcli.detailed.format(fastaext._des),
        usage=fastaext._us,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    ap_fastaext.set_defaults(func=fastaext.main)

    ap_pdbdl = subparsers.add_parser(
        'pdbdl',
        help='PDB Downloader',
        parents=[pdbdl.ap],
        add_help=False,
        description=libcli.detailed.format(pdbdl._des),
        usage=pdbdl._us,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    ap_pdbdl.set_defaults(func=pdbdl.main)

    ap_ssext = subparsers.add_parser(
        'ssext',
        help='PDB Downloader',
        parents=[ssext.ap],
        add_help=False,
        description=libcli.detailed.format(ssext._des),
        usage=ssext._us,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    ap_ssext.set_defaults(func=ssext.main)

    ap_segext = subparsers.add_parser(
        'segext',
        help='Segment extract',
        parents=[segext.ap],
        add_help=False,
        description=libcli.detailed.format(segext._des),
        usage=segext._us,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    ap_segext.set_defaults(func=segext.main)

    # prints help if not arguments are passed
    # if >2 prints help of subroutines.
    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    cmd = ap.parse_args()

    return cmd


def maincli():
    """
    Execute subroutine.

    Arguments are read from user command line input.
    """
    cmd = _load_args()
    cmd.func(**vars(cmd))


if __name__ == '__main__':

    maincli()
