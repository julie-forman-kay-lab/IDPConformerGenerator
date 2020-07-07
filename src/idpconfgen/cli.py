"""
IDP CONFORMER GENERATOR.

Generates IDP Conformers by combining observed backbone angles
    and random distribution.

USAGE:
    For help:
    >>> idpconfgen -h

    Using PDB Downloader:
    >>> idpconfgen cli_pdbdl

"""
import argparse
import sys

from idpconfgen import cli_build, cli_fastaext, cli_fetch
from idpconfgen import cli_pdbdownloader as cli_pdbdl
from idpconfgen import cli_ssext, cli_sscalc, cli_torsions, log
from idpconfgen.libs import libcli
from idpconfgen.logger import S


prog_, description_, usage_ = libcli.parse_doc_params(__doc__)

description = f"""
{description_}

Individual routines for DB creation:

    * {cli_pdbdl._name}
    * {cli_sscalc._name}
    * {cli_torsions._name}
    * {cli_build._name}

Other useful routines:

    * {cli_fetch._name}
    * {cli_ssext._name}
    * {cli_fastaext._name}
"""

ap = libcli.CustomParser(
    prog='idpconfgen',  # prog_,
    description=libcli.detailed.format(description),
    usage=usage_,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version(ap)

# routines from the main DB creation pipeline
subparsers = ap.add_subparsers(
    title='IDP Conformer Generator routines',
    help='Short description:',
    )

# argument parsers for main DB creation routines
libcli.add_subparser(subparsers, cli_pdbdl)
libcli.add_subparser(subparsers, cli_sscalc)
libcli.add_subparser(subparsers, cli_torsions)
libcli.add_subparser(subparsers, cli_build)

# argument parsers for secondary routines
libcli.add_subparser(subparsers, cli_fastaext)
libcli.add_subparser(subparsers, cli_fetch)
libcli.add_subparser(subparsers, cli_ssext)


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
    cmd.func(**vars(cmd))
    log.info(S('finished properly'))


if __name__ == '__main__':
    maincli()
