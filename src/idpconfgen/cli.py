"""
IDP CONFORMER GENERATOR

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

from idpconfgen.libs import libcli

from idpconfgen.cli_pdbdownloader import (
    ap as pdbdlap,
    main as pdbdlmain,
    )

# https://stackoverflow.com/questions/14988397
# https://stackoverflow.com/questions/4042452
def load_args():
    
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
    
    ap_pdbdl = subparsers.add_parser(
        'pdbdl',
        help='PDB Downloader',
        parents=[pdbdlap],
        add_help=False,
        description=libcli.detailed.format(pdbdlap.description),
        usage=pdbdlap.usage,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    ap_pdbdl.set_defaults(func=pdbdlmain)
   
    # prints help if not arguments are passed
    # if >2 prints help of subroutines.
    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    cmd = ap.parse_args()

    return cmd
    

def maincli():
    cmd = load_args() 
    cmd.func(**vars(cmd))


if __name__ == '__main__':
    
    maincli()
