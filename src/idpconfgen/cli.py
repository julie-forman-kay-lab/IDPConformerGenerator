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
    
    ap = argparse.ArgumentParser(
         prog=prog_,
         description=description_,
         usage=usage_,
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
        )
    ap_pdbdl.set_defaults(func=pdbdlmain)
   
    # prints help if not arguments are passed
    # < 3 to consider also subroutines.
    if len(sys.argv) < 3:
        ap.print_help()
        ap.exit()

    cmd = ap.parse_args()

    return cmd
    

def main():
    cmd = load_args() 
    cmd.func(cmd)    

if __name__ == '__main__':
    
    main()
