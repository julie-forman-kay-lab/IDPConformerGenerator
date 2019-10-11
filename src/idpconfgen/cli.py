"""
IDP CONFORMER GENERATOR

Calculates IDP Confomers
    from a database of angles
        and other stuff.

USAGE:
    Do this.
    >>> lalalal

"""
import argparse

from idpconfgen.cli_pdbdownloader import (
    ap as pdbdlap,
    main as pdbdlmain,
    )


# https://stackoverflow.com/questions/14988397/
def load_args():
    
    doclines = __doc__.split('\n')
    prog_ = doclines[1]
    des_ = '\n'.join(doclines[3:doclines.index('USAGE:')])
    usage_ = '\n' + '\n'.join(doclines[doclines.index('USAGE:') + 1:])

    ap = argparse.ArgumentParser(
         prog=prog_,
         description=des_,
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

    cmd = ap.parse_args()

    return cmd
    

def main(cmd):
    
    cmd.func(cmd)    

if __name__ == '__main__':
    
    cmd = load_args()
    main(cmd)
