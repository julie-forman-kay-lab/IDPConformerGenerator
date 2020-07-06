"""
Builds.

USAGE:
    $ idpconfgen build DB

"""
import argparse

from idpconfgen.libs import libcli
from idpconfgen.libs.libio import read_dictionary_from_disk


_name = 'build'
_help = 'Builds conformers from database.'


_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    'database',
    help='The IDPConfGen database.',
    )



def main(database):
    """."""

    db = read_dictionary_from_disk(database)

    return



if __name__ == "__main__":
    libcli.climain(ap, main)
