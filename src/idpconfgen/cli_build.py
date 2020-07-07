"""
Builds.

USAGE:
    $ idpconfgen build DB

"""
import argparse
import re
import random
from functools import partial

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libfilter import (
    aligndb,
    index_single_regex_overlap,
    )
from idpconfgen.libs.libtimer import timeme


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



def main(database, func=None):
    """."""

    db = read_dictionary_from_disk(database)
    ldb = len(db)
    _log = f'Read DB with {ldb} entries'
    log.info(_log)

    timed = partial(timeme, aligndb)
    pdbs, angles, dssp, resseq = timed(db)

    # esta parte tengo que ponerla de parametro externo
    regex = re.compile(r'(?=(L{6}))')
    timed = partial(timeme, index_single_regex_overlap)
    loops_6 = timed(dssp, regex)
    _log = f'Found {len(loops_6)} indexes for loops of 6'
    log.info(_log)


    from time import time
    start = time()
    for i in range(500_000):
        agl = random.choice(loops_6)
        agls = angles[agl, :]
    print(time() - start)

    return



if __name__ == "__main__":
    libcli.maincli(ap, main)
