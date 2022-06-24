"""# noqa: D205, D210, D400, E501
Interface for X-EISD logic to score and reweight conformer ensembles.
"""
import argparse
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    FileReaderIterator,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libparse import pop_difference_with_log
from idpconfgen.logger import S, T, init_files, report_on_crash

from idpconfgen.components.eisd import (
    eisd_modes,
    eisd_modules,
    eisd_optimization_types,
    default_type,
    default_mode,
    )

LOGFILESNAME = '.idpconfgen_eisd'

_name = 'eisd'
_help = 'Score and reweight generated ensembles.'  # noqa: E501

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-fpath',
    '--filepath',
    help=('Path to the folder containing experimental'
          'and back-calculation data files.'),
    type=str,
    required=True,
    default=None,
    )

ap.add_argument(
    '-eps',
    '--epochs',
    help='Number of epochs to run main optimization, defaults to 1000.',
    type=int,
    default=1000,
    required=True,
)

ap.add_argument(
    '-opt',
    '--optimization',
    help='If optimization is required, defaults to True.',
    action='store_true',
    required=False,
)

ap.add_argument(
    '-m',
    '--mode',
    help='Optimization mode on all, single, or dual experimental observables.',
    choices=eisd_modes,
    default=default_mode,
    required=False,
)

ap.add_argument(
    '--optimization-type',
    help='Type of optimzation done. Defaults to max.',
    choices=eisd_optimization_types,
    default=default_type,
    required=False,
)

libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)


def main():
    return



if __name__ == '__main__':
    libcli.maincli(ap, main)