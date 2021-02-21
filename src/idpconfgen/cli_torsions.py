"""
Extracts torsion angles from PDBs.

USAGE:
    $ idpconfgen torsions [PDBS]
"""
import argparse
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libhigherlevel import cli_helper_calc_torsions
from idpconfgen.libs.libio import (
    FileReaderIterator,
    read_dictionary_from_disk,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libparse import pop_difference_with_log
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = 'idpconfgen_torsion'

_name = 'torsions'
_help = 'Calculate torsion angles for PDB files.'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

libcli.add_argument_pdb_files(ap)
libcli.add_argument_db(ap)
libcli.add_argument_degrees(ap)
libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        database,
        degrees=True,
        func=None,
        ncores=1,
        ):
    """Perform main script logic."""
    log.info(T('Extracting torsion angles'))
    init_files(log, LOGFILESNAME)

    database_dict = read_dictionary_from_disk(database)

    log.info(T('reading input paths'))
    pdbs = FileReaderIterator(pdb_files, ext='.pdb')
    log.info(S('done'))

    consume = partial(starunpack, cli_helper_calc_torsions, degrees=degrees, decimals=10)

    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_prefix=_name,
        )

    execute_pool = pool_function(execute, pdbs, ncores=ncores)

    torsion_result = {
        Path(pdbid).stem: angles
        for pdbid, angles in execute_pool
        }

    pop_difference_with_log(database_dict, torsion_result)

    for key, value in torsion_result.items():
        # where value is a dictionary {'chi':, 'phi':, 'omega':}
        database_dict[key].update(value)

    save_dict_to_json(database_dict, output='torsions.json')
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
