"""
Extracts backbone torsion angles from PDBs.

PROTOCOL:

1. Reads backbone coordinates (N, CA, C) from PDB/mmCIF files.
2. Calculates torsion angles from the backbone.
3. Saves results to a JSON dictionary where keys are the input file
    names and the value is a dictionary containing three lists: 'OMEGA',
    'PHI', and 'PSI'.
4. If `source` JSON file is given, updates that file with the new
    information. Preexisting keys are deleted.

CONTROLLED CHECKS:

For each PDB/mmCIF analyzed, fails if:

1. The number of N, CA, and C atoms differ.
2. There are unexpected inconsistencies in the PDB/mmCIF files.
3. Any of the consecutive atoms are more than 2.1A apart.

Failed PDBs are registered in `.rpr_on_crash` files and ignored.

USAGE:
    $ idpconfgen torsions [PDBS]
    $ idpconfgen torsions [PDBS] -sc file.json
    $ idpconfgen torsions [PDBS] -sc file.json -n
    $ idpconfgen torsions [PDBS] -sc file.json -o mytorsions.json -n
    $ idpconfgen torsions [PDBS] -sc file.json -o mytorsions.json -n -deg
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
libcli.add_argument_source(ap)
libcli.add_argument_output(ap)
libcli.add_argument_degrees(ap)
libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        source=None,
        output=None,
        degrees=True,
        func=None,
        ncores=1,
        ):
    """Perform main script logic."""
    # validates before performing time consuming calculations
    if source and not source.endswith('.json'):
        raise ValueError('Source file should have `.json` extension.')

    output = output or 'torsions.json'
    if not output.endswith('.json'):
        raise ValueError('Output file should have `.json` extension.')

    log.info(T('Extracting torsion angles'))
    init_files(log, LOGFILESNAME)

    if source:
        database_dict = read_dictionary_from_disk(source)

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

    if source:

        pop_difference_with_log(database_dict, torsion_result)

        for key, value in torsion_result.items():
            # where value is a dictionary {'chi':, 'phi':, 'omega':}
            database_dict[key].update(value)

        save_dict_to_json(database_dict, output=output)

    else:
        save_dict_to_json(torsion_result, output=output)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
