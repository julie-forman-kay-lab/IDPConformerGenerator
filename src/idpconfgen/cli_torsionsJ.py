"""
Calculates J-torsion angles from PDB/mmCIF files.

Calculated torsion angles are between HN-CaHA.

The calculated angles will be stored in a JSON file, one entry per
structure, named 'torsionsJ_HN-CaHA.json'.

Structures are allowed to have unordered atoms as long as residues are
ordered.

If input structures lack protons at the N-terminal the flag `--hn_term`
should be set to 0, fails otherwise.

USAGE:
    $ idpconfgen torsionsJ [PDBS]
    $ idpconfgen torsionsJ file1.pdb file2.pdb
    $ idpconfgen torsionsJ folder_with_pdbs/*.pdb
    $ idpconfgen torsionsJ pdbs_file_list.list
    $ idpconfgen torsionsJ folder_with_cifs/*.cif
    $ idpconfgen torsionsJ [PDBS] -deg
    $ idpconfgen torsionsJ [PDBS] -deg -dec 3
    $ idpconfgen torsionsJ [PDBS] -deg -dec 3 --hnterm 0
"""
import argparse
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libhigherlevel import cli_helper_calc_torsionsJ
from idpconfgen.libs.libio import (
    FileReaderIterator,
    read_dictionary_from_disk,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libparse import pop_difference_with_log
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = 'idpconfgen_torsionJ'

_name = 'torsionsJ'
_help = 'Calculate HN-CaHA torsion angles from PDB files.'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)
libcli.add_argument_degrees(ap)
libcli.add_argument_nohnterm(ap)
libcli.add_argument_decimals(ap)
libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        degrees=True,
        decimals=3,
        no_hn_term=True,
        func=None,
        ncores=1,
        ):
    """Perform main script logic."""
    log.info(T('Extracting torsion angles'))
    init_files(log, LOGFILESNAME)

    log.info(T('reading input paths'))
    pdbs = FileReaderIterator(pdb_files, ext='.pdb')
    log.info(S('done'))

    consume = partial(
        cli_helper_calc_torsionsJ,
        degrees=degrees,
        decimals=decimals,
        hn_terminal=not no_hn_term,
        )

    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_prefix=_name,
        )

    execute_pool = pool_function(execute, pdbs, ncores=ncores)

    torsion_result = {
        pdbid.stem: list(angles)
        for pdbid, angles in execute_pool
        }

    save_dict_to_json(torsion_result, output='torsionsJ_HN-CaHA.json')
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
