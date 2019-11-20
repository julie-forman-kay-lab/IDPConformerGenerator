"""
Extracts secondary structure information from PDBs.

Uses an external third party software.

USAGE:
    $ idpconfgen ssext [PDBS]
"""
import argparse

from idpconfgen import Path, log
from idpconfgen.libs import libcli, libio, libmulticore, libparse
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = 'idpconfgen_ssext'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    'ss_cmd',
    help='The path to the DSSP executable file.',
    type=str,
    )

ap.add_argument(
    'pdbs',
    help='PDB file list.',
    nargs='+',
    )

ap.add_argument(
    '-o',
    '--output',
    help=(
        "The output file containing the PDBID and "
        "respective secondary structure information. "
        "Defaults to sys.stdout (prints to console)."
        ),
    type=Path,
    default=None,
    )

ap.add_argument(
    '-n',
    '--ncores',
    help='Number of cores to use',
    default=1,
    type=int,
    )


def _load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = _load_args()
    main(**vars(cmd))


def main(
        ss_cmd,
        pdbs,
        output=None,
        ncores=1,
        **kwargs,
        ):
    """Run main cli logic."""
    
    init_files(log, LOGFILESNAME)
    
    pdbs_paths = libio.read_path_bundle(pdbs)

    ss_ext_exec = libmulticore.JoinedResults(
        pdbs,
        ss_cmd,
        ncores=ncores,
        TaskMethod=libmulticore.DSSPTask,
        results_parser=libparse.DSSPParser.from_data_id_tuple,
        )

    ss_ext_exec.run()
    libparse.export_ss_from_DSSP(*ss_ext_exec.results)

    return


if __name__ == '__main__':

    maincli()

