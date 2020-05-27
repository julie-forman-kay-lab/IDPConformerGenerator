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

_name = 'ssext'
_help = 'Extract secondary structure profile.'

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
    '-r',
    '--reduced',
    help=(
        'Reduces nomenclature for secondary structure identity '
        'to \'L\', \'H\' and \'E\'.'
        ),
    action='store_true',
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
    """
    Execute main client function.

    Reads command line arguments and executes logic.
    """
    cmd = _load_args()
    main(**vars(cmd))


def main(
        ss_cmd,
        pdbs,
        output=None,
        ncores=1,
        reduced=False,
        **kwargs,
        ):
    """
    Run main cli logic.

    Parameters
    ----------
    ss_cmd : str or Path
        The command to run with subprocess module.

    pdbs : list
        A list of paths to PDB files or PDB file lists.

    output : string or Path, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.

    ncores : int
        The numbers of cores to use.
    """
    log.info(T('Extracting Secondary structure information'))
    init_files(log, LOGFILESNAME)

    log.info(T('reading input paths'))
    pdbs = libio.read_path_bundle(pdbs, ext='pdb')
    log.info(S('done'))

    log.info(T('preparing task execution'))
    log.info(S('for {} cores', ncores))

    manager = Manager()
    mdict = manager.dict()



    ss_ext_exec = libmulticore.JoinedResults(
        pdbs,
        ss_cmd,
        ncores=ncores,
        TaskMethod=libmulticore.DSSPTask,
        results_parser=libparse.DSSPParser.from_data_id_tuple,
        reduced=reduced,
        )
    log.info(S('executing...'))
    ss_ext_exec.run()

    log.info(T('exporting output to: {}', output))
    libparse.export_ss_from_DSSP(*ss_ext_exec.results, output=output)

    log.info(S('All done. Thanks!'))
    return


if __name__ == '__main__':
    maincli()
