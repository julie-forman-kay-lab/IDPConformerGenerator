import argparse

from idpconfgen import Path, log
from idpconfgen.libs import libcli, libio
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = 'idpconfgen_ssext'

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    'pdbs',
    help='PDB file list.',
    nargs='+',
    )

ap.add_argument(
    '-c',
    '--ss_cmd',
    help='The path to the DSSP executable file.',
    default=None,
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
        pdbs,
        ss_cmd,
        output=None,
        ncores=1,
        **kwargs,
        ):
    """Run main cli logic."""
    
    init_files()
    
    pdbs_paths = libio.read_paths(pdbs)

    ss_ext_exec = 

    return


if __name__ == '__main__':

    maincli()

