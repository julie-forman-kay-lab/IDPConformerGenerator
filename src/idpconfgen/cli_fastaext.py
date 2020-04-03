"""
Extracts FASTA information from PDBs.

USAGE:
    $ idpconfgen fastaext [PDBS]
"""
import argparse

from idpconfgen import Path, log
from idpconfgen.libs import libcli, libio, libmulticore, libparse
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = 'idpconfgen_fastaext'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

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
    '-o',
    '--output',
    help=(
        "The output file containing the PDBID and "
        "respective FASTA sequence information. "
        "Defaults to sys.stdout (prints to console)."
        ),
    type=Path,
    default=None,
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
        pdbs,
        output=None,
        **kwargs,
        ):
    """
    Run main cli logic.

    Parameters
    ----------
    pdbs : list
        A list of paths to PDB files or PDB file lists.

    output : string or Path, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.

    ncores : int
        The numbers of cores to use.
    """
    log.info(T('Extracting FASTA sequence information'))
    init_files(log, LOGFILESNAME)
     
    log.info(T('reading input paths'))
    pdbs = libio.read_path_bundle(pdbs)
    log.info(S('done'))
    pdbids = libpdb.PDBList(pdbs)
     
    out_data = []
    for pdbid, pdbfile in zip(pdbids, pdbs):
        structure = PDBDataFactory(Path(fin).read_text())
        structure.build()
        out_data.append('{}|{}'.format((pdbid, structure.fasta)))
     
    libio.write_text('\n'.join(out_data), output)
    return


if __name__ == '__main__':
    maincli()
