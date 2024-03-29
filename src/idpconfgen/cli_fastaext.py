"""
Extracts FASTA information from PDBs.

USAGE:
    $ idpconfgen fastaext [PDBS]
"""
import argparse
from functools import partial

from idpconfgen import Path, log
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import FileReaderIterator, save_dictionary
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libpdb import get_fasta_from_PDB
from idpconfgen.logger import T, init_files, report_on_crash


LOGFILESNAME = '.idpconfgen_fastaext'

_name = 'fastaext'
_help = 'Extract FASTA sequence from PDBs.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

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
    const='fasta.json',
    action=libcli.CheckExt({'.json'}),
    )

libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        output=None,
        ncores=1,
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
    """
    log.info(T('Extracting FASTA sequence information'))
    init_files(log, LOGFILESNAME)

    log.info(T('reading input paths'))
    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')

    execute = partial(
        report_on_crash,
        get_fasta_from_PDB,
        ROC_exception=IDPConfGenException,
        ROC_prefix=_name,
        )

    execute_pool = pool_function(execute, pdbs2operate, ncores=ncores)

    mdict = {i: j for i, j in execute_pool}

    save_dictionary(mdict, output)
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
