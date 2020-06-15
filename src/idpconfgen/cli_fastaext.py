"""
Extracts FASTA information from PDBs.

USAGE:
    $ idpconfgen fastaext [PDBS]
"""
import argparse
from multiprocessing import Manager

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.logger import S, T, init_files
from idpconfgen.libs.libhigherlevel import get_fastas
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libio import FileReaderIterator, save_dictionary


LOGFILESNAME = 'idpconfgen_fastaext'

_name = 'fastaext'
_help = 'Extract FASTA sequence from PDBs.'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
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
    #pdbs_ = sorted(libio.read_path_bundle(pdbs, ext='.pdb'), key=lambda x: x.stem)
    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')

    manager = Manager()
    mdict = manager.dict()

    pool_function(
        get_fastas,
        pdbs2operate,
        mdict=mdict,
        ncores=ncores,
        )

    #libio.write_text('\n'.join(f'{k}|{v}' for k,v in mdict.items()), output)
    save_dictionary(mdict, output)


if __name__ == '__main__':
    libcli.maincli(ap, main)
