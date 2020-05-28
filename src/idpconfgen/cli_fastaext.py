"""
Extracts FASTA information from PDBs.

USAGE:
    $ idpconfgen fastaext [PDBS]
"""
import argparse
import itertools as it
from multiprocessing import Manager

from idpconfgen import Path, log
from idpconfgen.libs import libcli, libio, libpdb, libstructure, libmulticore
from idpconfgen.logger import S, T, init_files


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
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    'pdbs',
    help='PDB file list.',
    nargs='+',
    )

ap.add_argument(
    '-o',
    '--output',
    help=( "The output file containing the PDBID and "
        "respective FASTA sequence information. "
        "Defaults to sys.stdout (prints to console)."
        ),
    type=Path,
    default=None,
    )

libcli.add_argument_ncores(ap)


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
    # tee is used to keep memory footprint low
    # though it would be faster to create a list from path_bundle
    #_, pdbs = it.tee(libio.read_path_bundle(pdbs, ext='.pdb'))
    pdbs_ = sorted(libio.read_path_bundle(pdbs, ext='.pdb'), key=lambda x: x.stem)
    #log.info(S('done'))
    #pdbids = libpdb.PDBList(pdbs_)
    #print(len(pdbids))
    #print(len(pdbs_))

    manager = Manager()
    mdict = manager.dict()

    libmulticore.pool_function(
        get_fastas,
        pdbs_,
        mdict=mdict,
        ncores=ncores,
        )

    libio.write_text('\n'.join(f'{k}|{v}' for k,v in mdict.items()), output)


def get_fastas(pdbfile, mdict=None):
    structure = libstructure.Structure(pdbfile)
    structure.build()
    fasta = structure.fasta
    assert len(fasta) == 1
    mdict[libpdb.PDBIDFactory(pdbfile)] = next(iter(fasta.values()))
    #out_data.append('{}|{}'.format(pdbid, next(iter(fasta.values()))))

    return


if __name__ == '__main__':
    maincli()
