"""
Extracts secondary structure information from PDBs.

Uses an external third party software.

USAGE:
    $ idpconfgen sscalc [PDBS]
"""
import argparse
import os
import subprocess
import pickle
import shutil
from multiprocessing import Manager, Pool


from idpconfgen import Path, log
from idpconfgen.libs import libcli, libio, libmulticore, libparse, libpdb
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = '.idpconfgen_sscalc'

_name = 'sscalc'
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
    '--complete',
    help='A previous DSSP DB file to complete with new entries.',
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
        ss_cmd,
        pdbs,
        output=None,
        ncores=1,
        reduced=False,
        complete=None,
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

    tmpdir = '__tmpsscalc__'

    log.info(T('reading input paths'))
    if pdbs[0].endswith('.tar'):
        pdbs2operate = libio.extract_from_tar(pdbs[0], output=tmpdir)
    else:
        pdbs2operate = list(libio.read_path_bundle(pdbs, ext='pdb'))
    log.info(S('done'))

    if complete:
        log.info(T(f'reading previous DSSP file: {complete}'))
        prev_dssp = libparse.read_stored_dict(complete)
        #prev_dssp = libparse.read_pipe_file(Path(complete).read_text())
        #pdbs2operate = list(filter(
        #    lambda x: x.stem not in prev_dssp.keys(),
        #    pdbs))
    else:
        prev_dssp = {}

    log.info(T('preparing task execution'))
    log.info(S('for {} cores', ncores))

    manager = Manager()
    mdict = manager.dict()

    try:
        libmulticore.pool_function(
            mkdssp,
            pdbs2operate,
            ss_cmd=ss_cmd,
            dssp_dict=mdict,
            reduced=reduced,
            ncores=ncores,
            )
    except Exception:
        log.error('FAILED')
    else:
        prev_dssp.update(mdict)
        libio.save_dictionary(prev_dssp, output=output)
        log.info(S('All done. Thanks!'))

    finally:
        shutil.rmtree(tmpdir)

    return



def mkdssp(pdb, ss_cmd, dssp_dict=None, reduced=False):

    cmd = [ss_cmd, '-i', os.fspath(pdb.resolve())]
    result = subprocess.run(cmd, capture_output=True)

    dssp, fasta, residues = libparse.parse_dssp(
        result.stdout.decode('utf-8'),
        reduced=reduced)
    dssp_dict[str(libpdb.PDBIDFactory(pdb))] = {
        'dssp': dssp,
        'fasta': fasta,
        'resids': residues,
        }


    #dssp_parser = libparse.DSSPParser(
    #    data=result.stdout.decode('utf8'),
    #    reduced=reduced,
    #    )

    #try:
    #    dssp_dict[str(libpdb.PDBIDFactory(pdb))] = {
    #        'dssp': ''.join(dssp_parser.ss),
    #        'fasta': ''.join(dssp_parser.fasta),
    #        'resids': ','.join(dssp_parser.resseq),
    #        }
    #except Exception:
    #    log.error(f'Error while saving to dict: {pdb}')

    return

if __name__ == '__main__':
    maincli()
