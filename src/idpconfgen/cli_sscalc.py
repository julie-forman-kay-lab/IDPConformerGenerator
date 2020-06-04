"""
Extracts secondary structure information from PDBs.

Uses an external third party software.

USAGE:
    $ idpconfgen sscalc [PDBS]
"""
import argparse
import shutil
from multiprocessing import Manager

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_path_bundle,
    read_dictionary_from_disk,
    save_dictionary,
    )
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libparse import mkdssp
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = '.idpconfgen_sscalc'
TMPDIR = '__tmpsscalc__'

_name = 'sscalc'
_help = 'Calculate secondary structure profile.'

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

libcli.add_parser_pdbs(ap)

ap.add_argument(
    '-o',
    '--output',
    help=(
        "A path to a new file where the PDBID seconda structure and FASTA"
        " information dictionary will be saved. "
        "Defaults to sscalc.json."
        ),
    type=Path,
    default='sscalc.json',
    action=libcli.CheckExt({'.json'}),
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

    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdbs, output=TMPDIR)
        _istarfile = True
    except TypeError:
        pdbs2operate = list(read_path_bundle(pdbs, ext='pdb'))
        _istarfile = False
    log.info(S('done'))

    if complete:
        log.info(T(f'reading previous DSSP file: {complete}'))
        prev_dssp = read_dictionary_from_disk(complete)
    else:
        prev_dssp = {}

    log.info(T('preparing task execution'))

    manager = Manager()
    mdict = manager.dict()

    try:
        pool_function(
            mkdssp,
            pdbs2operate,
            ncores=ncores,
            # kwargs for mkdssp function
            ss_cmd=ss_cmd,
            mdict=mdict,
            reduced=reduced,
            )
    except Exception as err:
        log.error('FAILED')
        raise err
    else:
        prev_dssp.update(mdict)
        save_dictionary(dict(sorted(prev_dssp.items())), output=output)
        log.info(S('All done. Thanks!'))
    finally:
        if _istarfile:
            shutil.rmtree(TMPDIR)

    return


if __name__ == '__main__':
    maincli()
