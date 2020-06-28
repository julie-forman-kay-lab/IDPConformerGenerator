"""
Extracts secondary structure information from PDBs.

Uses an external third party software.

USAGE:
    $ idpconfgen sscalc [PDBS]
"""
import argparse
import shutil
import traceback
from functools import partial

from idpconfgen import Path, log
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_dictionary_from_disk,
    read_path_bundle,
    save_dictionary,
    save_pairs_to_disk,
    )
from idpconfgen.libs.libmulticore import (
    consume_iterable_in_list,
    pool_function_in_chunks,
    )
from idpconfgen.libs.libparse import mkdssp_w_split
from idpconfgen.logger import S, T, init_files, report_on_crash


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

libcli.add_argument_cmd(ap)
libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-o',
    '--output',
    help=(
        "A path to the file where the PDBID secondary structure and FASTA"
        " information dictionary will be saved. "
        "Defaults to sscalc.json, requires \'.json\' extension."
        ),
    type=Path,
    default='sscalc.json',
    action=libcli.CheckExt({'.json'}),
    )

# can't use the libcli option because this one is different
ap.add_argument(
    '-d',
    '--destination',
    help=(
        'Destination folder where the split PDB files will be stored. '
        'Defaults to current working directory.'
        'Alternatively, you can provide a path to a .tar file '
        'where PDBs will be saved.'
        ),
    type=Path,
    default='sscalc_splitted.tar',
    )

ap.add_argument(
    '-u',
    '--update',
    help='Updates a previous SSCALC JSON file.',
    type=Path,
    default=None,
    )

ap.add_argument(
    '-tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )

libcli.add_argument_reduced(ap)
libcli.add_argument_chunks(ap)
libcli.add_argument_ncores(ap)


def main(
        cmd,
        pdb_files,
        chunks=1000,
        destination='sscalc_splitted.tar',
        func=None,  # here just to receive from main cli.py
        minimum=2,
        ncores=1,
        output='sscalc_output.json',
        reduced=False,
        tmpdir=TMPDIR,
        update=None,
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

    chunks : int, optional
        The number of items to process in memory before saving the
        results from the disk.

    reduced : Bool, optional
        Whether to reduce secondary structure information to H/E/L

    ncores : int
        The numbers of cores to use.
    """
    log.info(T('Extracting Secondary structure information'))
    init_files(log, LOGFILESNAME)

    # update is inspected here because it can raise errors and it is better
    # that such errors are spotted before performing the expensive calculations
    if update:
        previous = read_dictionary_from_disk(update)

    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except TypeError:
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))

    log.info(T('preparing task execution'))

    try:
        consume_func = partial(
            consume_iterable_in_list,
            mkdssp_w_split,
            cmd=cmd,
            reduced=reduced,
            minimum=minimum,
            )

        execute = partial(
            report_on_crash,
            consume_func,
            ROC_exception=IDPConfGenException,
            ROC_prefix=_name,
            )

        # generator
        execute_pool = pool_function_in_chunks(
            execute,
            pdbs2operate,              # items to process
            ncores=ncores,
            chunks=chunks,
            )

        # this implementation is very specific for this case
        # this is why I haven spread it into functions for now
        dssp_data = {}  # stores DSSP data to save at the end
        pdb_data = {}  # temporarily stores data to be writen to the disk
        for chunk in execute_pool:
            for result in chunk:
                for fname, dsspdict, pdb_split in result:
                    dssp_data[fname] = dsspdict

                    # notice the copy, this is needed for the .clear()
                    # to work later on
                    pdb_data[f'{fname}.pdb'] = pdb_split
            save_pairs_to_disk(pdb_data.items(), destination=destination)
            pdb_data.clear()  # clears the dictionary to release memory

    except Exception as err:
        log.error('FAILED')
        log.debug(traceback.format_exc())
        raise err

    finally:
        if _istarfile:
            shutil.rmtree(tmpdir)

    if update:
        save_dictionary(previous.update(dssp_data), output)
    else:
        save_dictionary(dssp_data, output)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
