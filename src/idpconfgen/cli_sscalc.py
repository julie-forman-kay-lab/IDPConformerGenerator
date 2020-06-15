"""
Extracts secondary structure information from PDBs.

Uses an external third party software.

USAGE:
    $ idpconfgen sscalc [PDBS]
"""
import argparse
import shutil

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_path_bundle,
    read_dictionary_from_disk,
    )
from idpconfgen.libs.libmulticore import pool_chunks_to_disk_and_data_at_the_end
from idpconfgen.libs.libparse import mkdssp_w_split
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

libcli.add_argument_cmd(ap)
libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-o',
    '--output',
    help=(
        "A path to the file where the PDBID secondary structure and FASTA"
        " information dictionary will be saved. "
        "Defaults to sscalc.json."
        ),
    type=Path,
    default='sscalc.json',
    action=libcli.CheckExt({'.json'}),
    )

libcli.add_argument_reduced(ap)
libcli.add_argument_chunks(ap)
#libcli.add_argument_update(ap)
libcli.add_argument_ncores(ap)


def main(
        cmd,
        pdb_files,
        chunks=5000,
        func=None,
        ncores=1,
        output='sscalc_output.json',
        reduced=False,
        minimum=2,
        split_destination='sscalc_splitted.tar',
        #update=False,
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
        pdbs2operate = extract_from_tar(pdb_files, output=TMPDIR, ext='.pdb')
        _istarfile = True
    except TypeError:
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))

    log.info(T('preparing task execution'))

    try:
        pool_chunks_to_disk_and_data_at_the_end(
            mkdssp_w_split,
            pdbs2operate,
            destination=split_destination,
            ncores=ncores,
            chunks=chunks,
            #mdata_source=prev_dssp,
            mdata_dest=output,
            # kwargs for mkdssp function
            cmd=cmd,
            reduced=reduced,
            minimum=minimum,
            )
    except Exception as err:
        log.error('FAILED')
        raise err
    finally:
        if _istarfile:
            shutil.rmtree(TMPDIR)
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
