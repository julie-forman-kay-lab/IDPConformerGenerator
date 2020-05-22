"""
Conformer Generator PDB Downloader.

Downloads structures from RCSB Databank provided a list of PDB
identifiers.

The PDB ID list can be given in the format of a file listing the PDBIDs
or as a list of arguments passed to the script call.

The following PDBID formats are allowed:

    - XXXX
    - XXXXY
    - XXXX_Y

where, XXXX is the PDB ID code and Y the chain identifier. Y can have
more then one character, for example, XXXXAB, will select chain 'AB'
of PDB ID XXXX. Digits are also allowed.

USAGE:
    $ icgpdbdl XXXX
    $ icgpdbdl XXXXY -d raw_pdbs
    $ icgpdbdl pdb.list -d raw_pdbs -u

"""
import argparse

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libparse import filter_pdb_for_db
from idpconfgen.libs.libdownload import download_structure
from idpconfgen.libs.libio import (
    read_path_bundle,
    concatenate_entries,
    glob_folder,
    make_destination_folder,
    )
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libpdb import PDBList
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = '.filter'

_name = 'pdb_filter'
_help = 'Filters PDBs for the database according to criteria.'
_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527



libcli.add_parser_pdbs(ap)
libcli.add_parser_destination_folder(ap)
libcli.add_argument_update(ap)
libcli.add_argument_ncores(ap)


def _load_args():
    cmd = ap.parse_args()
    return cmd


def main(
        pdbs,
        pdb_chains,
        destination=Path.cwd(),
        ncores=1,
        update=False,
        **kwargs
        ):
    """Run main script logic."""
    init_files(log, LOGFILESNAME)

    pdbs_paths = list(read_path_bundle(pdbs))
    pdbids_to_read = concatenate_entries(pdb_chains)
    pdb_chains = PDBList(pdbids_to_read)

    log.info(T('reading input PDB list'))
    log.info(S(f'from: {pdbs}'))
    #log.info(S(f'{str(pdblist)}'))
    log.info(S('done\n'))

    pdblist_destination = \
        PDBList(glob_folder(destination, '*.pdb'))
    log.info(T('reading destination folder'))
    log.info(S(f'from: {destination}'))
    log.info(S(f'{str(pdblist_destination)}'))
    log.info(S('done\n'))

    pdblist_comparison = pdb_chains.difference(pdblist_destination)
    log.info(T('Comparison between input and destination'))
    log.info(S(f'{str(pdblist_comparison)}'))
    log.info(S('done\n'))

    if update:

        dest = make_destination_folder(destination)

        pool_function(
            filter_pdb_for_db,
            pdbs_paths,
            ncores=ncores,
            # other kwargs for target function
            folder=dest,
            )

    #pdblist_updated = \
    #    PDBList(glob_folder(destination, '*.pdb'))
    #log.info(T('Reading UPDATED destination'))
    #log.info(S(f'{str(pdblist_updated)}'))
    #log.info(S('done\n'))

    return


def maincli():
    """Command-line interface entry point."""
    cmd = _load_args()
    main(**vars(cmd))


if __name__ == '__main__':
    maincli()
