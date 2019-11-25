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
from idpconfgen.libs import libcli, libfile, libio, libpdb
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = 'pdbdownloader'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    'pdblist',
    help='PDBID:CHAIN identifiers to download.',
    nargs='+',
    )

ap.add_argument(
    '-d',
    '--destination',
    help=(
        'Destination folder where PDB files will be stored.'
        ' Defaults to current working directory.'
        ),
    type=Path,
    default=Path.cwd(),
    )

ap.add_argument(
    '-u',
    '--update',
    help=(
        'Updates destination folder according to input PDB list. '
        'If not provided a comparison between input and destination '
        'is provided.'
        ),
    action='store_true',
    )

ap.add_argument(
    '-rn',
    '--record_name',
    help='The coordinate PDB record name. Default: ("ATOM",)',
    default=('ATOM',),
    action=libcli.ArgsToTuple,
    nargs='+',
    )

ap.add_argument(
    '-n',
    '--ncores',
    help='Number of cores to use.',
    type=int,
    default=1,
    )


def _load_args():
    cmd = ap.parse_args()
    return cmd


def main(
        pdblist,
        destination=None,
        update=False,
        record_name=('ATOM',),
        ncores=1,
        **kwargs
        ):
    """Run main script logic."""
    init_files(log, LOGFILESNAME)
    
    pdbids_to_read = libio.concatenate_entries(pdblist)
    pdblist = libpdb.PDBList(pdbids_to_read)

    log.info(T('reading input PDB list'))
    log.info(S(f'from: {pdblist}'))
    log.info(S(f'{str(pdblist)}'))
    log.info(S('done\n'))
    
    if destination:
        pdblist_destination = \
            libpdb.PDBList(libfile.glob_folder(destination, '*.pdb'))
        log.info(T('reading destination folder'))
        log.info(S(f'from: {destination}'))
        log.info(S(f'{str(pdblist_destination)}'))
        log.info(S('done\n'))
        
        pdblist_comparison = pdblist.difference(pdblist_destination)
        log.info(T('Comparison between input and destination'))
        log.info(S(f'{str(pdblist_comparison)}'))
        log.info(S('done\n'))
    
    if update:
        pdbdownloader = libpdb.PDBDownloader(
            pdblist_comparison,
            destination,
            record_name=record_name,
            ncores=ncores,
            )
        pdbdownloader.prepare_pdb_list()
        pdbdownloader.run()
        pdblist_updated = \
            libpdb.PDBList(libfile.glob_folder(destination, '*.pdb'))
        log.info(T('Reading UPDATED destination'))
        log.info(S(f'{str(pdblist_updated)}'))
        log.info(S('done\n'))
    
    return


def maincli():
    """Command-line interface entry point."""
    cmd = _load_args()
    main(**vars(cmd))


if __name__ == '__main__':
    maincli()
