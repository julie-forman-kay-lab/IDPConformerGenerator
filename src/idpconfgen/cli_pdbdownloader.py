"""
PDBDOWNLOADER

Downloads structures from RCSB Databank.

USAGE:
    >>> icgpdbdl CULL

"""
#import argparse
import sys

from idpconfgen import CustomParser, ArgsToTuple, log, Path
from idpconfgen.logger import init_files, T, S
from idpconfgen.libs import libcli, libpdb, libsearch


LOGFILESNAME = 'pdbdownloader'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = CustomParser(
    prog=_prog,
    description=_des,
    usage=_us,
    )
# https://stackoverflow.com/questions/24180527 

ap.add_argument(
    'pdblist',
    help='A list of PDBID:CHAIN to download.',
    type=Path,
    )

ap.add_argument(
    '-d',
    '--destination',
    help=(
        'Destination folder where PDB files will be stored.'
        'Defaults to current working directory.'
        ),
    type=Path,
    default=Path.cwd(),
    )

ap.add_argument(
    '-u',
    '--update',
    help='Updates destination folder according to input PDB list.',
    action='store_true',
    )

ap.add_argument(
    '-rn',
    '--record_name',
    help='The coordinate PDB record name.',
    default=('ATOM',),
    action=ArgsToTuple, 
    nargs='+',
    )

ap.add_argument(
    '-n',
    '--ncores',
    help='Number of cores to use.',
    type=int,
    default=1,
    )


def load_args():
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
    
    init_files(log, LOGFILESNAME)
    
    with pdblist.open('r') as fh:
        pdblist = libpdb.PDBList(fh.readlines())
    
    log.info(T('reading input PDB list'))
    log.info(S(f'from: {pdblist}'))
    log.info(S(f'{str(pdblist)}'))
    log.info(S('done\n'))
    
    if destination:
        pdblist_destination = \
            libpdb.PDBList(libsearch.glob_folder(destination, '*.pdb'))
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
        pdblist_updated = libpdb.PDBList(libsearch.glob_folder(destination, '*.pdb'))
        log.info(T('Reading UPDATED destination'))
        log.info(S(f'{str(pdblist_updated)}'))
        log.info(S('done\n'))
    
    return


def maincli():
    cmd = load_args()
    main(**vars(cmd))


if __name__ == '__main__':
    maincli()
