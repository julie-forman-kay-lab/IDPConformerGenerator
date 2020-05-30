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
import os
from multiprocessing import Manager
from io import StringIO, BytesIO
from os import SEEK_END
import tarfile

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libdownload import download_structure
from idpconfgen.libs.libio import (
    concatenate_entries,
    glob_folder,
    make_destination_folder,
    )
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libpdb import PDBList
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = '.pdbdownloader'

_name = 'pdbdl'
_help = 'Downloads filtered structures from RCSB.'
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
    help='The coordinate PDB record name. Default: ("ATOM", "HETATM")',
    default=('ATOM', 'HETATM'),
    action=libcli.ArgsToTuple,
    nargs='+',
    )

ap.add_argument(
    '-raw',
    help='Only parse chain ID',
    action='store_true',
    )


libcli.add_argument_replace(ap)
libcli.add_argument_ncores(ap)


def _load_args():
    cmd = ap.parse_args()
    return cmd


def main(
        pdblist,
        destination=None,
        update=False,
        record_name=('ATOM',),
        ncores=1,
        replace=False,
        raw=False,
        **kwargs
        ):
    """Run main script logic."""
    init_files(log, LOGFILESNAME)

    pdbids_to_read = concatenate_entries(pdblist)
    pdblist = PDBList(pdbids_to_read)

    log.info(T('reading input PDB list'))
    log.info(S(f'from: {pdblist}'))
    log.info(S(f'{str(pdblist)}'))
    log.info(S('done\n'))

    if destination:
        pdblist_destination = \
            PDBList(glob_folder(destination, '*.pdb'))
        log.info(T('reading destination folder'))
        log.info(S(f'from: {destination}'))
        log.info(S(f'{str(pdblist_destination)}'))
        log.info(S('done\n'))

        pdblist_comparison = pdblist.difference(pdblist_destination)
        log.info(T('Comparison between input and destination'))
        log.info(S(f'{str(pdblist_comparison)}'))
        log.info(S('done\n'))

    if update:

        dest = make_destination_folder(destination)

        if not replace:
            pdb2dl = pdblist_comparison
        else:
            pdb2dl = pdblist

        tar = tarfile.open(os.fspath(Path(dest, 'pdbdl.tar.gz')), mode='w:gz', compresslevel=9)
        tar.close()

        chunk = 10_000
        tasks = list(pdb2dl.name_chains_dict.items())
        for i in range(0, len(tasks), chunk):
            task = tasks[i: i + chunk]

            manager = Manager()
            mlist = manager.list()

            pool_function(
                download_structure,
                task,
                ncores=ncores,
                # other kwargs for target function
                folder=dest,
                record_name=record_name,
                renumber=True,
                raw=raw,
                mlist=mlist,
                )

            tar = tarfile.open(os.fspath(Path(dest, 'pdbdl.tar.gz')), mode='a:')

            for fout, _data in mlist:
                try:
                    sIO = BytesIO()
                    sIO.write('\n'.join(_data).encode())
                    info = tarfile.TarInfo(name=fout)
                    info.size=sIO.seek(0, SEEK_END)
                    sIO.seek(0)
                    tar.addfile(tarinfo=info, fileobj=sIO)
                except Exception:
                    log.error(f'failed for {fout}')
            tar.close()

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
