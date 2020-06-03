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
from idpconfgen.libs.libio import read_PDBID_from_folder


LOGFILESNAME = '.idpconfgen_pdbdl'

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
        'Alternatively, you can provide a path to a .tar file'
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
    help='Download the raw PDBs/mmCIFs, ignoring the CHAIN IDS.',
    action='store_true',
    )

ap.add_argument(
    '-chunks',
    help='Number of chunks to process in memory before saving to disk.',
    default=10_000,
    type=int,
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
        record_name=('ATOM', 'HETATM'),
        ncores=1,
        replace=False,
        raw=False,
        chunks=10_000,
        **kwargs
        ):
    """Run main script logic."""
    init_files(log, LOGFILESNAME)

    #
    log.info(T('reading input PDB list'))
    log.info(S(f'from: {pdblist}'))

    pdbids_to_read = concatenate_entries(pdblist)
    pdblist = PDBList(pdbids_to_read)

    log.info(S(f'{str(pdblist)}'))
    log.info(S('done\n'))

    #
    pdblist_destination = list_destination(destination)
    pdblist_comparison = pdblist.difference(pdblist_destination)
    log.info(T('Comparison between input and destination'))
    log.info(S(f'{str(pdblist_comparison)}'))
    log.info(S('done\n'))

    if update:

        if replace:
            pdb2dl = pdblist
        else:
            pdb2dl = pdblist_comparison

        DLER = {
            True: PDBFileDownloader,
            destination.suffix == '.tar': PDBFileTarDownloader,
            }

        downloader = DLER[True](destination)

        downloader.download(
            pdb2dl,
            ncores=ncores,
            chunks=chunks,
            record_name=record_name,
            )

        #pdblist_updated = \
        #    PDBList(glob_folder(destination, '*.pdb'))
        #log.info(T('Reading UPDATED destination'))
        #log.info(S(f'{str(pdblist_updated)}'))
        list_destination(destination)

    log.info(S('done\n'))
    return


def list_destination(destination):
    if destination.suffix == '.tar':
        pdblist_destination = read_PDBID_from_tar(destination)
    elif destination.is_dir:
        pdblist_destination = read_PDBID_from_folder(destination)

    log.info(T(f'destination {destination}'))
    log.info(S(str(pdblist_destination)))
    return pdblist_destination


def read_PDBID_from_tar(destination):
    try:
        tar = tarfile.open(destination.str(), 'r:*')
        p = PDBList(tar.getnames())
        tar.close()
    except FileNotFoundError:
        p = PDBList([])
    return p


class PDBFileDownloader:
    def __init__(self, destination):
        self.dest = make_destination_folder(destination)
        return

    def download(self, pdbs2dl, **kwargs):
        for result_dict in download_script(pdbs2dl, **kwargs):
            for fname, data in sorted(result_dict.items()):
                with open(Path(self.dest, fname), 'w') as fout:
                    fout.write('\n'.join(data))



class PDBFileTarDownloader:
    _exists = {True: 'a', False: 'w'}
    def __init__(self, destination):
        self.dests = destination.str()
        dest = tarfile.open(
            self.dests,
            mode=self._exists[destination.exists()],
            )
        dest.close()
        return


    def download(self, pdbs2dl, **kwargs):
        for result_dict in download_script(pdbs2dl, **kwargs):

            tar = tarfile.open(self.dests, mode='a:')

            for fout, _data in sorted(result_dict.items()):
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



def download_script(
        pdbids2dl,
        ncores=1,
        chunks=10_000,
        record_name=None,
        ):
    """
    """

    tasks = sorted(list(pdbids2dl.name_chains_dict.items()), key=lambda x: x[0])
    for i in range(0, len(tasks), chunks):
        task = tasks[i: i + chunks]

        manager = Manager()
        mdict = manager.dict()

        pool_function(
            download_structure,
            task,
            ncores=ncores,
            # other kwargs for target function
            record_name=record_name,
            mdict=mdict,
            )

        yield mdict


def maincli():
    """Command-line interface entry point."""
    cmd = _load_args()
    main(**vars(cmd))


if __name__ == '__main__':
    maincli()
