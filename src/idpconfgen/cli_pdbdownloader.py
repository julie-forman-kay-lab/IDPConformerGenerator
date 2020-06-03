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
from idpconfgen.libs.libdownload import get_pdbs_downloader
from idpconfgen.libs.libio import concatenate_entries, read_PDBID_from_source
from idpconfgen.libs.libpdb import PDBList
from idpconfgen.logger import S, T, init_files


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
    'pdbids',
    help='PDBID:CHAIN identifiers to download.',
    nargs='+',
    )

ap.add_argument(
    '-d',
    '--destination',
    help=(
        'Destination folder where PDB files will be stored. '
        'Defaults to current working directory.'
        'Alternatively, you can provide a path to a .tar file '
        'where PDBs will be saved.'
        ),
    type=Path,
    default=Path.cwd(),
    )

ap.add_argument(
    '-u',
    '--update',
    help=(
        'Updates destination folder/file according to input PDB list. '
        'If not provided, a comparison between input and destination '
        'is made.'
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
    '-c',
    '--chunks',
    help='Number of chunks to process in memory before saving to disk.',
    default=5_000,
    type=int,
    )

libcli.add_argument_ncores(ap)


def _load_args():
    cmd = ap.parse_args()
    return cmd


def main(
        pdbids,
        chunks=5_000,
        destination=None,
        ncores=1,
        record_name=('ATOM', 'HETATM'),
        update=False,
        **kwargs,
        ):
    """Run main script logic."""
    init_files(log, LOGFILESNAME)

    #
    log.info(T('reading input PDB list'))

    pdblist = PDBList(concatenate_entries(pdbids))

    log.info(
        f"{S(str(pdblist))}\n"
        f"{S('done')}\n"
        )

    #
    log.info(T('Comparison between input and destination'))
    pdblist_from_destination = read_PDBID_from_source(destination)
    pdblist_comparison = pdblist.difference(pdblist_from_destination)
    log.info(
        f"{S(f'{str(pdblist_comparison)}')}\n"
        f"{S('done')}\n"
        )

    if update:

        downloader = get_pdbs_downloader(destination)

        downloader(
            destination,
            sorted(list(pdb2dl.name_chains_dict.items()), key=lambda x: x[0]),
            ncores=ncores,
            chunks=chunks,
            record_name=record_name,
            )

        log.info(T('Reading UPDATED destination'))
        pdblist_updated = read_PDBID_from_source(destination)
        pdblist_up_comparison = pdblist.difference(pdblist_updated)
        log.info(S(f'{str(pdblist_up_comparison)}'))
        if len(pdblist_up_comparison) > 0:
            log.info(S('There are PDBIDs not downloaded'))
            log.info(S(
                'Those IDs have been registered in the '
                f'{LOGFILESNAME}.debug file.'
                ))
            log.debug(
                '\n'.join(str(_pdbid) for _pdbid in pdblist_up_comparison)
                )

    log.info(S('done\n'))
    return


def maincli():
    """Command-line interface entry point."""
    cmd = _load_args()
    main(**vars(cmd))


if __name__ == '__main__':
    maincli()
