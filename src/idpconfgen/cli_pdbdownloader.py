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
of PDB ID XXXX (for mmCIF cases); digits are also allowed. If no chainID
is provided, saves each chain of the PDB file separately.

Detailed procedures:
* PDBs/mmCIFs are saved parsed in PDB format.
* Known solvent and ligands are removed
* Considers only altLoc 'A' or ' '.
* Considers only elements composing aminoacids
* selects only the first model in multi MODEL structures
* renumbers atoms for saved chains
* parses output through pdb-tools `pdb_delinsert` filter.

Accepts TAR files as output destination.

DO NOT forget to use the `-u` parameter to perform the actual download.
Otherwise a simple comparison between source and destination is performed.

USAGE:
    $ idpconfgen pdbdl XXXX
    $ idpconfgen pdbdl XXXXY -d <FOLDER>
    $ idpconfgen pdbdl pdbid.list -d <FOLDER> -u
    $ idpconfgen pdbdl pdbid.list -d <DESTINATION TAR FILE> -u -d
"""
import argparse
from functools import reduce

from idpconfgen import Path, log
from idpconfgen.core.definitions import blocked_ids
from idpconfgen.libs import libcli
from idpconfgen.libs.libdownload import download_pdbs_from_ids
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

libcli.add_argument_pdbids(ap)
libcli.add_argument_destination_folder(ap)
libcli.add_argument_update(ap)
libcli.add_argument_ncores(ap)
libcli.add_argument_chunks(ap)


def _load_args():
    """Load arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    """Command-line interface entry point."""
    cmd = _load_args()
    main(**vars(cmd))


def main(
        pdbids,
        chunks=5_000,
        destination=None,
        func=None,
        ncores=1,
        record_name=('ATOM', 'HETATM'),
        update=False,
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
    log.info(T('Filtering input'))
    destination = destination or Path.cwd()
    log.info(
        f"{S(f'from destination: {destination}')}\n"
        f"{S('and other sources...')}"
        )

    # comparison block
    def diff(first, other):
        return first.difference(other)

    remove_from_input = [
        read_PDBID_from_source(destination),
        PDBList(blocked_ids),
        ]

    # yes, there are just two items in remove_from_input, why use reduce?
    # what if more are added in the future? :-P the engine is already created
    pdblist_comparison = reduce(diff, remove_from_input, pdblist)
    log.info(S(f'Found {str(pdblist_comparison)} to download'))
    #

    something_to_download = len(pdblist_comparison) > 0
    if something_to_download and update:

        download_pdbs_from_ids(
            destination,
            sorted(pdblist_comparison.name_chains_dict.items()),
            ncores=ncores,
            chunks=chunks,
            record_name=record_name,
            )

        log.info(T('Reading UPDATED destination'))
        pdblist_updated = read_PDBID_from_source(destination)
        pdblist_up_comparison = pdblist.difference(pdblist_updated)
        log.info(S(f'{str(pdblist_up_comparison)}'))
        if len(pdblist_up_comparison) > 0:
            log.info(S(
                'There are PDBIDs not downloaded\n.'
                'Those IDs have been registered in the '
                f'{LOGFILESNAME}.debug file.'
                ))
            log.debug('\n'.join(str(_id) for _id in pdblist_up_comparison))

    elif not something_to_download and update:
        log.info(S('There is nothing to download.'))
        log.info(S('All requested IDs are already at the destination folder.'))

    log.info(T('PDB Downloader finished'))
    return


if __name__ == '__main__':
    maincli()
