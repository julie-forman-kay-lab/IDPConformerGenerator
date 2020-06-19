"""
PDB/mmCIF Downloader.

Downloads structures from RCSB Databank for PDB formatted files of
individual chains.

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
* selects only the first model for multi MODEL structures
* renumbers atoms for saved chains
* passes through pdb-tools `pdb_delinsert` filter.

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

from idpconfgen.libs import libcli
from idpconfgen.libs.libdownload import download_structure
from idpconfgen.libs.libhigherlevel import download_pipeline


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


def main(*args, func=None, **kwargs):
    """Perform main logic."""
    f = download_pipeline(download_structure, LOGFILESNAME)
    f(*args, **kwargs)


if __name__ == '__main__':
    libcli.maincli(ap, main)
