"""
Fetchs structures from RCSB.

Structure files are saved to the local dist in its raw format,
that is, without any filtering or parsing.

Downloads structures from RCSB Databank provided a list of PDB
identifiers. The PDB ID list can be given in the format of a file
listing the PDBIDs or as a list of arguments passed to the script call.

The following PDBID formats are allowed:

    - XXXX

where, XXXX is the PDB ID code.

USAGE:
    $ idpconfgen fetch XXXX
    $ idpconfgen fetch XXXX -d destination_folder
    $ idpconfgen fetch pdb.list -d destination_folder -u
    $ idpconfgen fetch pdb.list -d file.tar -u -n
"""
import argparse

from idpconfgen.libs import libcli
from idpconfgen.libs.libdownload import fetch_raw_PDBs
from idpconfgen.libs.libhigherlevel import download_pipeline


LOGFILESNAME = '.raw_pdbdl'

_name = 'fetch'
_help = 'Fetch structures from RCSB.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )


libcli.add_argument_pdbids(ap)
libcli.add_argument_destination_folder(ap)
libcli.add_argument_update(ap)
libcli.add_argument_ncores(ap)
libcli.add_argument_chunks(ap)


# the func=None receives the `func` attribute from the main CLI interface
# defined at cli.py
def main(*args, func=None, **kwargs):
    """Perform main logic."""
    f = download_pipeline(fetch_raw_PDBs)
    f(*args, **kwargs)


if __name__ == '__main__':
    libcli.maincli(ap, main)
