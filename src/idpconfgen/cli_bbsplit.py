"""
Split PDB into segments.

USAGE:
    do something
"""
import argparse

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libhigherlevel import split_backbone_segments
from idpconfgen.libs.libio import read_dictionary_from_disk, FileReaderIterator
from idpconfgen.libs.libmulticore import pool_chunks_to_disk_and_data_at_the_end
from idpconfgen.logger import init_files, S, T


LOGFILESNAME = '.idpconfgen_segsplit'

_name = 'bbsplit'
_help = 'Splits PDBs and data into backbone segments.'
_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_parser_pdbs(ap)

ap.add_argument(
    '--dssp',
    help='The DSSP file as saved by IDPConfGen SSCALC CLI.',
    default=None,
    )

ap.add_argument(
    '--dssp-dest',
    help='The DSSP file after split.',
    default='dssp_split.json',
    )

libcli.add_parser_destination_folder(ap)
libcli.add_argument_ncores(ap)


def main(
        pdbs,
        destination,
        chunks=5000,
        dssp=None,
        dssp_dest=None,
        ncores=1,
        **kwargs,
        ):
    """Perform main logic of the file."""
    init_files(log, LOGFILESNAME)

    pdbs2operate = FileReaderIterator(pdbs, ext='.pdb')

    sscalc_data = read_dictionary_from_disk(dssp) if dssp else None

    pool_chunks_to_disk_and_data_at_the_end(
        split_backbone_segments,
        pdbs2operate,
        destination,
        sscalc_data=sscalc_data,
        mdata_dest=dssp_dest,
        minimum=2,
        ncores=ncores,
        chunks=chunks,
        )
