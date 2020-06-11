"""
Split PDB into segments.

USAGE:
    do something
"""
import shutil
import argparse
import os
import time
from collections import defaultdict

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libhigherlevel import split_backbone_segments
from idpconfgen.libs.libio import read_dictionary_from_disk, FileReaderIterator, extract_from_tar, save_dictionary, save_pairs_to_disk
from idpconfgen.libs.libmulticore import pool_chunks_to_disk_and_data_at_the_end
from idpconfgen.logger import init_files, S, T
from idpconfgen.libs.libtimer import ProgressWatcher


TMPDIR = '__tmpsscalc__'
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

libcli.add_argument_pdb_files(ap)

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

ap.add_argument(
    '-c',
    '--chunks',
    help='Number of chunks to process in memory before saving to disk.',
    default=5_000,
    type=int,
    )

ap.add_argument(
    '-m',
    '--minimum',
    help='The minimum length of the protein chunk.',
    default=2,
    )

libcli.add_argument_destination_folder(ap)
libcli.add_argument_ncores(ap)


def main(
        pdbs,
        destination,
        chunks=5000,
        dssp=None,
        dssp_dest=None,
        ncores=1,
        minimum=2,
        **kwargs,
        ):
    """Perform main logic of the file."""
    init_files(log, LOGFILESNAME)

    pdbs2operate = FileReaderIterator(pdbs, ext='.pdb')
    #extract_from_tar(pdbs, output=TMPDIR)
    #pdbs2operate = FileReaderIterator([TMPDIR], ext='.pdb')

    sscalc_data = read_dictionary_from_disk(dssp) if dssp else None

    mdata = defaultdict(dict)
    with ProgressWatcher(pdbs2operate) as PW:
        for i in range(0, len(pdbs2operate), chunks):
            mfiles = {}
            task = pdbs2operate[i: i + chunks]
            for t in task:
                split_backbone_segments(
                    t,
                    sscalc_data = sscalc_data,
                    mfiles=mfiles,
                    mdata=mdata,
                    minimum=minimum,
                    )
                PW.increment()
            save_pairs_to_disk(dict(sorted(mfiles.items())).items(), destination=destination)
    save_dictionary(dict(sorted(mdata.items())), dssp_dest)



    #pool_chunks_to_disk_and_data_at_the_end(
    #    split_backbone_segments,
    #    pdbs2operate,
    #    destination,
    #    sscalc_data=sscalc_data,
    #    mdata_dest=dssp_dest,
    #    minimum=2,
    #    ncores=ncores,
    #    chunks=chunks,
    #    )
