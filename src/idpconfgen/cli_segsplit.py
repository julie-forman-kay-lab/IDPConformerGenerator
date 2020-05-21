"""
Split PDB into segments

USAGE:
    do something
"""
import argparse
import traceback
from multiprocessing import Manager

from idpconfgen.libs import libcli
from idpconfgen.libs.libio import read_path_bundle
from idpconfgen.logger import init_files, S, T
from idpconfgen.libs.libparse import read_pipe_file, group_consecutive_ints
from idpconfgen import log, Path
from idpconfgen.libs.libstructure import Structure, structure_to_pdb, write_PDB, col_resSeq
from idpconfgen.libs.libtimer import ProgressBar
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libio import make_destination_folder

LOGFILESNAME = '.pdb_segsplit'

_name = 'pdb_segsplit'
_help = 'Splits PDB into segments'
_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    'dssp',
    help='The DSSP file as saved by IDPConfGen SSEXT CLI',
    )

libcli.add_parser_pdbs(ap)
libcli.add_parser_destination_folder(ap)
libcli.add_argument_ncores(ap)


def main(pdbs, dssp, destination=None , ncores=1, **kwargs):
    destination = make_destination_folder(destination)

    init_files(log, LOGFILESNAME)

    pdbs_paths = list(read_path_bundle(pdbs))

    dssp_data = read_pipe_file(Path(dssp).read_text())

    dssp_output_dict = {}
    with ProgressBar(len(pdbs_paths)) as pb:
        for pdbid in pdbs_paths:

            s = Structure(pdbid)
            s.build()

            residues = list(set(int(i) for i in s.filtered_atoms[:, col_resSeq]))

            segments = group_consecutive_ints(residues)

            above_2 = filter(
                lambda x: x.stop - x.start > 2,
                segments,
                )

            for i, seg in enumerate(above_2):
                pdbout = f'{pdbid.stem}_seg{i}'
                fout_seg = Path(destination, pdbout).with_suffix('.pdb')


                s.add_filter(lambda x: x[col_resSeq] in residues[seg])

                try:
                    s.write_PDB(fout_seg)
                except EXCPTS.IDPConfGenException as err:
                    log.error(S('* Something went wrong with {}: {}', pdbid, repr(err)))
                    log.debug(traceback.format_exc())

                s.pop_last_filter()

                dssp_output_dict[pdbout] = dssp_data[pdbid.stem][seg]
            pb.increment()

    Path(destination, 'dssp_segs.dssp').write_text(
        '\n'.join(f'{k}|{v}' for k, v in dssp_output_dict.items())
        )

