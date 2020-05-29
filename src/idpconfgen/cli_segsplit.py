"""
Split PDB into segments

USAGE:
    do something
"""
import argparse
import traceback
from multiprocessing import Manager
import pickle

from idpconfgen.libs import libcli
from idpconfgen.libs.libio import read_path_bundle
from idpconfgen.logger import init_files, S, T
from idpconfgen.libs.libparse import read_pipe_file, group_consecutive_ints, identify_backbone_gaps, get_slice
from idpconfgen import log, Path
from idpconfgen.libs.libstructure import Structure, structure_to_pdb, write_PDB, col_resSeq, col_record
from idpconfgen.libs.libtimer import ProgressBar
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libio import make_destination_folder, glob_folder
from idpconfgen.libs import libmulticore

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

    #already_done = set(f'{i.name}_{i.chain}' for i in PDBList(libio.glob_folder('*.pdb')))

    #pdb_filtered = [pp for pp in pdb_paths if pp.stem in already_done]

    dssp_data = read_pipe_file(Path(dssp).read_text())

    manager = Manager()
    dssp_output_dict = manager.dict()

    libmulticore.pool_function(
        split_segs,
        pdbs_paths,
        dssps=dssp_data,
        minimum=2,
        dssp_out=dssp_output_dict,
        destination=destination,
        ncores=ncores,
        )

    Path(destination, 'dssp_segs.dssp').write_text(
        '\n'.join(f'{k}|{v}' for k, v in dssp_output_dict.items())
        )


def split_segs(pdbdata, dssps, minimum=2, dssp_out=None, destination=''):
    s = Structure(pdbdata)
    s.build()

    print(len(dssps[pdbdata.stem]))

    residues = [int(i) for i in dict.fromkeys(s.filtered_atoms[:, col_resSeq])]
    print(residues)

    # returns slices
    segments = group_consecutive_ints(residues)

    # removes ligands
    above_2 = filter(
        lambda x: x.stop - x.start > minimum,
        segments,
        )

    seg_counter = 0
    for seg in above_2:
        s.add_filter(lambda x: int(x[col_resSeq]) in residues[seg])

        if set(s.filtered_atoms[:, col_record]) == {'HETATM'}:
            s.pop_last_filter()
            continue


        # split regions with missing backbone
        s.add_filter_backbone(minimal=True)

        try:
                #identify_backbone_gaps(s.filtered_atoms)
            backbone_segs_in_resSeq_sets = \
                list(filter(
                    lambda x: len(x) > minimum,
                    get_slice(s.filtered_atoms)
                    ))
        except Exception as err:
            log.error(traceback.format_exc())
            log.error(repr(err))
            log.error(f'error in {pdbdata}')
            seg_counter += 1
            continue
        finally:
            s.pop_last_filter()


        for resSeq_set in backbone_segs_in_resSeq_sets:

            print(resSeq_set)

            s.add_filter(lambda x: x[col_resSeq] in resSeq_set)

            pdbout = f'{pdbdata.stem}_seg{seg_counter}'
            fout_seg = Path(destination, pdbout).with_suffix('.pdb')

            try:
                s.write_PDB(fout_seg, renumber=True)
            except EXCPTS.EmptyFilterError as err:
                log.error(f'Empty filter for: {repr(err)}')
            except EXCPTS.IDPConfGenException as err:
                log.error(S('* Something went wrong with {}: {}', pdbdata, repr(err)))
                log.debug(traceback.format_exc())

            s.pop_last_filter()

            dssp_slicing = slice(
                residues.index(int(resSeq_set[0])),
                residues.index(int(resSeq_set[-1])) + 1,
                None,
                )

            print(dssp_slicing)

            print(len(resSeq_set))
            print(len(dssps[pdbdata.stem][dssp_slicing]))
            print(dssps[pdbdata.stem][dssp_slicing])

            if not len(resSeq_set) == len(dssps[pdbdata.stem][dssp_slicing]):
                log.error(traceback.format_exc())
                log.error(f'error in {pdbdata}')
                seg_counter += 1
                continue
            dssp_out[pdbout] = dssps[pdbdata.stem][dssp_slicing]
            seg_counter += 1

