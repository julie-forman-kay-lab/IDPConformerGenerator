"""
Split PDB into segments

USAGE:
    do something
"""
import argparse
import traceback
from collections import defaultdict
from multiprocessing import Manager
import pickle

from idpconfgen.libs import libcli
from idpconfgen.libs.libio import read_path_bundle, read_dictionary_from_disk, FileReaderIterator
from idpconfgen.logger import init_files, S, T
from idpconfgen.libs.libparse import read_pipe_file, group_consecutive_ints, identify_backbone_gaps, get_segments_based_on_backbone_continuity, group_runs
from idpconfgen import log, Path
from idpconfgen.libs.libstructure import Structure, structure_to_pdb, write_PDB, col_resSeq, col_record
from idpconfgen.libs.libtimer import ProgressBar
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libio import make_destination_folder, glob_folder
from idpconfgen.libs.libmulticore import pool_chunks_to_disk_and_data_at_the_end

LOGFILESNAME = '.idpconfgen_segsplit'

_name = 'bbsplit'
_help = 'Splits data into backbone segments'
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
    help='The DSSP file as saved by IDPConfGen SSCALC CLI',
    default=None,
    )

libcli.add_parser_destination_folder(ap)
libcli.add_argument_ncores(ap)


def main(
        pdbs,
        dssp=None,
        destination=None,
        ncores=1,
        chunks=5000,
        **kwargs,
        ):
    """
    """

    init_files(log, LOGFILESNAME)

    pdbs2operate = FileReaderIterator(pdbs, ext='.pdb')

    # read dssp from disk (can be json or pickle)
    dssp_data = read_dictionary_from_disk(dssp) if dssp else None

    # needs to account for tar files
    #destination = make_destination_folder(destination)

    pool_chunks_to_disk_and_data_at_the_end(
        split_segments,
        pdbs2operate,
        destination,
        data=dssp_data,
        minimum=2,
        ncores=ncores,
        chunks=chunks,
        )


def split_segments(
        pdbname,
        pdbdata,
        mfiles,
        sscalc_data=None,
        mdict=None,
        minimum=2,
        ):
    """
    """
    # this might be slightly slower but it definitively more modular
    # and testable
    # `structure_segments` are in residue number (str)
    structure, structure_segments = backbone_split(pdbdata, minimum=2)

    splitted_segments = \
        split_structure_in_segments(structure, structure_segments)

    assert len(structure_segments) == len(splitted_segments), pdbname
    for i, txt in enumerate(splitted_segments):
        mfiles[f'{pdbname}_seg{i}'] = txt

    # here I have to put them in the mfiles
    try:
        dssp_segments = split_sscalc_data(sscalc_data, splitted_segments)
        assert len(dssp_segments) == len(structure_segments), pdbname
    except TypeError:  # dssp_data is None
        pass
    else:
        mdict.update(dssp_segments)


def split_sscalc_data(pdbname, data, segments):
    """
    """
    pdbdata = data[pdbname]
    structural_data = (k for k in pdbdata.keys() if k != 'resids')
    residues = structural_data['resids'].split(',')

    splitdata = defaultdict(dict)
    for i, segment in enumerate(segments):

        assert all(s.isdigit() for s in segment), pdbname

        key = f'{pdbname}_seg{i}'
        for datatype in structural_data:
            splitdata[key][datatype] = \
                ''.join(
                    c
                    for res, c in zip(residues, pdbdata[datatype])
                    if res in segment
                    )

        splitdata[key]['resids'] = ','.join(segment)

    return splitdata


def split_structure_in_segments(structure, residue_segments):
    """
    Parameters
    ----------
    structure : :class:`libs.libstructure.Structure`

    residue_seguments: list of str that are digits

    Yields
    ------
    str (PDB text data) for each segment
    """
    for segment in residue_segments:
        structure.clear_filters()
        structure.add_filter(lambda x: x[col_resSeq] in segment)
        txt = '\n'.join(structure.get_PDB())
        yield txt



def backbone_split(pdbdata, minimum=2):
    """
    Split PDBs in continuous backbone chunks.

    Parameters
    ----------
    pdbdata : str
        Text of the PDB file.

    minimum : int
        Minimum backbone to accept.
    """
    hetatm_set = {'HETATM'}

    s = Structure(pdbdata)
    s.build()

    # this will ignore the iCode
    # please use pdb-tool pdb_delinsertions before this step
    # PDBs downloaded with IDPConfGen already correct for these
    # see libs.libparse.delete_insertions
    residues = [int(i) for i in dict.fromkeys(s.filtered_atoms[:, col_resSeq])]

    consecutive_residue_segments = group_runs(residues)

    # removes segments shorter than the minimum allowed
    above_minimum = filter(
        lambda segment: len(segment) > minimum,
        consecutive_residue_segments,
        )

    #seg_counter = 0
    residue_segments = []
    for segment in above_minimum:
        s.clear_filters()
        s.add_filter(lambda x: int(x[col_resSeq]) in segment)

        # discard segments that might be composed only of HETATM
        tmp_filtered = s.filtered_atoms
        if set(tmp_filtered[:, col_record]) == hetatm_set:
        #    log.debug(
        #        'Found a segment with only HETATM.\n'
        #        'You may wish to record those HETATM residue codes:\n'
        #        f'{set(tmp_filtered[:, col_resName])}\n'
        #        )
            continue

        # identify backbone gaps
        s.add_filter_backbone(minimal=True)
        #try:
        residue_segments.extend(
            list(filter(
                lambda segment: len(segment) > minimum,
                get_segments_based_on_backbone_continuity(s.filtered_atoms),
                ))
            )
        #except Exception as err:
        #    log.error(traceback.format_ext())
        #    log.error(repr(err))
        #    log.error(f'error in {pdbname}')
        #    continue
        #finally:
        #s.pop_last_filter()

        # this loop will go at least once
        #for resSeq_set in consecutive_backbone_segs:
        #    s.add_filter(lambda x: x[col_resSeq] in resSeq_set)

        #    pdbout = f'{pdbname}_seg{seg_counter}'
        #    mfiles[pdbout] = '\n'.join(s.get_PDB())

        #    s.pop_last_filter()
        #    seg_counter += 1
    s.clear_filters()
    return s, residue_segments












#def split_segs(pdbdata, dssps=None, minimum=2, dssp_out=None, destination=''):
#    """
#    """
#    s = Structure(pdbdata)
#    s.build()
#
#    # this will ignore the iCode
#    residues = [int(i) for i in dict.fromkeys(s.filtered_atoms[:, col_resSeq])]
#
#    # returns slices
#    #segments = group_consecutive_ints(residues)
#    segments = group_runs(residues)
#
#    # removes ligands
#    above_2 = filter(
#        lambda x: len(x) > minimum,
#        segments,
#        )
#
#    seg_counter = 0
#    for seg in above_2:
#        s.add_filter(lambda x: int(x[col_resSeq]) in seg)
#
#        if set(s.filtered_atoms[:, col_record]) == {'HETATM'}:
#            s.pop_last_filter()
#            continue
#
#
#        # split regions with missing backbone
#        s.add_filter_backbone(minimal=True)
#
#        try:
#                #identify_backbone_gaps(s.filtered_atoms)
#            backbone_segs_in_resSeq_sets = \
#                list(filter(
#                    lambda x: len(x) > minimum,
#                    get_slice(s.filtered_atoms)
#                    ))
#        except Exception as err:
#            log.error(traceback.format_exc())
#            log.error(repr(err))
#            log.error(f'error in {pdbdata}')
#            seg_counter += 1
#            continue
#        finally:
#            s.pop_last_filter()
#
#
#        for resSeq_set in backbone_segs_in_resSeq_sets:
#
#            s.add_filter(lambda x: x[col_resSeq] in resSeq_set)
#
#            pdbout = f'{pdbdata.stem}_seg{seg_counter}'
#            fout_seg = Path(destination, pdbout).with_suffix('.pdb')
#
#            try:
#                s.write_PDB(fout_seg)
#            except EXCPTS.EmptyFilterError as err:
#                log.error(f'Empty filter for: {repr(err)}')
#            except EXCPTS.IDPConfGenException as err:
#                log.error(S('* Something went wrong with {}: {}', pdbdata, repr(err)))
#                log.debug(traceback.format_exc())
#
#            s.pop_last_filter()
#
#
#            # these are aligned
#            dssp_data = dssps[pdbdata.stem]
#            _fasta = dssp_data['fasta']
#            _dssp = dssp_data['dssp']
#            _res = dssp_data['resids'].split(',')
#
#
#            ffasta = []
#            ddssp = []
#            rres = []
#            for f, d, r in zip(_fasta, _dssp, _res):
#                if r in resSeq_set:
#                    ffasta.append(f)
#                    ddssp.append(d)
#                    rres.append(r)
#
#            dssp_out[pdbout] = {
#                'fasta': ''.join(ffasta),
#                'dssp': ''.join(ddssp),
#                'residues': ','.join(rres),
#                }
#
#            #
#            #dssp_slicing = slice(
#            #    residues.index(int(resSeq_set[0])),
#            #    residues.index(int(resSeq_set[-1])) + 1,
#            #    None,
#            #    )
#            #print(dssp_slicing)
#
#            #print(len(resSeq_set))
#            #print(len(dssps[pdbdata.stem][dssp_slicing]))
#            #print(dssps[pdbdata.stem][dssp_slicing])
#
#            #if not len(resSeq_set) == len(dssps[pdbdata.stem][dssp_slicing]):
#            #    log.error(traceback.format_exc())
#            #    log.error(f'error in {pdbdata}')
#            #    seg_counter += 1
#            #    continue
#            #dssp_out[pdbout] = dssps[pdbdata.stem][dssp_slicing]
#            seg_counter += 1
#
