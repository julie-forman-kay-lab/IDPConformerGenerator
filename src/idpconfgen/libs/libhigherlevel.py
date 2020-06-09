"""
Higher level functions.

Function which operate with several libraries
and are defined here to avoid circular imports.
"""
from collections import defaultdict

from idpconfgen import Path
from idpconfgen.libs.libstructure import Structure, eval_bb_gaps, get_PDB_from_residues, col_resSeq
from idpconfgen.libs.libparse import group_runs
from idpconfgen.libs.libtimer import record_time
from idpconfgen.libs.libpdb import atom_resSeq


def split_backbone_segments(
        pdbid,
        mfiles=None,
        mdata=None,
        minimum=2,
        sscalc_data=None,
        ):
    """
    """
    pdbname = Path(pdbid[0]).stem
    pdbdata = pdbid[1]
    pdb_dssp_data = sscalc_data[pdbname]
    residues = pdb_dssp_data['resids'].split(',')
    dssp_residues = [int(i) for i in residues]
    residues_groups = [[str(i) for i in seg] for seg in group_runs(dssp_residues)]
    structural_data = [k for k in pdb_dssp_data.keys() if k != 'resids']
    #splitdata = defaultdict(dict)
    for i, segment in enumerate(residues_groups):
        key = f'{pdbname}_seg{i}'
        for datatype in structural_data:
            mdata[key][datatype] = \
                ''.join(
                    c
                    for res, c in zip(residues, pdb_dssp_data[datatype])
                    if res in segment
                    )
        mdata[key]['resids'] = ','.join(segment)

        #structure = Structure(pdbdata)
        #structure.build()
        #structure.add_filter(lambda x: x[col_resSeq] in segment)
        #_mfiles[f'{key}.pdb'] = '\n'.join(structure.get_PDB())
        segset = set(segment)
        mfiles[f'{key}.pdb'] = b'\n'.join(line for line in pdbdata.split(b'\n') if line[atom_resSeq].strip() in segset)


def _split_backbone_segments(
        pdbid,
        mfiles,
        sscalc_data=None,
        mdata=None,
        minimum=2,
        ):
    """
    Parameters
    ----------
    pdbid : tuple (str, str)

    mfiles : `multiprocessing.Manager.dict`

    mdata : `multiprocessing.Manager.dict`, optional

    minimum : int, optional
        The minimum size of each segment.
        Defaults to 2.

    sscalc_data : dict, optional
    """
    pdbname = Path(pdbid[0]).stem
    pdbdata = pdbid[1]
    structure = Structure(pdbdata)
    structure.build()

    # this might be slightly slower but it definitively more modular
    # and testable
    # `structure_segments` are in residue number (str)
    structure_segments = eval_bb_gaps(structure, minimum=2)
    assert isinstance(structure_segments, list)

    # here I have to put them in the mfiles
    if sscalc_data:
        dssp_segments, structure_segments = split_sscalc_data(
            pdbname,
            sscalc_data,
            structure_segments,
            )
        assert len(dssp_segments) == len(structure_segments), pdbname
        mdata.update(dssp_segments)

    splitted_segments = \
        list(get_PDB_from_residues(structure, structure_segments))

    assert len(structure_segments) == len(splitted_segments), pdbname
    for i, txt in enumerate(splitted_segments):
        mfiles[f'{pdbname}_seg{i}'] = txt


def split_sscalc_data(pdbname, data, segments):
    """
    """
    pdb_dssp_data = data[pdbname]
    structural_data = [k for k in pdb_dssp_data.keys() if k != 'resids']

    # residues as identified by DSSP
    residues = pdb_dssp_data['resids'].split(',')
    assert all(s.isdigit() for s in residues), pdbname

    splitdata = defaultdict(dict)
    segs_in_dssp = []

    for i, segment in enumerate(segments):

        assert all(s.isdigit() for s in segment), pdbname
        # this is done with list instead of sets because I want to keep order
        # two reduction step

        # reduces residues in DSSP to the residues in the segment
        reduction = [res for res in residues if res in segment]

        # some residues can be in the DSSP and not on the segment
        # residues in the segment must be reduced accordingly
        segs_in_dssp.append([res for res in segment if res in reduction])
        assert len(reduction) == len(segs_in_dssp[-1])

        key = f'{pdbname}_seg{i}'
        for datatype in structural_data:
            splitdata[key][datatype] = \
                ''.join(
                    c
                    for res, c in zip(residues, pdb_dssp_data[datatype])
                    if res in reduction
                    )

        assert all(
            len(value) == len(reduction)
            for value in splitdata[key].values()
            ), pdbname

        splitdata[key]['resids'] = ','.join(reduction)

    return splitdata, segs_in_dssp

