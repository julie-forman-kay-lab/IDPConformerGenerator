"""Implementations that have been droped or placed apart for now."""

def split_residue_segments(atom_lines):
    """
    Splits a sequence of atoms into consecutive residue segments.

    Parameters
    ----------
    atom_lines : np.ndarray of shape (N, col*)
        A :class:`Structure.data_array` like.
    """
    segments = [[]]
    prev = int(atom_lines[0, col_resSeq])
    for line in atom_lines:
        curr = int(line[col_resSeq])
        if 0 <= curr - prev < 2:
            segments[-1].append(line)
        else:
            segments.append([])
            segments[-1].append(line)
        prev = curr
    return segments


def eval_bb_gaps(structure, ignore_hetatms_segments=True, minimum=2):
    """
    Split PDBs in continuous backbone chunkstructure.

    Parameters
    ----------
    pdbdata : str
        Text of the PDB file.

    minimum : int
        Minimum backbone to accept.
    """
    hetatm_set = {'HETATM'}

    # removes segments shorter than the minimum allowed
    # creates a list because the next loop will operate on the filters
    # and .consecutive_residues does uses filters
    above_minimum = list(filter(
        lambda segment: len(segment) > minimum,
        structure.consecutive_residues,
        ))

    #seg_counter = 0
    residue_segments = []
    for segment in above_minimum:
        structure.clear_filters()
        structure.add_filter(lambda x: int(x[col_resSeq]) in segment)

        # discard segments that might be composed only of HETATM
        if ignore_hetatms_segments:
            tmp_filtered = structure.filtered_atoms
            if set(tmp_filtered[:, col_record]) == hetatm_set:
                continue

        # identify backbone gaps
        structure.add_filter_backbone(minimal=True)
        residue_segments.extend(
            list(filter(
                lambda segment: len(segment) > minimum,
                get_segments_based_on_backbone_continuity(structure.filtered_atoms),
                ))
            )
    structure.clear_filters()
    return residue_segments


def get_segments_based_on_backbone_continuity(atoms):
    """
    Split backbones in chunks of continuity.

    Consideres only minimal backbone

    Have to explain big thing here.
    """
    ref = _minimal_bb_atoms
    start = 0
    idx = 0
    slices = []
    max_size = atoms.shape[0]

    bb = atoms[:, col_name]
    resis = np.core.defchararray.add(atoms[:, col_resSeq], atoms[:, col_iCode])

    while idx < max_size:
        a = list(bb[idx: idx + 3])
        bb_continuity_lost = a != ref or len(set(resis[idx: idx + 3])) > 1
        if bb_continuity_lost:
            slices.append(slice(start, idx, None))
            idx += 1
            start = idx
            while idx <= max_size:
                a = list(bb[idx: idx + 3])
                bb_continuity_restored = \
                    a == ref and len(set(resis[idx: idx + 3])) == 1
                if bb_continuity_restored:
                    start = idx
                    idx += 3
                    break  # back to the main while
                idx += 1
        else:
            idx += 3
    else:
        slices.append(slice(start, idx + 1, None))

    return [list(dict.fromkeys(atoms[seg, col_resSeq])) for seg in slices]

def get_PDB_from_residues(structure, residue_segments):
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
        yield list(structure.get_PDB())


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
    To do.
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

