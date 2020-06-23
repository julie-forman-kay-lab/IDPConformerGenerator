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

def download_pdbs_to_folder(destination, items, func=None, **kwargs):
    """
    Download PDBs to folder.

    Uses :func:`idpconfgen.libs.libmulticore.pool_function_in_chunks`
    """
    dest = make_destination_folder(destination)
    for chunk in pool_function_in_chunks(consume_iterable_in_list, items, func, **kwargs):
        for result in chunk:
            for fname, data in result:
                with open(Path(dest, fname), 'w') as fout:
                    fout.write(data)


def download_pdbs_to_tar(destination, items, func=None, **kwargs):
    """
    Download PDBs to tarfile.

    Uses :func:`idpconfgen.libs.libmulticore.pool_function_in_chunks`
    """
    # append 'a' here combos with the
    # read_PDBID_from_source(destination), in cli_pdbdownloader
    _exists = {True: 'a', False: 'w'}
    dests = str(destination)
    with tarfile.open(dests, mode=_exists[destination.exists()]) as tar:
        for chunk in pool_function_in_chunks(consume_iterable_in_list, items, func, **kwargs):
            for result in chunk:
                for fout, data in result:
                    save_file_to_tar(tar, fout, data)


def download_dispacher(func, destination, *args, **kwargs):
    """Dispaches the appropriate download env based on `destination`."""
    dispacher = download_dispacher
    suffix = destination.suffix
    return dispacher[suffix](destination, *args, func=func, **kwargs)



# NOT USED
def filter_structure(pdb_path, **kwargs):
    """
    Download a PDB/CIF structure chains.

    Parameters
    ----------
    pdbid : tuple of 2-elements
        0-indexed, the structure ID at RCSB.org;
        1-indexed, a list of the chains to download.

    **kwargs : as for :func:`save_structure_by_chains`.
    """
    pdbid = PDBIDFactory(pdb_path)
    pdbname = pdbid.name
    chains = pdbid.chain

    save_structure_by_chains(
        pdb_path,
        pdbname,
        chains=chains,
        **kwargs,
        )

download_dispacher = {
    '': download_pdbs_to_folder,
    '.tar': download_pdbs_to_tar,
    }


# NOT USED
def mkdssp_w_split_w_save(*args, destination=None, **kwargs):
     for fname, split_dssp, split_pdb_bytes in mkdssp_w_split(*args, **kwargs):
        save_pairs_to_disk(((f'{fname}.pdb', split_pdb_bytes),), destination=destination)
        yield fname, split_dssp

@contextmanager
def try_to_write(data, fout):
    """Context to download."""
    try:
        yield
    except Exception as err:
        log.debug(traceback.format_exc())
        log.error(S('error found for {}: {}', fout, repr(err)))
        erp = Path(fout.myparents(), 'errors')
        erp.mkdir(parents=True, exist_ok=True)
        p = Path(erp, fout.stem).with_suffix('.structure')
        if p.exists():
            return
        else:
            try:
                p.write_bytes(data)
            except TypeError:
                p.write_text(data)


def eval_chain_case(chain, chain_set):
    """
    Evalute chain case.

    Tests whether chain exists in chain_set in both UPPER and LOWER cases.

    Returns
    -------
    str
        The found chain.
    """
    if chain in chain_set:
        return chain
    else:
        for i in range(len(chain) + 1):
            cl = f'{chain[:i]}{chain[i:].lower()}'
            if cl in chain_set:
                return cl

    raise ValueError(f'Neither {chain}/{cl} are in set {chain_set}')


def identify_backbone_gaps(atoms):
    """Atoms is expected already only minimal backbone."""
    # this is a tricky implementation, PRs welcomed :-)
    if not set(atoms[:, col_name]).issubset({'N', 'CA', 'C'}):
        raise EXCPTS.PDBFormatError(errmsg='Back bone is not subset')

    resSeq, slices = zip(*group_by(atoms[:, col_resSeq]))
    print(slices)
    if not all(isinstance(i, slice) for i in slices):
        raise TypeError('Expected slices found something else')

    # calculates the length of each slice
    # each slice corresponds to a residue
    slc_len = np.array([slc.stop - slc.start for slc in slices])
    print(slc_len)

    # identifies indexes where len slice differ from three.
    # considering that atoms contains only minimal backbone atoms
    # (N, CA, C) this would mean that some atoms are missing
    aw = np.argwhere(slc_len != 3)

    # splits the slices according to the indexes where != 3
    s = np.split(slices, aw[:, 0])

    # np.split keeps the index where it was split, contrarily to str.split
    # here I filter out those unwanted indexes
    # parsed_split should contain a list of list with only the slices
    # for residues with 3 backbone atoms
    # the first index is a special case because of np.split
    parsed_split = \
        [list(i) for i in s[0:1] if i[1:].size > 0] \
        + [list(i[1:]) for i in s[1:] if i[1:].size > 0]

    # concatenate residue slices in segments
    segment_slices = [slice(i[0].start, i[-1].stop) for i in parsed_split]

    # returns a list of array views
    # return [atoms[seg, :] for seg in segment_slices]
    return [
        list(dict.fromkeys(atoms[seg, col_resSeq]))
        for seg in segment_slices
        ]

# DEPRECATED
def group_consecutive_ints(sequence):
    """Groupd consecutive ints."""
    prev = sequence[0]
    start = 0
    slices = []
    for i, integer in enumerate(sequence[1:], start=1):
        if integer - prev not in (0, 1):
            slices.append(slice(start, i))
            start = i
        prev = integer
    else:
        slices.append(slice(start, i + 1))
    return slices


def list_index_to_array(
        list_,
        start=None,
        stop=None,
        step=None,
        sObj=None,
        index=None,
        ):
    """
    Extract slices of strings in lists to an array.

    At least one the named parameters must be provided.

    In case more than one parameter set is provided:
        the `start`, `stop` or `step` trio have the highest priority,
        followed by `sObj` and last `index`.

    Raises
    ------
    ValueError
        If no named parameter is provided.

    ValueError
        If sObj is provided and is not a ``slice`` object.
    """
    if any((i is not None for i in (start, stop, step))):
        index = slice(start, stop, step)

    elif sObj is not None:
        if isinstance(sObj, slice):
            index = sObj
        else:
            raise ValueError(
                'sObj must be slice type: {} given.'.format(type(sObj))
                )

    elif index is not None:
        index = index

    else:
        raise ValueError('At least one of the index args most be provided.')

    len_of_slice = len(list_[0][index])

    array = np.empty(
        len(list_),
        dtype='<U{}'.format(len_of_slice),
        )

    for i, line in enumerate(list_):
        array[i] = line[index]

    return array


# DEPRECATED
def _concatenate_ss_from_dsspparsers(dsspparsers):
    output = []
    for dsspparser in dsspparsers:
        output.append(
            '{}|{}'.format(
                dsspparser.pdbid,
                ''.join(dsspparser.ss),
                )
            )

    output.sort(key=lambda x: x.split('|')[0])

    return output


# DEPECRATED
def read_pipe_file(text):
    """
    Read IDPConfGen pipe file.

    Pipe files have the following structure:

    CODE|DATA
    CODE|DATA
    (...)

    Parameters
    ----------
    text : str
        The crude text to read.

    Returns
    -------
    dict {str: str}
        Of codes and data.
    """
    lines = text.split('\n')
    elines = filter(bool, lines)
    return {c: s for line in elines for c, s in (line.split('|'),)}


# DEPRECATED
@libcheck.kwargstype((str, Path, type(None)))
def export_ss_from_DSSP(*dssp_task_results, output='dssp.database'):
    """
    Export secondary structure information from :class:`DSSPParser` results.

    Parameters
    ----------
    dssp_task_results : multiple :class:`DSSPParser`
    """
    output_data = _concatenate_ss_from_dsspparsers(dssp_task_results)

    output_data = '\n'.join(output_data) + '\n'

    if output:
        opath = Path(output)
        opath.myparents().mkdir(parents=True, exist_ok=True)
        opath.write_text(output_data)
    else:
        sys.stdout.write(output_data)

