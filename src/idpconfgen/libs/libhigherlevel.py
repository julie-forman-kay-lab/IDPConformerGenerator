"""
Higher level functions.

Function which operate with several libraries
and are defined here to avoid circular imports.
"""
import re
from contextlib import suppress
from functools import partial, reduce

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core.definitions import (
    aa3to1,
    bgeo_CaC,
    bgeo_CaCNp1,
    bgeo_CaCO,
    bgeo_Cm1NCa,
    bgeo_CNp1,
    bgeo_CO,
    bgeo_NCa,
    bgeo_NCaC,
    blocked_ids,
    )
from idpconfgen.core.exceptions import IDPConfGenException, PDBFormatError
from idpconfgen.libs.libcalc import (
    calc_angle_njit,
    calc_torsion_angles,
    rrd10_njit,
    )
from idpconfgen.libs.libio import (
    concatenate_entries,
    read_PDBID_from_source,
    save_pairs_to_disk,
    )
from idpconfgen.libs.libmulticore import (
    consume_iterable_in_list,
    flat_results_from_chunk,
    pool_function_in_chunks,
    )
from idpconfgen.libs.libparse import group_by
from idpconfgen.libs.libpdb import PDBList, atom_name, atom_resSeq
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resName,
    cols_coords,
    )
from idpconfgen.logger import S, T, init_files, report_on_crash


# USED OKAY!
def download_pipeline(func, logfilename='.download'):
    """
    Context pipeline to download PDB/mmCIF files.

    Exists because fetching and download filtered PDBs shared the same
    operational contextm, only differing on the function which
    orchestrates the download process.

    Parameters
    ----------
    func : function
        The actual function that orchestrates the download operation.

    logfilename : str, optional
        The common stem name of the log files.
    """
    LOGFILESNAME = logfilename

    def main(
            pdbids,
            chunks=5_000,
            destination=None,
            ncores=1,
            update=False,
            **funckwargs,
            ):
        """Run main script logic."""
        init_files(log, LOGFILESNAME)

        #
        log.info(T('reading input PDB list'))

        pdblist = PDBList(concatenate_entries(pdbids))

        log.info(
            f"{S(str(pdblist))}\n"
            f"{S('done')}\n"
            )

        #
        log.info(T('Filtering input'))
        destination = destination or Path.cwd()
        log.info(
            f"{S(f'from destination: {destination}')}\n"
            f"{S('and other sources...')}"
            )

        # comparison block
        def diff(first, other):
            return first.difference(other)

        remove_from_input = [
            read_PDBID_from_source(destination),
            PDBList(blocked_ids),
            ]

        # yes, there are just two items in remove_from_input, why use reduce?
        # what if more are added in the future? :-P
        # the engine is already created
        pdblist_comparison = reduce(diff, remove_from_input, pdblist)
        log.info(S(f'Found {str(pdblist_comparison)} to download'))
        #

        something_to_download = len(pdblist_comparison) > 0
        if something_to_download and update:

            # the function to be used in multiprocessing
            consume_func = partial(consume_iterable_in_list, func, **funckwargs)

            # context to perform a dedicated report in case function fails
            # partial is mandatory because decorators won't pickle
            # in this way, crashes are reported for files, crashed files
            # ignored, and the process continues
            execute = partial(
                report_on_crash,
                consume_func,
                ROC_exception=Exception,
                ROC_prefix='download_pipeline',
                )

            # convinient partial
            execute_pool = partial(
                pool_function_in_chunks,
                execute,
                list(pdblist_comparison.name_chains_dict.items()),
                ncores=ncores,
                chunks=chunks,
                )

            flat_results_from_chunk(
                execute_pool,
                save_pairs_to_disk,
                destination=destination,
                )

            log.info(T('Reading UPDATED destination'))
            pdblist_updated = read_PDBID_from_source(destination)
            pdblist_up_comparison = pdblist.difference(pdblist_updated)
            log.info(S(f'{str(pdblist_up_comparison)}'))
            if len(pdblist_up_comparison) > 0:
                log.info(S(
                    'There are PDBIDs not downloaded\n.'
                    'Those IDs have been registered in the '
                    f'{LOGFILESNAME}.debug file.'
                    ))
                log.debug('\n'.join(str(_id) for _id in pdblist_up_comparison))

        elif not something_to_download and update:
            log.info(S('There is nothing to download.'))
            log.info(S(
                'All requested IDs are already at '
                'the destination folder.'
                ))

        log.info(T('PDB Downloader finished'))
        return

    return main


# USED OKAY
def extract_secondary_structure(
        pdbid,
        ssdata,
        atoms='all',
        minimum=0,
        structure='all',
        ):
    """
    Extract secondary structure elements from PDB data.

    Parameters
    ----------
    pdbid : tuple
        Where index 0 is the PDB id code, for example `12AS.pdb`, or
        `12AS_A`, or `12AS_A_seg1.pdb`.
        And, index 1 is the PDB data itself in bytes.

    ssdata : dict
        Dictionary containing the DSSP information. Must contain a key
        equal to `Path(pdbid).stem, where `dssp` key contains the DSSP
        information for that PDB.

    atoms : str or list of str or bytes, optional
        The atom names to keep.
        Defaults to `all`.

    minimum : int, optional
        The minimum size a segment must have in order to be considered.
        Defaults to 0.

    structure : str or list of chars
        The secondary structure character to separate.
        Multiple can be given in the form of a list.
    """
    # caused problems in the past
    assert atoms != ['all']
    assert structure != ['all']

    pdbname = Path(pdbid[0]).stem
    pdbdata = pdbid[1].split(b'\n')

    # gets PDB computed data from dictionary
    try:
        pdbdd = ssdata[pdbname]
    except KeyError:
        pdbdd = ssdata[f'{pdbname}.pdb']

    if structure == 'all':
        ss_to_isolate = set(pdbdd['dssp'])
    else:
        ss_to_isolate = set(structure)

    # general lines filters
    line_filters = []
    LF_append = line_filters.append
    LF_pop = line_filters.pop

    # in atoms filter
    if atoms != 'all':
        with suppress(AttributeError):  # it is more common to receive str
            atoms = [c.encode() for c in atoms]
        line_filters.append(lambda x: x[atom_name].strip() in atoms)

    dssp_slices = group_by(pdbdd['dssp'])
    # DR stands for dssp residues
    DR = [c for c in pdbdd['resids'].encode().split(b',')]

    for ss in ss_to_isolate:

        ssfilter = (slice_ for char, slice_ in dssp_slices if char == ss)
        minimum_size = (s for s in ssfilter if s.stop - s.start >= minimum)

        for counter, seg_slice in enumerate(minimum_size):

            LF_append(lambda x: x[atom_resSeq].strip() in DR[seg_slice])  # noqa: E501, B023

            pdb = b'\n'.join(
                line for line in pdbdata if all(f(line) for f in line_filters)
                )
            LF_pop()

            yield f'{pdbname}_{ss}_{counter}.pdb', pdb

            counter += 1


def get_torsionsJ(
        fdata,
        decimals=5,
        degrees=False,
        hn_terminal=True,
        hn_labels=('H', 'H1'),
        proline_value=np.nan,
        ):
    """
    Calculate HN-CaHA torsion angles from a PDB/mmCIF file path.

    Needs atom labels: H or H1, N, CA, HA or HA2 (Glycine).

    Parameters
    ----------
    decimals : int
        The decimal number to round the result.

    degrees : bool
        Whether or not to return values as degrees. If `False` returns
        radians.

    hn_terminal : bool
        If the N-terminal has no hydrogens, flag `hn_terminal` should be
    provided as `False`, and the first residue will be discarded.
    If `True` expects N-terminal to have `H` or `H1`.

    Returns
    -------
    np.ndarray
        The NH-CaHA torsion angles for the whole protein.
        Array has the same length of the protein if N-terminal has H,
        otherwise has length of protein minus 1.

    Notes
    -----
    Not optimized for speed. Not slow either.
    """
    # reads the structure file
    structure = Structure(fdata)
    structure.build()
    data = structure.data_array

    # to adjust data to calc_torsion_angles(), we consider the CD of Prolines
    # later we will DELETE those entries
    protons_and_proline = np.logical_or(
        np.isin(data[:, col_name], hn_labels),
        np.logical_and(data[:, col_resName] == 'PRO', data[:, col_name] == 'CD')
        )

    hn_idx = 0 if hn_terminal else 1

    # some PDBs may not be sorted, this part sorts atoms properly before
    # performing calculation
    hs = data[protons_and_proline, :]

    n = data[data[:, col_name] == 'N', :][hn_idx:, :]
    ca = data[data[:, col_name] == 'CA', :][hn_idx:, :]
    ha = data[np.isin(data[:, col_name], ('HA', 'HA2')), :][hn_idx:, :]

    # expects N-terminal to have `H` or `H1`
    assert hs.shape == n.shape == ca.shape == ha.shape, (
        'We expected shapes to be equal. '
        'A possible reason is that the presence/absence of protons in '
        'the N-terminal does not match the flag `hn_terminal`. '
        f'Shapes found are as follow: {hs.shape, n.shape, ca.shape, ha.shape}'
        )

    n_data = np.hstack([hs, n, ca, ha]).reshape(hs.shape[0] * 4, hs.shape[1])

    coords = n_data[:, cols_coords].astype(np.float64)

    # notice that calc_torsion_angles() is designed to accepted sequential
    # atoms for which torsions can be calculated. In this particular case
    # because of the nature of `n_data`, the only torsion angles that will have
    # physical meaning are the ones referrent to HN-CaHA, which are at indexes
    # 0::4
    torsions = calc_torsion_angles(coords)[0::4]

    # not the fastest approach
    # increase performance when needed
    tfasta = structure.fasta

    # assumes there is a single chain
    fasta = list(tfasta.values())[0][hn_idx:]
    pro_idx = [m.start() for m in re.finditer('P', fasta)]

    # assigns nan to proline positions
    torsions[pro_idx] = proline_value

    if degrees:
        torsions = np.degrees(torsions)

    return np.round(torsions, decimals)


def get_torsions(fdata, degrees=False, decimals=3, validate=True):
    """
    Calculate torsion angles from structure.

    Corrects for labels not sorted according to (N, CA, C). But does not
    correct for non-sorted residues (do these exist?).

    Calculates torsion angles with `:func:libcalc.calc_torsion_angles`.

    Parameters
    ----------
    fdata : str, bytes or Path
        A path to the structure file, or the string representing
        the file.
        In fact, accepts any type `:class:libstructure.Structure` would
        accept.

    degrees : bool, optional
        Whether to return torsion angles in degrees or radians.

    decimals : int, optional
        The number of decimals to return.
        Defaults to 3.

    validate : bool
        Validates coordinates and labels according to general
        expectations. Expectations are as described by functions:
        - :func:libhigherlevel.validate_backbone_labels_for_torsion
        - :func:libhigherlevel.validate_coords_for_backbone_torsions

        If `False`, does not perform validation. Be sure to provide PDBs
        what will output meaningful results.

    Return
    ------
    np.ndarray
        As described in `:func:libcalc.calc_torsion_angles()` but with
        `decimals` and `degrees` applied.

    See Also
    --------
    libhigherlevel.validate_backbone_labels_for_torsion
    libhigherlevel.validate_coords_for_backbone_torsions
    libcalc.calc_torsion_angles
    """
    structure = Structure(fdata)
    structure.build()
    structure.add_filter_backbone(minimal=True)

    data = structure.filtered_atoms
    names = data[:, col_name]
    coords_raw = structure.coords

    n_mask = names == 'N'
    ca_mask = names == 'CA'
    c_mask = names == 'C'

    n = coords_raw[n_mask, :]
    ca = coords_raw[ca_mask, :]
    c = coords_raw[c_mask, :]

    n_names = names[n_mask]
    ca_names = names[ca_mask]
    c_names = names[c_mask]

    try:
        # https://stackoverflow.com/questions/5347065
        labels = np.empty((n_names.size * 3), dtype=ca_names.dtype)
        labels[0::3] = n_names
        labels[1::3] = ca_names
        labels[2::3] = c_names
    except ValueError as err:
        errmsg = (
            'Labels do not match expectation. '
            'Some labels possibly missing.'
            )
        raise IDPConfGenException(errmsg) from err

    try:
        coords = np.empty((n.shape[0] * 3, 3), dtype=np.float64)
        coords[0::3, :] = n
        coords[1::3, :] = ca
        coords[2::3, :] = c
    except ValueError as err:
        errmsg = (
            'Coordinates do not match expectation. '
            'Some possibly missing.'
            )
        raise IDPConfGenException(errmsg) from err

    coords_distances = np.linalg.norm(coords[:-1, :] - coords[1:, :], axis=1)
    assert coords_distances.size == coords.shape[0] - 1
    if np.any(coords_distances > 2.1):
        errmsg = (
            'Chain is broken. Found distance above 2.1A for '
            'consecutive atoms.'
            )
        raise IDPConfGenException(errmsg)

    # DEPRECATE / REMOVE
    # if validate:
    #    # validates structure data
    #    # rare are the PDBs that produce errors, still they exist.
    #    # errors can be from a panoply of sources, that is why I decided
    #    # not to attempt correcting them further and instead ignore and report.
    #    validation_labels = validate_backbone_labels_for_torsion(labels)
    #    if not validation_labels:
    #        errmsg = (
    #            'Found errors on backbone label consistency: '
    #            f'{validation_error}\n'
    #            )
    #        err = IDPConfGenException(errmsg)
    #        # if running through cli_torsions, `err` will be catched and
    #        # reported by logger.report_on_crash
    #        raise err

    #    validation_coords = validate_coords_for_backbone_torsions(coords)
    #    if not validation_coords:
    #        errmsg = (
    #            'Found errors on coords consistency: '
    #            f'{validation_error}\n'
    #            )
    #        err = IDPConfGenException(errmsg)
    #        raise err

    torsions = calc_torsion_angles(coords)

    if degrees:
        torsions = np.degrees(torsions)

    return np.round(torsions, decimals)


def cli_helper_calc_torsions(fname, fdata, **kwargs):
    """
    Help `cli_torsion` to operate.

    Returns
    -------
    dict
        key: `fname`
        value: -> dict, `phi`, `phi`, `omega` -> list of floats
    """
    torsions = get_torsions(fdata, **kwargs)
    CA_C, C_N, N_CA = get_separate_torsions(torsions)
    return fname, {'phi': N_CA, 'psi': CA_C, 'omega': C_N}


def cli_helper_calc_torsionsJ(fdata_tuple, **kwargs):
    """Help cli_torsionsJ.py."""
    return fdata_tuple[0], get_torsionsJ(fdata_tuple[1], **kwargs)


def get_separate_torsions(torsions_array):
    """
    Separate torsion angles according to the protein backbone concept.

    Considers torsion angles for bonds in between atom pairs:
        - CA - C
        - C - N
        - N - CA

    Backbone obeys the order: N-CA-C-N-CA-C(...)

    And the first value corresponds to a CA-C pair, because the
    first N-CA pair of the protein backbone has no torsion angle.
    """
    assert torsions_array.ndim == 1
    assert torsions_array.size % 3 == 0

    CA_C = torsions_array[::3].tolist()
    C_N = torsions_array[1::3].tolist()
    N_CA = torsions_array[2::3].tolist()

    assert len(CA_C) == len(C_N) == len(N_CA)
    return CA_C, C_N, N_CA


# DEPRECATE / REMOVE
def validate_coords_for_backbone_torsions(coords, minimum=2):
    """
    Validate coords for torsions.

    Does NOT validate with values have physical sense. Validates only
    valid input for `:func:libcalc.calc_torsion_angles`.

    *Validations performed:*

        * coords are two-dimensional arrays
        * coords have 3 values in the first dimension (XYZ) (shape[1])
        * number of coords is multiple of 3.

    Returns
    -------
    str
        A string explaining the error if an error is found.
        An empty string if no error is found.
    """
    if len(coords.shape) != 2:
        return 'Expected two dimensions {len(coords.shape)} found.'

    if coords.shape[1] != 3:
        return (
            'Expected 3 values for XYZ coordinates. '
            f'Found {coords.shape[1]}'
            )

    if coords.shape[0] % 3:
        return (
            'Expected multiple of 3 (N, CA, C), some coordinates are '
            'missing.'
            )

    return ''


# DEPRECATE / REMOVE
def validate_backbone_labels_for_torsion(labels, minimum=2):
    """
    Validate labels for torsion angle calculation.

    Assumes labels are aligned with their corresponding coordinates.
    Yet, coordinates have no scope in this function.

    Excepts only the mininal backbone labels, these are: N, CA, and C.

    Parameters
    ----------
    labels : np.array of shape (N,) or alike
        Where N % 3 equals 0.

    minimum : int
        The minimum number of residues to consider valid.
    """
    if len(labels) / 3 < minimum:
        return 'Segment too short.'

    if labels[0] != 'N':
        return 'The first atom is not N, it should be!'

    if len(labels) % 3:
        return 'Number of backbone atoms is not module of 3.'

    if set(labels) != {'N', 'C', 'CA'}:
        return 'There are atoms other than N, C and CA.'

    return ''


def get_bond_geos(fdata):
    """
    Calculate bond angles from structure.

    Parameters
    ----------
    fdata : data to :pyclass:`idpconfgen.libstructure.Structure`

    degrees : bool, optional
        Defaults to False.

    decimals : int, optional
        Defaults to 3.
    """
    ALL = np.all
    CEQ = np.char.equal
    CO_LABELS = np.array(['CA', 'C', 'O', 'CA', 'C', 'O', 'CA'])
    NORM = np.linalg.norm

    s = Structure(fdata)
    s.build()
    s.add_filter_backbone(minimal=True)

    if s.data_array[0, col_name] != 'N':
        raise PDBFormatError(
            'PDB does not start with N. '
            f'{s.data_array[0, col_name]} instead.'
            )

    N_CA_C_coords = s.coords
    s.clear_filters()

    s.add_filter(lambda x: x[col_name] in ('CA', 'C', 'O'))
    CA_C_O_coords = s.coords
    co_minimal_names = s.filtered_atoms[:, col_name]

    bgeo_results = {
        bgeo_Cm1NCa: [],
        bgeo_NCaC: [],
        bgeo_CaCNp1: [],
        bgeo_CaCO: [],
        bgeo_NCa: [],
        bgeo_CaC: [],
        bgeo_CNp1: [],
        bgeo_CO: [],
        }

    for i in range(1, len(N_CA_C_coords) - 7, 3):
        idx = list(range(i, i + 7))

        # calc bend angles
        c = N_CA_C_coords[idx]
        Cm1_N = c[1] - c[2]
        Ca_N = c[3] - c[2]
        N_Ca = c[2] - c[3]
        C_Ca = c[4] - c[3]
        Ca_C = c[3] - c[4]
        Np1_C = c[5] - c[4]
        assert Cm1_N.shape == (3,)

        # calc bond lengths
        # need float for json.dump else float32
        NCa = float(NORM(N_Ca))
        CaC = float(NORM(Ca_C))
        CNp1 = float(NORM(Np1_C))

        # the angles here are already corrected to the format needed by the
        # builder, which is (pi - a) / 2
        Cm1_N_Ca = (np.pi - calc_angle_njit(Cm1_N, Ca_N)) / 2
        N_Ca_C = (np.pi - calc_angle_njit(N_Ca, C_Ca)) / 2
        Ca_C_Np1 = (np.pi - calc_angle_njit(Ca_C, Np1_C)) / 2

        co_idx = np.array(idx) - 1
        c = CA_C_O_coords[co_idx]
        Ca_C = c[3] - c[4]
        O_C = c[5] - c[4]

        CO = float(NORM(O_C))

        if not ALL(CEQ(co_minimal_names[co_idx], CO_LABELS)):
            log.info(S(
                'Found not matching labels '
                f'{",".join(co_minimal_names[co_idx])}'
                ))
            continue
        Ca_C_O = calc_angle_njit(Ca_C, O_C) / 2

        bgeo_results[bgeo_Cm1NCa].append(Cm1_N_Ca)
        bgeo_results[bgeo_NCaC].append(N_Ca_C)
        bgeo_results[bgeo_CaCNp1].append(Ca_C_Np1)
        bgeo_results[bgeo_CaCO].append(Ca_C_O)

        bgeo_results[bgeo_NCa].append(NCa)
        bgeo_results[bgeo_CaC].append(CaC)
        bgeo_results[bgeo_CNp1].append(CNp1)
        bgeo_results[bgeo_CO].append(CO)

    if len(set(map(len, bgeo_results.values()))) != 1:
        raise AssertionError('something wrong here, this is a bug')

    return bgeo_results


def cli_helper_calc_bgeo(fname, fdata, **kwargs):
    """
    Help `cli_bgeodb` to operate.

    Returns
    -------
    dict
        key: `fname`
        value -> dict, `Ca_C_Np1`, `Ca_C_O`, `Cm1_N_Ca`, `N_Ca_C`
    """
    bond_geometries = get_bond_geos(fdata, **kwargs)
    return fname, bond_geometries


def read_trimer_torsion_planar_angles(pdb, bond_geometry):
    """
    Create a trimer/torsion library of bend/planar angles.

    Given a PDB file:

    1) reads each of its trimers, and for the middle residue:
    2) Calculates phi/psi and rounds them to the closest 10 degree bin
    3) assign planar angles found for that residue to the trimer/torsion key.
    4) the planar angles are converted to the format needed by cli_build,
       which is that of (pi - angle) / 2.
    5) updates that information in `bond_gemetry`.

    Created key:values have the following form in `bond_geometry` dict::

        {
            'AAA:10,-30': {
                'Cm1_N_Ca': [],
                'N_Ca_C': [],
                'Ca_C_Np1': [],
                'Ca_C_O': [],
                }
            }


    Parameters
    ----------
    pdb : any input of `libstructure.Structure`
        The PDB/mmCIF file data.

    bond_geometry : dict
        The library dictionary to update.

    Returns
    -------
    None
    """
    ALL = np.all
    CEQ = np.char.equal
    TORSION_LABELS = np.array(['CA', 'C', 'N', 'CA', 'C', 'N', 'CA'])
    CO_LABELS = np.array(['CA', 'C', 'O', 'CA', 'C', 'O', 'CA'])
    aa3to1['MSE'] = 'M'  # seleno methionine

    s = Structure(pdb)
    s.build()
    s.add_filter_backbone(minimal=True)

    if s.data_array[0, col_name] != 'N':
        raise PDBFormatError(
            'PDB does not start with N. '
            f'{s.data_array[0, col_name]} instead.'
            )

    bb_minimal_names = s.filtered_atoms[:, col_name]
    bb_residue_names = s.filtered_atoms[:, col_resName]

    N_CA_C_coords = s.coords
    s.clear_filters()

    s.add_filter(lambda x: x[col_name] in ('CA', 'C', 'O'))
    CA_C_O_coords = s.coords
    co_minimal_names = s.filtered_atoms[:, col_name]

    # calc torsion angles
    for i in range(1, len(N_CA_C_coords) - 7, 3):

        idx = list(range(i, i + 7))

        _trimer = bb_residue_names[idx][0::3]

        try:
            trimer = ''.join(aa3to1[_t] for _t in _trimer)
        except KeyError:
            log.info(S(
                'trimer '
                f"{','.join(_trimer)}"
                ' not found. Skipping...'
                ))
            continue

        assert len(trimer) == 3
        del _trimer

        if not ALL(CEQ(bb_minimal_names[idx], TORSION_LABELS)):
            log.info(S(
                'Found non-matching labels: '
                f'{",".join(bb_minimal_names[idx])}'
                ))
            continue

        # selects omega, phi, and psi for the central residue
        rad_tor = np.round(calc_torsion_angles(N_CA_C_coords[idx])[1:3], 10)
        ptorsions = [rrd10_njit(_) for _ in rad_tor]

        assert len(ptorsions) == 2
        for angle in ptorsions:
            assert -180 <= angle <= 180, 'Bin angle out of expected range.'

        # TODO: better key
        tuple_key = trimer + ':' + ','.join(str(_) for _ in ptorsions)

        # calc bend angles
        c = N_CA_C_coords[idx]
        Cm1_N = c[1] - c[2]
        Ca_N = c[3] - c[2]
        N_Ca = c[2] - c[3]
        C_Ca = c[4] - c[3]
        Ca_C = c[3] - c[4]
        Np1_C = c[5] - c[4]
        assert Cm1_N.shape == (3,)

        # the angles here are already corrected to the format needed by the
        # builder, which is (pi - a) / 2
        Cm1_N_Ca = (np.pi - calc_angle_njit(Cm1_N, Ca_N)) / 2
        N_Ca_C = (np.pi - calc_angle_njit(N_Ca, C_Ca)) / 2
        Ca_C_Np1 = (np.pi - calc_angle_njit(Ca_C, Np1_C)) / 2

        _ = bond_geometry[tuple_key].setdefault(bgeo_Cm1NCa, [])
        _.append(Cm1_N_Ca)

        _ = bond_geometry[tuple_key].setdefault(bgeo_NCaC, [])
        _.append(N_Ca_C)

        _ = bond_geometry[tuple_key].setdefault(bgeo_CaCNp1, [])
        _.append(Ca_C_Np1)

        co_idx = np.array(idx) - 1

        if not ALL(CEQ(co_minimal_names[co_idx], CO_LABELS)):
            log.info(S(
                'Found not matching labels '
                f'{",".join(co_minimal_names[co_idx])}'
                ))
            continue

        c = CA_C_O_coords[co_idx]
        Ca_C = c[3] - c[4]
        O_C = c[5] - c[4]

        Ca_C_O = calc_angle_njit(Ca_C, O_C)
        _ = bond_geometry[tuple_key].setdefault('Ca_C_O', [])
        _.append(Ca_C_O / 2)

    return


def convert_bond_geo_lib(bond_geo_db):
    """
    Convert bond geometry library to bond type first hierarchy.

    Restructure the output of `read_trimer_torsion_planar_angles` such
    that the main keys are the bond types, followed by the main residue,
    the pairs in the trimer, and, finally, the torsion angles.

    Parameters
    ----------
    bond_geo_db : dict
        The output of `read_PDBID_from_source`.

    Returns
    -------
    dict
    """
    converted = {}
    for key in bond_geo_db:
        trimer, torsion = key.split(':')
        letters = list(trimer)

        for at in bond_geo_db[key].keys():

            atkey = converted.setdefault(at, {})

            main = atkey.setdefault(letters[1], {})
            pairs = main.setdefault(''.join(letters[::2]), {})
            tor = pairs.setdefault(torsion, [])
            tor.extend(bond_geo_db[key][at])

    return converted


def bgeo_reduce(bgeo):
    """Reduce BGEO DB to trimer and residue scopes."""
    dres = {}
    dpairs = {}
    for btype in bgeo.keys():
        dres_ = dres.setdefault(btype, {})
        dpairs_ = dpairs.setdefault(btype, {})

        for res in bgeo[btype].keys():
            resangs = dres_.setdefault(res, [])
            dpairs__ = dpairs_.setdefault(res, {})

            for pairs in bgeo[btype][res].keys():
                respairs = dpairs__.setdefault(pairs, [])

                for tor in bgeo[btype][res][pairs].keys():
                    resangs.extend(bgeo[btype][res][pairs][tor])
                    respairs.extend(bgeo[btype][res][pairs][tor])

    return dpairs, dres
