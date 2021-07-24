"""Tools for conformer building operations."""
import itertools as it
import re
from collections import Counter, defaultdict, namedtuple
from functools import partial
from itertools import cycle

import numpy as np
from libfuncpy import flatlist, make_iterable, pass_
from numba import njit

# import idpcpp, imported locally at init_faspr_sidechains
from idpconfgen import log
from idpconfgen.core.build_definitions import (
    bonds_equal_3_inter,
    bonds_le_2_inter,
    distances_C_Np1,
    distances_CA_C,
    distances_N_CA,
    )
from idpconfgen.core.definitions import (
    bgeo_CaCNp1,
    bgeo_Cm1NCa,
    bgeo_NCaC,
    faspr_dun2010bbdep_path,
    )
from idpconfgen.libs.libcalc import (
    calc_all_vs_all_dists_njit,
    multiply_upper_diagonal_raw_njit,
    sum_upper_diagonal_raw_njit,
    )
from idpconfgen.libs.libenergyij import (
    energycalculator_ij,
    init_coulomb_calculator,
    init_lennard_jones_calculator,
    )
from idpconfgen.libs.libfilter import (
    aligndb,
    regex_forward_with_overlap,
    regex_search,
    )
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libparse import (
    get_mers,
    translate_seq_to_3l,
    remove_empty_keys,
    )
from idpconfgen.libs.libtimer import ProgressCounter, timeme


# See TODO at init_faspr_sidechains
# Variables related to the sidechain building process.
# # Variables for the FASPR algorithm.
# faspr_sc = idpcpp.faspr_sidechains
# faspr_dun2010_bbdep_str = str(faspr_dun2010bbdep_path)


ConfMasks = namedtuple(
    'ConfMaks',
    [
        'bb3',
        'bb4',
        'NHs',
        'COs',
        'Hterm',
        'OXT1',
        'OXT2',
        'cterm',
        'non_Hs',
        'non_Hs_non_OXT',
        'H1_N_CA_CB',
        ]
    )


ConfLabels = namedtuple(
    'ConfLabels',
    [
        'atom_labels',
        'res_nums',
        'res_labels',
        ]
    )
"""
Contain label information for a protein/conformer.

Named indexes
-------------
atom_labels
res_nums
res_labels
"""


def build_regex_substitutions(
        s,  # str
        options,  # dict
        pre_treatment=list,
        post_treatment=''.join,
        ):  # -> str
    """
    Build character replacements in regex string.

    Example
    -------
    >>> build_regex_substitutions('ASD', {'S': 'SE'})
    'A[SE]D'

    >>> build_regex_substitutions('ASDS', {'S': 'SE'})
    'A[SE]D[SE]'

    >>> build_regex_substitutions('ASDS', {})
    'ASDS'

    Parameters
    ----------
    s : regex string

    options : dict
        Dictionary of char to multichar substitutions

    pre_treatment : callable, optional
        A treatment to apply in `s` before substitution.
        `pre_treatment` must return a list-like object.
        Default: list because it expects s to be a string.

    post_treatment : callable, optional
        A function to apply on the resulting list-like object
        before returning.
        Default: ''.join, to return a string.
    """
    s_treat = pre_treatment(s)
    for i, char in enumerate(s_treat):
        if char in options:
            s_treat[i] = f'[{options[char]}]'

    return post_treatment(s_treat)


def init_confmasks(atom_labels):
    """
    Create a ConfMask object (namedtuple).

    ConfMask is a named tuple which attributes are integer masks for the
    respective groups.

    Parameters
    ----------
    atom_labels : array-like
        The atom names of the protein.

    Returns
    -------
    namedtuple
        ConfMasks object.

    Notes
    -----
    ConfMask attributes map to the following atom groups:

    bb3 : N, CA, C
    bb4 : N, CA, C, O
    NHs : amide protons
    Hterm : N-terminal protons
    OXT1 : O atom of C-terminal carboxyl group
    OXT2 : OXT atom of the C-terminal carboxyl group
    cterm : (OXT2, OXT1)
    non_Hs : all but hydrogens
    non_Hs_non_OXT : all but hydrogens and the only OXT atom
    H1_N_CA_CB : these four atoms from the first residue
                 if Gly, uses HA3.
    """
    bb3 = np.where(np.isin(atom_labels, ('N', 'CA', 'C')))[0]
    assert len(bb3) % 3 == 0

    bb4 = np.where(np.isin(atom_labels, ('N', 'CA', 'C', 'O')))[0]
    assert len(bb4) % 4 == 0

    NHs = np.where(atom_labels == 'H')[0]

    # last O belongs to the C-term carboxyl, we don't want it in the carbonyl
    # mask
    COs = np.where(atom_labels == 'O')[0][:-1]

    OXT1 = np.where(atom_labels == 'O')[0][-1]
    # the actual OXT atom, this is only one O of the C-term carboxyl pair
    OXT2 = np.where(atom_labels == 'OXT')[0][0]  # int instead of list

    cterm = [OXT2, OXT1]

    rr = re.compile(r'H+')
    hs_match = np.vectorize(lambda x: bool(rr.match(x)))
    non_Hs = np.where(np.logical_not(hs_match(atom_labels)))[0]
    non_Hs_non_OXT = non_Hs[:-1]

    # used to rotate the N-terminal Hs to -60 degrees to  HA during
    # the building process
    _H1_idx = np.where(atom_labels == 'H1')[0]
    assert len(_H1_idx) == 1

    # of the first residue
    _N_CA_idx = np.where(np.isin(atom_labels, ('N', 'CA')))[0][:2]
    assert len(_N_CA_idx) == 2, _N_CA_idx

    _final_idx = np.where(np.isin(atom_labels, ('CB', 'HA3')))[0][0:1]
    assert len(_final_idx) == 1

    H1_N_CA_CB = list(_H1_idx) + list(_N_CA_idx) + list(_final_idx)
    assert len(H1_N_CA_CB) == 4

    Hterm = np.where(np.isin(atom_labels, ('H1', 'H2', 'H3')))[0]
    assert len(Hterm) == 3

    conf_mask = ConfMasks(
        bb3=bb3,
        bb4=bb4,
        NHs=NHs,
        COs=COs,
        Hterm=Hterm,
        OXT1=OXT1,
        OXT2=OXT2,
        cterm=cterm,
        non_Hs=non_Hs,
        non_Hs_non_OXT=non_Hs_non_OXT,
        H1_N_CA_CB=H1_N_CA_CB,
        )

    return conf_mask


def init_conflabels(*args, **kwargs):
    """
    Create atom and residue labels from sequence.

    Parameters
    ----------
    *args, **kwargs
        Whichever `:func:create_conformer_labels` accepts.

    Returns
    -------
    namedtuple
        ConfLabels named tuple populated according to input sequence.

    See Also
    --------
    create_conformer_labels()
    ConfLabels
    """
    return ConfLabels(*create_conformer_labels(*args, **kwargs))


def create_conformer_labels(
        input_seq,
        atom_names_definition,
        transfunc=translate_seq_to_3l
        ):
    """
    Create all atom/residue labels model based on an input sequence.

    The labels are those expected for a all atom model PDB file. Hence,
    residue labels are repeated as needed in order to exist one residue
    label/number per atom.

    Parameters
    ----------
    input_seq : str
        The protein input sequence in 1-letter code format.

    atom_names_definition : dict
        Keys are residue identity and values are list/tuple of strings
        identifying atoms. Atom names should be sorted by the desired
        order.

    transfunc : func
        Function used to translate 1-letter input sequence to 3-letter
        sequence code.

    Returns
    -------
    tuple (atom labels, residue numbers, residue labels)
        Each is a np.ndarray of types: '<U4', np.int, and '<U3' and
        shape (N,) where N is the number of atoms.
        The three arrays have the same length.
    """
    input_seq_3_letters = transfunc(input_seq)
    # /
    # prepares data based on the input sequence
    # considers sidechain all-atoms
    atom_labels = np.array(
        make_list_atom_labels(
            input_seq,
            atom_names_definition,
            )
        )
    num_atoms = len(atom_labels)

    # /
    # per atom labels
    residue_numbers = np.empty(num_atoms, dtype=np.int)
    residue_labels = np.empty(num_atoms, dtype='<U3')

    # generators
    _res_nums_gen = gen_residue_number_per_atom(atom_labels, start=1)
    _res_labels_gen = \
        gen_3l_residue_labels_per_atom(input_seq_3_letters, atom_labels)

    # fills empty arrays from generators
    _zipit = zip(range(num_atoms), _res_nums_gen, _res_labels_gen)
    for _i, _num, _label in _zipit:
        residue_numbers[_i] = _num
        residue_labels[_i] = _label

    # maniatic cleaning from pre-function isolation
    del _res_labels_gen, _res_nums_gen, _zipit

    # ensure
    assert len(residue_numbers) == num_atoms
    assert len(residue_labels) == num_atoms, (len(residue_labels), num_atoms)
    # ?
    return atom_labels, residue_numbers, residue_labels


def make_list_atom_labels(input_seq, atom_labels_dictionary):
    """
    Make a list of the atom labels for an `input_seq`.

    Considers the N-terminal to be protonated H1 to H3.
    Adds also 'OXT' terminal label.

    Parameters
    ----------
    input_seq : str
        1-letter amino-acid sequence.

    atom_labels_dictionary : dict
        The ORDERED atom labels per residue.

    Returns
    -------
    list
        List of consecutive atom labels for the protein.
    """
    labels = []
    LE = labels.extend

    first_residue_atoms = atom_labels_dictionary[input_seq[0]]

    # the first residue is a special case, we add here the three protons
    # for consistency with the forcefield
    # TODO: parametrize? multiple forcefields?
    for atom in first_residue_atoms:
        if atom == 'H':
            LE(('H1', 'H2', 'H3'))
        else:
            labels.append(atom)

    for residue in input_seq[1:]:
        LE(atom_labels_dictionary[residue])

    labels.append('OXT')

    assert Counter(labels)['N'] == len(input_seq)
    assert labels[-1] == 'OXT'
    assert 'H1' in labels
    assert 'H2' in labels
    assert 'H3' in labels
    return labels


def gen_residue_number_per_atom(atom_labels, start=1):
    """
    Create a list of residue numbers based on atom labels.

    This is a contextualized function, not an abstracted one.
    Considers `N` to be the first atom of the residue.

    Yields
    ------
    ints
        The integer residue number per atom label.
    """
    assert atom_labels[0] == 'N', atom_labels[0]

    # creates a seamless interface between human and python 0-indexes
    start -= 1
    for al in atom_labels:
        if al == 'N':
            start += 1
        yield start


def gen_3l_residue_labels_per_atom(
        input_seq_3letter,
        atom_labels,
        ):
    """
    Generate residue 3-letter labels per atom.

    Parameters
    ----------
    input_seq_3letter : list of 3letter residue codes
        Most not be a generator.

    atom_labels : list or tuple of atom labels
        Most not be a generator.

    Yields
    ------
    String of length 3
         The 3-letter residue code per atom.
    """
    _count_labels = Counter(atom_labels)['N']
    _len = len(input_seq_3letter)
    assert _count_labels == _len, (_count_labels, _len)

    counter = -1
    for atom in atom_labels:
        if atom == 'N':
            counter += 1
        yield input_seq_3letter[counter]


def get_cycle_distances_backbone():
    """
    Return an inifinite iterator of backbone atom distances.

    Sampling, in order, distances between atom pairs:
        - N - Ca, used for OMEGA
        - Ca - C, used for PHI
        - C - N(+1), used for PSI
    """
    return cycle((
        distances_N_CA,  # used for OMEGA
        distances_CA_C,  # used for PHI
        distances_C_Np1,  # used for PSI
        ))


# deactivated after using bend library BGEO
# def get_cycle_bend_angles():
#     """
#     Return an infinite iterator of the bend angles.
#     """
#     return cycle((
#         build_bend_angles_Cm1_N_CA,  # used for OMEGA
#         build_bend_angles_N_CA_C,  # used for PHI
#         build_bend_angles_CA_C_Np1,  # used for PSI
#         ))


def get_cycle_bond_type():
    """
    Return an infinite interator of the bond types.

    Labels returns are synced with bgeo library.
    See `core.definitions.bgeo_*`.
    """
    return cycle((
        bgeo_Cm1NCa,  # used for OMEGA
        bgeo_NCaC,  # used for PHI
        bgeo_CaCNp1,  # used for PSI
        ))


#def read_db_to_slices(database, dssp_regexes, ncores=1):
#    """Create database base of slice and angles."""
#    # reads db dictionary from disk
#    db = read_dictionary_from_disk(database)
#    log.info(f'Read DB with {len(db)} entries')
#
#    # reads and prepares IDPConfGen data base
#    timed = partial(timeme, aligndb)
#    pdbs, angles, dssp, resseq = timed(db)
#
#    # searchs for slices in secondary structure, according to user requests
#    timed = partial(timeme, regex_search, ncores=ncores)
#    dssp_regexes = \
#        [dssp_regexes] if isinstance(dssp_regexes, str) else dssp_regexes
#    slices = []
#    for dssp_regex_string in dssp_regexes:
#        slices.extend(timed(dssp, dssp_regex_string))
#    log.info(f'Found {len(slices)} indexes for {dssp_regexes}')
#
#    return slices, angles


def read_db_to_slices_given_secondary_structure(database, ss_regexes):
    """
    Read slices in the DB that belong to a single secondary structure.

    Concatenates primary sequences accordingly.
    """
    log.info(f'ss_regexes, {ss_regexes}')
    db = read_dictionary_from_disk(database)
    timed = partial(timeme, aligndb)
    _, angles, dssp, resseq = timed(db)

    _ = (regex_search(dssp, _regex) for _regex in ss_regexes)
    slices = list(it.chain.from_iterable(_))
    seqs = [resseq[slc] for slc in slices]

    primary = '|'.join(seqs)

    omega, phi, psi = [], [], []

    oe, he, se = omega.extend, phi.extend, psi.extend
    oa, ha, sa = omega.append, phi.append, psi.append
    nan = np.nan

    for s in slices:
        oe(angles[s, 0])
        he(angles[s, 1])
        se(angles[s, 2])
        oa(nan)
        ha(nan)
        sa(nan)

    omega.pop()
    phi.pop()
    psi.pop()

    _omega = np.array(omega)
    _phi = np.array(phi)
    _psi = np.array(psi)

    seq_angles = np.array([_omega, _phi, _psi]).T

    assert seq_angles.shape[0] == len(primary)

    del omega, phi, psi
    # TODO:
    # save table to a file
    return primary, seq_angles


def prepare_slice_dict(
        primary,
        input_seq,
        mers_size=(1, 2, 3, 4, 5),
        res_tolerance=None,
        ncores=1,
        ):
    """
    Prepare a dictionary mapping chunks to slices in `primary`.

    Parameters
    ----------
    primary : str
        A concatenated version of all primary sequences in the database.
        In the form of "QWERY|IPASDF", etc.

    input_seq : str
        The 1-letter code amino-acid sequence of the conformer to construct.

    mers_size : iterable
        A iterable of integers denoting the size of the chunks to search
        for. Defaults from 1 to 5.

    res_tolerance : dict
        A dictionary mapping residue tolerances, for example:
        {"A": "AIL"}, noting Ala can be replaced by Ile and Leu in the
        search (this is a dummy example).

    ncores : int
        The number of processors to use.

    Return
    ------
    dict
        A dict with the given mapping. First key-leve of the dict
        is the length of the chunks, hence, integers.
        The second key level are the residue chunks found in the `primary`.
        A chunk in input_seq but not in `primary` is removed from the
        dict.
    """
    res_tolerance = res_tolerance or {}
    mers_size = make_iterable(mers_size)

    log.info('preparing regex xmers')

    xmers = (get_mers(input_seq, i) for i in mers_size)
    xmers_flat = flatlist(xmers)
    slice_dict = defaultdict(dict)

    with ProgressCounter(suffix='Searching for xmers: ') as PW:
        for mer in xmers_flat:
            lmer = len(mer)
            altered_mer = build_regex_substitutions(mer, res_tolerance)

            merregex = f'(?=({altered_mer}))'

            slice_dict[lmer][altered_mer] = \
                regex_forward_with_overlap(primary, merregex)

            # if no slices were found
            if not slice_dict[lmer][altered_mer]:
                slice_dict[lmer].pop(altered_mer)

            # this is a trick to find the sequence that are proceeded
            # by Proline. Still, these are registered in the "lmer" size
            # without consider the proline addition.
            # this is a combo with get_adjancent_angles
            merregex_P = f'(?=({altered_mer}P))'
            altered_mer_P = f'{altered_mer}_P'
            slice_dict[lmer][altered_mer_P] = \
                regex_forward_with_overlap(primary, merregex_P)

            # if no entrey was found
            if not slice_dict[lmer][altered_mer_P]:
                slice_dict[lmer].pop(altered_mer_P)
            # for _s in slice_dict[lmer][mer]:
                # assert '|' not in input_seq[_s]
            PW.increment()

    return slice_dict


# Other functions should have the same API:
# parameters = input_seq
def init_faspr_sidechains(
        input_seq,
        # faspr_dun2010db_spath=faspr_dun2010_bbdep_str,
        # faspr_func=faspr_sc,
        ):
    """
    Instantiate dedicated function environment for FASPR sidehchain calculation.

    Examples
    --------
    >>> calc_faspr = init_fastpr_sidechains('MASFRTPKKLCVAGG')
    >>> # a (N, 3) array with the N,CA,C,O coordinates
    >>> coords = np.array( ... )
    >>> calc_faspr(coords)

    Parameters
    ----------
    input_seq : str
        The FASTA sequence of the protein for which this function will
        be used.

    Returns
    -------
    np.ndarray (M, 3)
        Heavy atom coordinates of the protein sequence.
    """
    # TODO:
    # this is here because tox is not able to detect idpcpp module.
    # this is a turnaround to allow tests to pass.
    # currently tests to not test this function.
    import idpcpp
    faspr_func = idpcpp.faspr_sidechains
    faspr_dun2010_bbdep_str = str(faspr_dun2010bbdep_path)

    def compute_faspr_sidechains(coords):
        """Do calculation."""
        return faspr_func(coords, input_seq, faspr_dun2010_bbdep_str)

    return compute_faspr_sidechains


def prepare_energy_function(
        atom_labels,
        residue_numbers,
        residue_labels,
        forcefield,
        lj_term=True,
        coulomb_term=False,
        energy_threshold_type='nansum',
        **kwnull,
        ):
    """."""
    # this mask identifies covalently bonded pairs and pairs two bonds apart
    bonds_le_2_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        forcefield.bonds_le2_intra,
        bonds_le_2_inter,
        base_bool=False,
        )

    # this mask identifies pairs exactly 3 bonds apart
    bonds_exact_3_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        forcefield.bonds_eq3_intra,
        bonds_equal_3_inter,
        )

    # /
    # assemble energy function
    energy_func_terms = []

    if lj_term:

        acoeff, bcoeff = create_LJ_params_raw(
            atom_labels,
            residue_numbers,
            residue_labels,
            forcefield.forcefield,
            )

        # 0.2 as 0.4
        _lj14scale = float(forcefield.forcefield['lj14scale'])
        acoeff[bonds_exact_3_mask] *= _lj14scale * 0.2
        bcoeff[bonds_exact_3_mask] *= _lj14scale * 0.2
        acoeff[bonds_le_2_mask] = np.nan
        bcoeff[bonds_le_2_mask] = np.nan

        lf_calc = init_lennard_jones_calculator(
            acoeff,
            bcoeff,
            postf=energy_threshold_type,
            )
        energy_func_terms.append(lf_calc)
        log.info('prepared lj')

    if coulomb_term:

        charges_ij = create_Coulomb_params_raw(
            atom_labels,
            residue_numbers,
            residue_labels,
            forcefield.forcefield,
            )

        charges_ij[bonds_exact_3_mask] *= float(forcefield.forcefield['coulomb14scale'])  # noqa: E501
        charges_ij[bonds_le_2_mask] = np.nan

        coulomb_calc = init_coulomb_calculator(charges_ij, postf=energy_threshold_type)
        energy_func_terms.append(coulomb_calc)
        log.info('prepared Coulomb')

    # in case there are iji terms, I need to add here another layer
    calc_energy = energycalculator_ij(
        calc_all_vs_all_dists_njit,
        energy_func_terms,
        )
    log.info('done preparing energy func')
    return calc_energy


def create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        bonds_intra,
        bonds_inter,
        base_bool=False,
        ):
    """
    Create bool mask array identifying the pairs X bonds apart in ij pairs.

    Given `bonds_intra` and `bonds_inter` criteria, idenfities those ij
    atom pairs in N*(N-1)/2 condition (upper all vs all diagonal) that
    agree with the described bonds.

    Inter residue bonds are only considered for consecutive residues.

    Paramters
    ---------
    atom_labels : iterable, list or np.ndarray
        The protein atom labels. Ex: ['N', 'CA, 'C', 'O', 'CB', ...]

    residue_numbers : iterable, list or np.ndarray
        The protein residue numbers per atom in `atom_labels`.
        Ex: [1, 1, 1, 1, 1, 2, 2, 2, 2, ...]

    residue_labels : iterable, list or np.ndarray
        The protein residue labels per atom in `atom_labels`.
        Ex: ['Met', 'Met', 'Met', ...]

    Depends
    -------
    `gen_ij_pairs_upper_diagonal`
    `gen_atom_pair_connectivity_masks`
    """
    atom_labels_ij_gen = gen_ij_pairs_upper_diagonal(atom_labels)
    residue_numbers_ij_gen = gen_ij_pairs_upper_diagonal(residue_numbers)
    residue_labels_ij_gen = gen_ij_pairs_upper_diagonal(residue_labels)

    bonds_indexes_gen = gen_atom_pair_connectivity_masks(
        residue_labels_ij_gen,
        residue_numbers_ij_gen,
        atom_labels_ij_gen,
        bonds_intra,
        bonds_inter,
        )

    num_ij_pairs = len(atom_labels) * (len(atom_labels) - 1) // 2
    other_bool = not base_bool
    bonds_mask = np.full(num_ij_pairs, base_bool)

    for idx in bonds_indexes_gen:
        bonds_mask[idx] = other_bool

    return bonds_mask


def gen_ij_pairs_upper_diagonal(data):
    """
    Generate upper diagonal ij pairs in tuples.

    The diagonal is not considered.

    Yields
    ------
    tuple of length 2
        IJ pairs in the form of N*(N-1) / 2.
    """
    for i in range(len(data) - 1):
        for j in range(i + 1, len(data)):
            yield (data[i], data[j])


def gen_atom_pair_connectivity_masks(
        res_names_ij,
        res_num_ij,
        atom_names_ij,
        connectivity_intra,
        connectivity_inter,
        ):
    """
    Generate atom pair connectivity indexes.

    Given atom information for the ij pairs and connectivity criteria,
    yields the index of the ij pair if the pair is connected according
    to the connectivity criteria.

    For example, if the ij pair is covalently bonded, or 3 bonds apart,
    etc.

    Parameters
    ----------
    res_names_ij
    res_num_ij,
    atom_names_ij, iterables of the same length and synchronized information.

    connectivity_intra,
    connectivity_inter, dictionaries mapping atom labels connectivity

    Depends
    -------
    `are_connected`
    """
    zipit = zip(res_names_ij, res_num_ij, atom_names_ij)
    counter = 0
    for (rn1, _), (n1, n2), (a1, a2) in zipit:

        found_connectivity = are_connected(
            int(n1),
            int(n2),
            rn1,
            a1,
            a2,
            connectivity_intra,
            connectivity_inter,
            )

        if found_connectivity:
            yield counter

        counter += 1


def create_LJ_params_raw(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        ):
    """Create ACOEFF and BCOEFF parameters."""
    sigmas_ii = extract_ff_params_for_seq(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        'sigma',
        )

    epsilons_ii = extract_ff_params_for_seq(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        'epsilon',
        )

    num_ij_pairs = len(atom_labels) * (len(atom_labels) - 1) // 2
    # sigmas
    sigmas_ij_pre = np.empty(num_ij_pairs, dtype=np.float64)
    sum_upper_diagonal_raw_njit(np.array(sigmas_ii), sigmas_ij_pre)
    #
    # epsilons
    epsilons_ij_pre = np.empty(num_ij_pairs, dtype=np.float64)
    multiply_upper_diagonal_raw_njit(
        np.array(epsilons_ii),
        epsilons_ij_pre,
        )
    #

    # mixing rules
    epsilons_ij = epsilons_ij_pre ** 0.5
    # mixing + nm to Angstrom converstion
    # / 2 and * 10
    sigmas_ij = sigmas_ij_pre * 5

    acoeff = 4 * epsilons_ij * (sigmas_ij ** 12)
    bcoeff = 4 * epsilons_ij * (sigmas_ij ** 6)

    return acoeff, bcoeff


def create_Coulomb_params_raw(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        ):
    """."""
    charges_i = extract_ff_params_for_seq(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        'charge',
        )

    num_ij_pairs = len(atom_labels) * (len(atom_labels) - 1) // 2
    charges_ij = np.empty(num_ij_pairs, dtype=np.float64)
    multiply_upper_diagonal_raw_njit(charges_i, charges_ij)
    charges_ij *= 0.25  # dielectic constant

    return charges_ij


def extract_ff_params_for_seq(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        param,
        ):
    """
    Extract a parameter from forcefield dictionary for a given sequence.

    See Also
    --------
    create_conformer_labels

    Parameters
    ----------
    atom_labels, residue_numbers, residue_labels
        As returned by `:func:create_conformer_labels`.

    forcefield : dict

    param : str
        The param to extract from forcefield dictionary.
    """
    params_l = []
    params_append = params_l.append

    zipit = zip(atom_labels, residue_numbers, residue_labels)
    for atom_name, res_num, res_label in zipit:

        # adds C to the terminal residues
        if res_num == residue_numbers[-1]:
            res = 'C' + res_label
            was_in_C_terminal = True
            assert res.isupper() and len(res) == 4, res

        elif res_num == residue_numbers[0]:
            res = 'N' + res_label
            was_in_N_terminal = True
            assert res.isupper() and len(res) == 4, res

        else:
            res = res_label

        # TODO:
        # define protonation state in parameters
        if res_label.endswith('HIS'):
            res_label = res_label[:-3] + 'HIP'

        try:
            # force field atom type
            atype = force_field[res][atom_name]['type']

        # TODO:
        # try/catch is here to avoid problems with His...
        # for this purpose we are only using side-chains
        except KeyError:
            raise KeyError(tuple(force_field[res].keys()))

        params_append(float(force_field[atype][param]))

    assert was_in_C_terminal, \
        'The C terminal residue was never computed. It should have.'
    assert was_in_N_terminal, \
        'The N terminal residue was never computed. It should have.'

    assert isinstance(params_l, list)
    return params_l


def are_connected(n1, n2, rn1, a1, a2, bonds_intra, bonds_inter):
    """
    Detect if a certain atom pair is bonded accordind to criteria.

    Considers only to the self residue and next residue
    """
    # requires
    assert isinstance(n1, int) and isinstance(n2, int), (type(n1), type(n2))
    assert all(isinstance(i, str) for i in (rn1, a1, a2)), \
        (type(i) for i in (rn1, a1, a2))
    assert all(isinstance(i, dict) for i in (bonds_intra, bonds_inter)), \
        (type(i) for i in (bonds_intra, bonds_inter))

    answer = (
        (n1 == n2 and a2 in bonds_intra[rn1][a1])
        or (
            n1 + 1 == n2
            and (
                a1 in bonds_inter  # this void KeyError
                and a2 in bonds_inter[a1]
                )
            )
        )

    assert isinstance(answer, bool)
    return answer


def create_sidechains_masks_per_residue(
        residue_numbers,
        atom_labels,
        backbone_atoms,
        ):
    """
    Create a map of numeric indexing masks pointing to side chains atoms.

    Create separate masks per residue.

    Parameters
    ----------
    residue_numbers : np.ndarray, shape (N,)
        The atom residue numbers of the protein.

    atom_labels : np.ndarray, shape (N,)
        The atom labels of the protein.

    backbone_atoms : list or tuple
        The labels of all possible backbone atoms.

    Returns
    -------
    list of tuples of length 2
        List indexes refer to protein residues, index 0 is residue 1.
        Per residue, a tuple of length 2 is given. Tuple index 0 are
        the indexes of that residue sidechain atoms mapped to an array
        of the `atom_labels` and `residue_numbers` characteristics.
        The tuple index 1 is an array of length M, where M is the number
        of sidechain atoms for that residue, defaults to np.nan.
    """
    assert residue_numbers.size == atom_labels.size
    assert type(backbone_atoms) in (list, tuple, np.ndarray)

    ss = []
    ssa = ss.append

    is_backbone = np.isin(atom_labels, backbone_atoms)
    is_sidechain = np.logical_not(is_backbone)

    float64 = np.float64
    full = np.full
    logical_and = np.logical_and
    nan = np.nan
    where = np.where

    # sorted here is mandatory because ss indexes must follow the residue
    # numbering
    for resnum in sorted(list(set(residue_numbers))):
        is_residue = residue_numbers == resnum
        bool_result = logical_and(is_residue, is_sidechain)
        ssa((
            where(bool_result)[0],
            full((np.sum(bool_result), 3), nan, dtype=float64)
            ))

    return ss


# njit available
def get_indexes_from_primer_length(
        sequence,
        plen,
        current_residue,
        ):
    """Get sequence chunk based on position and length."""
    if plen == 1:
        return current_residue
    elif plen == 2:
        return sequence[current_residue: current_residue + 2]
    elif plen == 3:
        return sequence[current_residue - 1: current_residue + 3]
    elif plen == 4:
        return sequence[current_residue - 1: current_residue + 4]
    elif plen == 5:
        return sequence[current_residue - 2: current_residue + 5]
    elif plen == 6:
        return sequence[current_residue - 2: current_residue + 6]
    elif plen == 7:
        return sequence[current_residue - 3: current_residue + 7]


compute_sidechains = {
    'faspr': init_faspr_sidechains,
    }

get_idx_primer_njit = njit(get_indexes_from_primer_length)
