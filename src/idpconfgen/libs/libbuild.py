"""Tools for conformer building operations."""
import re
from collections import namedtuple
from itertools import cycle

import numpy as np

import idpcpp

from idpconfgen.core.definitions import faspr_dun2010bbdep_path
from idpconfgen.core.build_definitions import (
    distances_CA_C,
    distances_C_Np1,
    distances_N_CA,
    )


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


def init_conflabels(*args, **kwargs):
    return ConfLabels(create_conformer_labels(*args, **kwargs))



# This is a class because the moment of instantiation is thought to be
# different from the moment on which the attributes are populated
class ConformerLabels:

    __slots__ = [
        'atom_labels',
        'res_nums',
        'res_labels',
        'atom_masks',
        'num_atoms',
        'energy_func',
        '_al',
        ]

    def __init__(self, atom_names):

        self.atom_labels = None
        self.res_nums = None
        self.res_labels = None
        self.atom_masks = None
        self.num_atoms = None
        self._al = atom_names

        return

    def define_labels(self, input_seq):
        self.atom_labels, self.res_nums, self.res_labels = \
            create_conformer_labels(input_seq, topology.atom_labels)

        self.num_atoms = len(self.atoms_labels)
        self.atom_masks = init_confmasks(self.atom_labels)
        del _al



# TODO: correct for HIS/HIE/HID/HIP
def translate_seq_to_3l(input_seq):
    return [
        'HIP' if _res == 'H' else aa1to3[_res]
        for _res in input_seq
        ]


def create_conformer_labels(
        input_seq,
        atom_labels_definition,
        transfunc=translate_seq_to_3l):
    """
    Create all atom labels model based on an input sequence.

    The use of `input_seq_3_letters` is still experimental. Is likely
    to be refactored.
    """
    input_seq_3_letters = transfunc(input_seq)
    # /
    # prepares data based on the input sequence
    # considers sidechain all-atoms
    atom_labels = np.array(
        make_list_atom_labels(
            input_seq,
            atom_labels_definition,
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


#def create_conf_labels(input_seq):
#    """
#    Populate the global variables for all-atom representation labels.
#
#    Parameters
#    ----------
#    input_seq : str
#        The input sequence of the protein to be built.
#
#    Returns
#    -------
#    None
#        Global variables as changed in place.
#    """
#    input_seq_3l = translate_seq_to_3l(input_seq)
#    atom_labels, res_nums, res_labels = create_conformer_labels(
#        input_seq,
#        input_seq_3l,
#        )
#    atom_masks = init_confmasks(atom_labels)
#
#    return (atom_labels, res_nums, res_labels, atom_masks)



# Variables related to the sidechain building process.
# # Variables for the FASPR algorithm.
faspr_sc = idpcpp.faspr_sidechains


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


def get_cycle_distances_backbone():
    """
    Return an inifinite iterator of backbone atom distances.
    """
    return cycle((
            distances_N_CA,  # used for OMEGA
            distances_CA_C,  # used for PHI
            distances_C_Np1,  # used for PSI
            ))


# deactivated after using bend library BGEO
#def get_cycle_bend_angles():
#    """
#    Return an infinite iterator of the bend angles.
#    """
#    return cycle((
#        build_bend_angles_Cm1_N_CA,  # used for OMEGA
#        build_bend_angles_N_CA_C,  # used for PHI
#        build_bend_angles_CA_C_Np1,  # used for PSI
#        ))


def get_cycle_bond_type():
    """Return an infinite interator of the bond types."""
    return cycle((
        'Cm1_N_Ca',
        'N_Ca_C',
        'Ca_C_Np1',
        ))


# Other functions should have the same APIs:
# parameters = input_seq
def init_faspr_sidechains(
        input_seq,
        faspr_dun2010db_spath=str(faspr_dun2010bbdep_path),
        faspr_func=faspr_sc,
        ):
    """
    Initiates dedicated function environment for FASPR sidehchain calculation.

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
    def compute_faspr_sidechains(coords):
        """Does calculation."""
        return faspr_func(coords, input_seq, faspr_dun2010db_spath)

    return compute_faspr_sidechains


def prepare_energy_function(
        atom_labels,
        residue_numbers,
        residue_labels,
        topology,
        lj_term=True,
        coulomb_term=False,
        **kwnull,
        ):
    # /
    # Prepare Topology
    #ff14SB = read_ff14SB_params()
    #ambertopology = AmberTopology(add_OXT=True, add_Nterminal_H=True)

    # this mask will deactivate calculations in covalently bond atoms and
    # atoms separated 2 bonds apart
    bonds_le_2_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        topology.bonds_le2_intra,
        bonds_le_2_inter,
        base_bool=False,
        )

    bonds_exact_3_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        topology.bonds_eq3_intra,
        bonds_equal_3_inter,
        )

    # /
    # assemble energy function
    energy_func_terms = []
    if lj_term:

        ap_acoeff, ap_bcoeff = create_LJ_params_raw(
            atom_labels,
            residue_numbers,
            residue_labels,
            topology.forcefield,
            )

        ap_acoeff[bonds_exact_3_mask] *= float(ff14SB['lj14scale']) * 0.2  # was 0.4
        ap_bcoeff[bonds_exact_3_mask] *= float(ff14SB['lj14scale']) * 0.2
        ap_acoeff[bonds_le_2_mask] = np.nan
        ap_bcoeff[bonds_le_2_mask] = np.nan

        lf_calc = init_lennard_jones_calculator(ap_acoeff, ap_bcoeff)
        energy_func_terms.append(lf_calc)
        print('prepared lj')

    if coulomb_term:

        ap_charges_ij = create_Coulomb_params_raw(
            atom_labels,
            residue_numbers,
            residue_labels,
            topology.forcefield,
            )

        ap_charges_ij[bonds_exact_3_mask] *= float(ff14SB['coulomb14scale'])
        ap_charges_ij[bonds_le_2_mask] = np.nan

        coulomb_calc = init_coulomb_calculator(charges_ij)
        energy_func_term.append(coulomb_calc)
        print('prepared Coulomb')

    # in case there are iji terms, I need to add here another layer
    calc_energy = energycalculator_ij(
        calc_all_vs_all_dists,
        energy_func_terms,
        )
    print('done preparing energy func')
    return calc_energy
    # ?


compute_sidechains = {
    'faspr': init_faspr_sidechains,
    }


