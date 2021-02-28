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
        'non_Hs',
        'non_Hs_non_OXT',
        'H1_N_CA_CB',
        ]
    )


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


compute_sidechains = {
    'faspr': init_faspr_sidechains,
    }


