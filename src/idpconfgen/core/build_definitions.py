"""Definitions for the building process."""
import xml.etree.ElementTree as ET
from collections import defaultdict
from copy import copy
from math import pi
from pathlib import Path
from pprint import pprint
from statistics import mean, stdev

import numpy as np
from scipy import spatial

from idpconfgen.libs.libstructure import Structure, col_name
from idpconfgen.core.definitions import aa1to3, aa3to1


pdist = spatial.distance.pdist
_filepath = Path(__file__).resolve().parent  # folder
_sidechain_template_files = sorted(list(
    _filepath.joinpath('sidechain_templates', 'pdb_names').glob('*.pdb')))
amber_pdbs = sorted(list(
    _filepath.joinpath('sidechain_templates', 'amber_names').glob('*.pdb')))
_amber14sb = _filepath.joinpath('data', 'protein.ff14SB.xml')

backbone_atoms = ('N', 'C', 'CA', 'O', 'OXT', 'H', 'H1', 'H2', 'H3')

# amino-acids atom labels
# from: http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
# PDB column
# Taken from PDB entry 6I1B REVDAT 15-OCT-92.


def _read_labels(pdbs):
    """Read atom labels from residue template files."""
    labels_d = {}
    for pdb in pdbs:
        s = Structure(pdb)
        s.build()
        a_labels = tuple(s.data_array[:, col_name])
        pdb_name = pdb.stem.upper()
        pdb_1letter = aa3to1[pdb_name]

        # labels are duplicated for 1-letter and 3-letter codes to avoid
        # double dictionary lookup in other implementations
        labels_d[pdb_name] = a_labels
        labels_d[pdb_1letter] = a_labels
    return labels_d


# support figure, for the different histidine protonation states.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3639364/figure/Fig1/
atom_names_pdb = _read_labels(_sidechain_template_files)
atom_names_amber = _read_labels(amber_pdbs)

def read_ff14SB_params():
    """
    Read Amber protein.ff14SB.xml parameters to a dictionary.

    `protein.ff14SB.xml` is in `src.idpconfgen.core.data` folder.

    Dictionary structure:

    dict_keys(
        ['protein-C',
        'protein-CA',
        (... atom types ...)
        'protein-3C',
        'protein-C8',
        'ALA',
        'ARG',
        (... residues ...)
        'NTYR',
        'NVAL',
        'coulomb14scale',   <- note these two keys
        'lj14scale'])

    Atoms types have the structure:

    'protein-OH': {
        'class': 'OH',
        'element': 'O',
        'epsilon': '0.8803136',
        'mass': '16.0',
        'sigma': '0.3066473387839048',
        }

    Residues have the structure:

    'VAL': {
        'C': {'charge': '0.5973', 'type': 'protein-C'},
        'CA': {'charge': '-0.0875', 'type': 'protein-CX'},
        'CB': {'charge': '0.2985', 'type': 'protein-3C'},
        'CG1': {'charge': '-0.3192', 'type': 'protein-CT'},
        'CG2': {'charge': '-0.3192', 'type': 'protein-CT'},
        'H': {'charge': '0.2719', 'type': 'protein-H'},
        (...)
        }

    Returns
    -------
    dict
    """

    with open(_amber14sb, 'r') as fin:
        ff14sb = ET.fromstring(fin.read())

    ff14SB_params = defaultdict(dict)

    for child in ff14sb:
        if child.tag == 'AtomTypes':
            for atomtype in child:
                key = atomtype.attrib['name']
                ff14SB_params[key].update(atomtype.attrib)
                ff14SB_params[key].pop('name')

        elif child.tag == 'Residues':
            for residue in child:
                for atom in filter(lambda x: x.tag == 'Atom', residue):
                    key = residue.attrib['name']
                    atom_name = atom.attrib['name']

                    atom_par = ff14SB_params[key].setdefault(atom_name, {})
                    atom_par.update(atom.attrib)
                    ff14SB_params[key][atom_name].pop('name')

        elif child.tag == 'NonbondedForce':
            ff14SB_params.update(child.attrib)
            for atom in child:
                if atom.tag == 'Atom':
                    key = atom.attrib['type']
                    ff14SB_params[key].update(atom.attrib)
                    ff14SB_params[key].pop('type')

    return ff14SB_params


def generate_residue_template_topology(
        pdb_files,
        residue_labels,
        add_OXT=True,
        add_Nterminal_H=True,
        ):
    """
    Generate topology for the residue templates.

    Residue templates are stored in folder core/sidechain_templates/

    Generated topology represents all vs. all. Meaning, if atom A is
    covalently bond to atom B, B appears in A records and A in B records,
    as well.

    Returns
    -------
    dict
        A dictionary with the topology.
    """
    # Creates topoloty of amino acids templates, that is, a dictionary of
    # all vs. all covalent bond pairs
    res_covalent_bonds = defaultdict(dict)
    for pdb in pdb_files:

        pdbname = pdb.stem.upper()

        s = Structure(pdb)
        s.build()
        coords = s.coords
        atoms = s.data_array[:, 2]

        # atom names of residue templates must be sorted equally to the
        # atom_names dictionary, or vice-versa.
        _list_names = tuple(s.data_array[:, col_name])
        _atoms = residue_labels[pdbname]
        assert _list_names == _atoms, (
            'Atom names in `atom_names` dictionary differ from the '
            f'atom names in {pdbname} residue template.'
            )

        all_dists = pdist(coords)
        # note that 1.6 AA wont capture the S bonds which are 1.8 A away
        # the Cys and Met special cases are treated after the loop
        cov_bonds = all_dists <= 1.6

        # all vs. all atom pairs
        atom_pairs = (
            (a, b)
            for i, a in enumerate(atoms, start=1)
            for b in atoms[i:]
            )

        current_pdb = res_covalent_bonds[pdbname]
        for (a, b), is_covalent in zip(atom_pairs, cov_bonds):
            if is_covalent:
                connects_a = current_pdb.setdefault(a, [])
                connects_a.append(b)
                connects_b = current_pdb.setdefault(b, [])
                connects_b.append(a)

    # special cases: CYS and MET
    res_covalent_bonds['CYS']['CB'].append('SG')
    res_covalent_bonds['CYS']['SG'] = []
    res_covalent_bonds['CYS']['SG'].extend(('CB', 'HG'))
    res_covalent_bonds['CYS']['HG'] = []
    res_covalent_bonds['CYS']['HG'].append('SG')

    res_covalent_bonds['MET']['CG'].append('SD')
    res_covalent_bonds['MET']['SD'] = []
    res_covalent_bonds['MET']['SD'].extend(('CG', 'CE'))
    res_covalent_bonds['MET']['CE'].append('SD')

    # asserts all atoms are considered
    for k1, v1 in res_covalent_bonds.items():
        assert len(v1) == len(residue_labels[k1]), k1

        # add OXT connectivity
        if add_OXT:
            add_OXT_to_connectivity(v1)

        ## added 'OXT' connectivity
        #for atom, connects in v1.items():
        #    if 'O' in connects:
        #        connects.append('OXT')

        ## this should be only 'C'
        #v1['OXT'] = copy(v1['O'])
        if k1 == 'PRO':
            continue

        if add_Nterminal_H:
            add_Nterm_H_connectivity(v1)

    return res_covalent_bonds


def add_Nterm_H_connectivity(connectivity_dict):
    """
    Adds protons for Nterm connectivity.

    Adds H1, H2, and H3 protons to N connectivity.
    This maintains compatibility with XML forcefields obtained
    from the OpenMM project.
    """
    assert 'H' in connectivity_dict, connectivity_dict
    for atom, list_of_connects in connectivity_dict.items():
        if 'H' in list_of_connects:
            list_of_connects.extend(('H1', 'H2', 'H3'))

    for h in ('H1', 'H2', 'H3'):
        connectivity_dict[h] = copy(connectivity_dict['H'])


def add_OXT_to_connectivity(connectivity_dict):
    """Add OXT connectivity to residue."""
    connectivity_dict['OXT'] = copy(connectivity_dict['O'])
    for _atom in connectivity_dict['OXT']:
        connectivity_dict[_atom].append('OXT')

    assert 'OXT' in connectivity_dict
    assert set(connectivity_dict['OXT']).issubset(connectivity_dict.keys())
    assert connectivity_dict['OXT'] == connectivity_dict['O']


def _recursive_bonds_apart(
        init_cov_res,
        bonds_apart,
        cov_bonded,
        counter,
        max_bonds,
        ):
    """
    Expand the list of bonded atoms recursively.

    Parameters
    ----------
    init_cov_res : dict
        A dictionary representing a residues (or structure) covalent
        bonds for each atom of the structure.

    bonds_apart : list
        A list that will be filled with the new atoms according to the
        `max_bond` criteria. Normally this list contains the atoms
        already covalently bond to the atom being processed.
        `bonds_apart` list is modified in place.

    cov_bonded : list
        The initial atoms covalentely bond to the atom being processed.

    max_bonds : int
        The maximum number of bonds to consider.

    Returns
    -------
    None
        Modifies `bonds_apart` in place.
    """
    if counter > max_bonds:
        return

    for atom in cov_bonded:
        cov_bonded_next = copy(init_cov_res[atom])
        bonds_apart.extend(set(cov_bonded_next).difference(set(bonds_apart)))
        _recursive_bonds_apart(
            init_cov_res,
            bonds_apart,
            cov_bonded_next,
            counter + 1,
            max_bonds)


def expand_topology_bonds_apart(cov_bond_dict, bonds_apart):
    """
    Expand a topogy dictionary of covalent bonds to X bonds apart.

    The resulting topology is defined by the specified number of bonds apart.

    Parameters
    ----------
    cov_bond_dict : dict
        A dictionary with covalent bonds topology of a structure.
        The dictionary must be all vs. all with repetitions, meaing,
        if A is covalently bond to B, B should be noted in A and A should
        be noted in B, as well.
        `cov_bond_dict` can be generated previously with
        :func:`generate_residue_template_topology`.

    bonds_apart : int
        The number of bonds apart to consider expanding.

    Returns
    -------
    dict
        The expanded topology dictionary, all vs. all with repetition.
    """
    # the size of the data here is small so I have prioritize readability
    # and sanity check instead of speed
    expanded_topology = {}

    for res in cov_bond_dict:

        res_d = expanded_topology.setdefault(res, {})

        for atom in cov_bond_dict[res]:

            atoms_X_bonds_apart = res_d.setdefault(atom, [])

            cov_bonded = copy(cov_bond_dict[res][atom])

            atoms_X_bonds_apart.extend(cov_bonded)

            _recursive_bonds_apart(
                cov_bond_dict[res],
                atoms_X_bonds_apart,
                cov_bonded,
                2,  # starts at 2 cause the 1st iteration is the 2nd bond
                bonds_apart,
                )

            # the self atom is added to the list unavoidable
            atoms_X_bonds_apart.remove(atom)

    return expanded_topology



def topology_3_bonds_apart(covalent_bond_dict):
    """
    Map atom connectivity EXACTLY 3 bonds apart.

    See Amber20 manual Figure 14.1.

    Parameters
    ----------
    covalent_bond_dict : dict
        Per residue covalent connectivity.

    Returns
    -------
    dict
        Per residue 3 bonds apart connectivity.
    """
    x_bonds_apart = {}
    prevs = set()

    for res, residue_atoms in covalent_bond_dict.items():
        res_d = x_bonds_apart.setdefault(res, {})

        # for all atoms in the residue
        for atom in residue_atoms:
            xba = res_d.setdefault(atom, set())
            prevs.clear()
            prevs.add(atom)

            # for CA, H in N... one bond
            for subatom1 in set(residue_atoms[atom]).difference(prevs):
                prevs.add(subatom1)

                # for CB, HA in CA... two bonds
                for subatom2 in set(residue_atoms[subatom1]).difference(prevs):
                    prevs.add(subatom2)

                    # for HB1, HB2 in CB... three bonds
                    for subatom3 in set(residue_atoms[subatom2]).difference(prevs):
                        xba.add(subatom3)

    return x_bonds_apart


# interresidue exact 3 bonds connectivity
bonds_equal_1_inter = {'C': ['N']}

bonds_equal_3_inter = {
    'N': ['N'],
    'CA': ['H', 'CA'],
    'C': ['HA', 'HA2', 'HA3', 'CB', 'C'],
    'O': ['CA', 'H'],
    'HA': ['N'],
    'HA2': ['N'],
    'HA3': ['N'],
    'CB': ['N'],
    }

# /
# Less than or equal to:
# mar 15 dic 2020 15:50:48 EST

bonds_le_2_inter = {
    'C': ['CA', 'H', 'N'],
    'CA': ['N'],
    'O': ['N'],
    }

# interresidue 3-bonds connectivity
_C_3_connectivities = [
    'N',  # 1 bond apart
    'H', 'CA',  # 2 bonds apart
    'HA', 'HA2', 'HA3', 'CB', 'C',  # 3 bonds apart
    ]
_CA_O_3_connectivities = ['N', 'H', 'CA']
_N_HA_HA2_HA3_CB_3_connectivities = ['N']

# less or equal to:
bonds_le_3_inter = {
    'C': _C_3_connectivities,
    #
    'N': _N_HA_HA2_HA3_CB_3_connectivities,
    'HA': _N_HA_HA2_HA3_CB_3_connectivities,
    'HA2': _N_HA_HA2_HA3_CB_3_connectivities,
    'HA3': _N_HA_HA2_HA3_CB_3_connectivities,
    'CB': _N_HA_HA2_HA3_CB_3_connectivities,
    #
    'CA': _CA_O_3_connectivities,
    'O': _CA_O_3_connectivities,
    #
    }

# interresidue 4-bonds connectivity
#C_4_connectivities = [
#    'N',  # 1 bond apart
#    'H', 'CA',  # 2 bond apart
#    'HA', 'CB', 'C',  # 3 bond apart
#    'O', '1HB', '2HB', '3HB', 'CG',  # 4 bond apart
#    ]
#CA_O_4_connectivities = ['N', 'H', 'CA', 'HA', 'C', 'CB']
#N_HA_CB_4_connectivities = ['N', 'H', 'CA']

#    'C': C_4_connectivities,
#    #
#    'N': N_HA_CB_4_connectivities,
#    'HA': N_HA_CB_4_connectivities,
#    'CB': N_HA_CB_4_connectivities,
#    #
#    'CA': CA_O_4_connectivities,
#    'O': CA_O_4_connectivities,
#    #
#    }


class Amber14SBForceField:

    __slots__ = [
        'atom_names',
        'bonds_eq3_intra',
        'bonds_le2_intra',
        'forcefield',
        'res_topology',
        ]

    _state = None

    def __new__(cls, *args, **kwargs):

        if cls._state:
            return cls._state

        elif cls._state is None:
            cls._state = super().__new__(cls)
            return cls._state

    def __init__(self, *args, **kwargs):

        self.atom_names = atom_names_amber

        self.forcefield = read_ff14SB_params()

        self.res_topology = generate_residue_template_topology(
            amber_pdbs,
            atom_names_amber,
            **kwargs)

        self.bonds_eq3_intra = topology_3_bonds_apart(self.res_topology)
        self.bonds_le2_intra = expand_topology_bonds_apart(self.res_topology, 2)


forcefields = {
    'Amberff14SB': Amber14SBForceField,
    }



#_amber_res_topology = None
#def generate_amber_residue_topology(**kwargs):
#    global _amber_res_topology
#
#    if not _amber_res_topology:
#        _amber_res_topology = generate_residue_template_topology(
#            amber_pdbs,
#            atom_names_amber,
#            **kwargs)
#
#    return _amber_res_topology




inter_residue_connectivities = {
    2: bonds_le_2_inter,
    3: bonds_le_3_inter,
    #4: inter_4_connect,
    }

# bend angles are in radians
# bend angle for the CA-C-O bond was virtually the same, so will be computed
# as a single value. Changes in CA-C-Np1 create variations in Np1-C-O according
# beacuse CA-C-O is fixed
bend_angles_N_CA_C = {
    'A': 4365251754739323 / 2251799813685248,  # 1.939 0.036
    'R': 8719677895929259 / 4503599627370496,  # 1.936 0.039
    'N': 4379328358183431 / 2251799813685248,  # 1.945 0.046
    'D': 8725333487606339 / 4503599627370496,  # 1.937 0.045
    'C': 272188429672133 / 140737488355328,    # 1.934 0.043
    'E': 4367135841953607 / 2251799813685248,  # 1.939 0.038
    'Q': 2182067275023153 / 1125899906842624,  # 1.938 0.039
    'G': 4450758630407041 / 2251799813685248,  # 1.977 0.045
    'H': 8728063705254555 / 4503599627370496,  # 1.938 0.044
    'I': 1077205649192195 / 562949953421312,   # 1.914 0.041
    'L': 8709813928320865 / 4503599627370496,  # 1.934 0.04
    # 'K': 273144179869987 / 140737488355328, # Example by THG templates
    'K': 8725282730666569 / 4503599627370496,  # 1.937 0.039
    'M': 4357918721408073 / 2251799813685248,  # 1.935 0.039
    'F': 4350621570448131 / 2251799813685248,  # 1.932 0.042
    'P': 4431159883507209 / 2251799813685248,  # 1.968 0.039
    'S': 2185976912710845 / 1125899906842624,  # 1.942 0.041
    'T': 8701671394739109 / 4503599627370496,  # 1.932 0.042
    'W': 4357706976096223 / 2251799813685248,  # 1.935 0.041
    'Y': 2177418951092023 / 1125899906842624,  # 1.934 0.042
    'V': 4310842035620741 / 2251799813685248,  # 1.914 0.04
    }

build_bend_angles_N_CA_C = {
    key: (pi - v) / 2
    for key, v in bend_angles_N_CA_C.items()
    }

bend_angles_CA_C_Np1 = {
    'A': 4588994859787807 / 2251799813685248,  # 2.038 0.022
    'R': 4588049090836895 / 2251799813685248,  # 2.038 0.022
    'N': 4589163477109627 / 2251799813685248,  # 2.038 0.024
    'D': 4591476087195707 / 2251799813685248,  # 2.039 0.023
    'C': 4584298882571901 / 2251799813685248,  # 2.036 0.024
    'E': 4591360788509533 / 2251799813685248,  # 2.039 0.022
    'Q': 2295079292392017 / 1125899906842624,  # 2.038 0.022
    'G': 2292570018001035 / 1125899906842624,  # 2.036 0.026
    'H': 2293539098323251 / 1125899906842624,  # 2.037 0.024
    'I': 1146740789727809 / 562949953421312,   # 2.037 0.021
    'L': 2295421027511699 / 1125899906842624,  # 2.039 0.021
    'K': 2294523703162747 / 1125899906842624,  # 2.038 0.022
    'M': 71693445620547 / 35184372088832,      # 2.038 0.023
    'F': 4583159124932161 / 2251799813685248,  # 2.035 0.023
    'P': 4587486592844193 / 2251799813685248,  # 2.037 0.027
    'S': 2293301456194919 / 1125899906842624,  # 2.037 0.024
    'T': 1146400995781303 / 562949953421312,   # 2.036 0.023
    'W': 1146548728446729 / 562949953421312,   # 2.037 0.024
    'Y': 4582697136860539 / 2251799813685248,  # 2.035 0.023
    'V': 4584841207534447 / 2251799813685248,  # 2.036 0.021
    }

build_bend_angles_CA_C_Np1 = {
    key: (pi - v) / 2
    for key, v in bend_angles_CA_C_Np1.items()
    }

bend_angles_Cm1_N_CA = {
    'A': 4768579151967919 / 2251799813685248,  # 2.118 0.028
    'R': 1192445900065887 / 562949953421312,   # 2.118 0.028
    'N': 1193332907551887 / 562949953421312,   # 2.12 0.03
    'D': 4771817124476497 / 2251799813685248,  # 2.119 0.029
    'C': 4773336800981679 / 2251799813685248,  # 2.12 0.03
    'E': 596106910867665 / 281474976710656,    # 2.118 0.028
    'Q': 2384417995863009 / 1125899906842624,  # 2.118 0.029
    'G': 1194549700095835 / 562949953421312,   # 2.122 0.028
    'H': 4771903805828599 / 2251799813685248,  # 2.119 0.03
    'I': 4770277841981895 / 2251799813685248,  # 2.118 0.029
    'L': 2383227038060725 / 1125899906842624,  # 2.117 0.029
    'K': 1193021713030775 / 562949953421312,   # 2.119 0.029
    'M': 4767963234571985 / 2251799813685248,  # 2.117 0.029
    'F': 4771480556811017 / 2251799813685248,  # 2.119 0.03
    'P': 2386863444043781 / 1125899906842624,  # 2.12 0.03
    'S': 4772350172472667 / 2251799813685248,  # 2.119 0.03
    'T': 4772846148285813 / 2251799813685248,  # 2.12 0.03
    'W': 4773546458813579 / 2251799813685248,  # 2.12 0.032
    'Y': 4773380997634081 / 2251799813685248,  # 2.12 0.031
    'V': 596711317736503 / 281474976710656,    # 2.12 0.029
    }

build_bend_angles_Cm1_N_CA = {
    key: (pi - v) / 2
    for key, v in bend_angles_Cm1_N_CA.items()
    }

# distances are in angstroms
distances_N_CA = {
    'A': 6579089706805643 / 4503599627370496,  # 1.461 0.012
    'R': 822248550821425 / 562949953421312,    # 1.461 0.012
    'N': 3288466758786951 / 2251799813685248,  # 1.46 0.012
    'D': 6581461030414551 / 4503599627370496,  # 1.461 0.012
    'C': 6573024845758137 / 4503599627370496,  # 1.46 0.012
    'E': 6578019054232101 / 4503599627370496,  # 1.461 0.012
    'Q': 6577391610142879 / 4503599627370496,  # 1.46 0.013
    'G': 6551432980914649 / 4503599627370496,  # 1.455 0.013
    'H': 6576300419750417 / 4503599627370496,  # 1.46 0.013
    'I': 1644634584789097 / 1125899906842624,  # 1.461 0.012
    'L': 6576809707492807 / 4503599627370496,  # 1.46 0.012
    'K': 1644714268825877 / 1125899906842624,  # 1.461 0.015
    'M': 6579479788753381 / 4503599627370496,  # 1.461 0.013
    'F': 6574810553972951 / 4503599627370496,  # 1.46 0.012
    'P': 6602907863490325 / 4503599627370496,  # 1.466 0.011
    'S': 6575297426758937 / 4503599627370496,  # 1.46 0.012
    'T': 6573005843763805 / 4503599627370496,  # 1.46 0.012
    'W': 1643809779698609 / 1125899906842624,  # 1.46 0.012
    'Y': 3287342856106863 / 2251799813685248,  # 1.46 0.012
    'V': 6577083632702971 / 4503599627370496,  # 1.46 0.012
    }

average_distance_N_CA = mean(distances_N_CA.values())
std_distance_N_CA = stdev(distances_N_CA.values())

distances_CA_C = {
    'A': 6864424746753997 / 4503599627370496,  # 1.524 0.012
    'R': 6864915463757277 / 4503599627370496,  # 1.524 0.012
    'N': 1716505987820907 / 1125899906842624,  # 1.525 0.012
    'D': 3436265496864801 / 2251799813685248,  # 1.526 0.012
    'C': 1714936829082835 / 1125899906842624,  # 1.523 0.012
    'E': 3433630145140679 / 2251799813685248,  # 1.525 0.012
    'Q': 6864676438206105 / 4503599627370496,  # 1.524 0.012
    'G': 3412585135594759 / 2251799813685248,  # 1.515 0.011
    'H': 6858488656304525 / 4503599627370496,  # 1.523 0.013
    'I': 858607386435583 / 562949953421312,    # 1.525 0.012
    'L': 3431008760147005 / 2251799813685248,  # 1.524 0.012
    'K': 858360898168883 / 562949953421312,    # 1.525 0.015
    'M': 3430952602650511 / 2251799813685248,  # 1.524 0.012
    'F': 1714968063284791 / 1125899906842624,  # 1.523 0.012
    'P': 6865394859536305 / 4503599627370496,  # 1.524 0.012
    'S': 3431691276368069 / 2251799813685248,  # 1.524 0.012
    'T': 6865564932278851 / 4503599627370496,  # 1.524 0.012
    'W': 6859935246968493 / 4503599627370496,  # 1.523 0.013
    'Y': 6859340134216747 / 4503599627370496,  # 1.523 0.012
    'V': 6867874764281747 / 4503599627370496,  # 1.525 0.012
    }

average_distance_CA_C = mean(distances_CA_C.values())
std_distance_CA_C = stdev(distances_CA_C.values())

distances_C_Np1 = {
    'A': 5993805121385571 / 4503599627370496,  # 1.331 0.01
    'R': 2996218862401521 / 2251799813685248,  # 1.331 0.009
    'N': 5993831729486957 / 4503599627370496,  # 1.331 0.01
    'D': 1498633271940567 / 1125899906842624,  # 1.331 0.01
    'C': 2996353095777321 / 2251799813685248,  # 1.331 0.01
    'E': 2996409892006783 / 2251799813685248,  # 1.331 0.01
    'Q': 5992661913157337 / 4503599627370496,  # 1.331 0.009
    'G': 5991694366764827 / 4503599627370496,  # 1.33 0.01
    'H': 5993608389065211 / 4503599627370496,  # 1.331 0.01
    'I': 5994057580226143 / 4503599627370496,  # 1.331 0.009
    'L': 2996685407441657 / 2251799813685248,  # 1.331 0.009
    'K': 5992148258082287 / 4503599627370496,  # 1.331 0.01
    'M': 5993036042256113 / 4503599627370496,  # 1.331 0.01
    'F': 5992325017189537 / 4503599627370496,  # 1.331 0.009
    'P': 2994612945437999 / 2251799813685248,  # 1.33 0.01
    'S': 5992949209364285 / 4503599627370496,  # 1.331 0.01
    'T': 187268519548059 / 140737488355328,    # 1.331 0.01
    'W': 2996589405775763 / 2251799813685248,  # 1.331 0.009
    'Y': 2996393655695127 / 2251799813685248,  # 1.331 0.009
    'V': 2996598670315555 / 2251799813685248,  # 1.331 0.009
    }

average_distance_C_Np1 = mean(distances_C_Np1.values())
std_distance_C_Np1 = stdev(distances_C_Np1.values())

build_bend_CA_C_OXT = (pi - (2 * pi / 3)) / 2
build_bend_CA_C_O = 2.102 / 2
distance_C_OXT = 1.27
distance_C_O = 5556993099130213 / 4503599627370496

# NH atom:
distance_H_N = 0.9
build_bend_H_N_C = np.radians(114) / 2


# side chain template coordinates


def _get_structure_coords(path_):
    s = Structure(path_)
    s.build()
    coords = s.coords.astype(np.float64)
    sidechains = np.where(
        np.logical_not(np.isin(s.data_array[:, col_name], backbone_atoms))
        )[0]

    return coords, sidechains


sidechain_templates = {
    pdb.stem.upper(): _get_structure_coords(pdb)
    for pdb in _sidechain_template_files
    }

# these template coordinates were created using Chimera-X daily given
# a N-terminal at 0,0,0 and a CA along the X axis.
n_terminal_h_coords_at_origin = np.array([
    [ 0.087,  0.76 ,  1.245],
    [ 0.   ,  0.   ,  0.   ],
    [ 1.231,  0.229, -0.869],
    [-0.789,  1.241,  1.392],
    [ 0.837,  1.435,  1.189],
    [ 0.258,  0.13 ,  2.015]
    ])
