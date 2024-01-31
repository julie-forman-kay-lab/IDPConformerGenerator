"""Functions for recognizing and processing multiple protein chains."""
from difflib import SequenceMatcher
from itertools import combinations
from math import ceil, floor

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.SASA import ShrakeRupley

from idpconfgen import Path
from idpconfgen.core.definitions import aa3to1, pk_aa_dict
from idpconfgen.ldrs_helper import consecutive_grouper
from idpconfgen.libs.libparse import split_consecutive_groups
from idpconfgen.libs.libpdb import get_fasta_from_PDB
from idpconfgen.libs.libstructure import (
    Structure,
    col_chainID,
    col_name,
    col_resName,
    col_resSeq,
    col_segid,
    cols_coords,
    )


def process_multichain_pdb(fld_struc, input_seq):
    """
    Search multiple chains in a PDB file.

    Parameters
    ----------
    fld_struc : np.ndarray from IDPConformerGenerator.Structure
        Essentially the .data_array
    
    input_seq : dict
        Dictionary of all the sequences given in the .fasta file

    Returns
    -------
    fld_chainseq : dict
        Dictionary of chains, their matching input sequence, and
        the respective chain structure as an array from the
        template PDB
    """
    fld_resseq = fld_struc[:, col_resSeq]
    fld_resname = fld_struc[:, col_resName]
    fld_chain = fld_struc[:, col_chainID]
    # Check if missing chain ID but there should be segment ID
    if not all(c for c in fld_chain):
        fld_chain = fld_struc[:, col_segid]
    unique_chains = set(fld_chain)
    fld_chainseq = {}
    
    for chain in unique_chains:
        fld_chainseq[chain] = []
    for i, res in enumerate(fld_resseq):
        name = fld_resname[i]
        try:
            next_res = fld_resseq[i + 1]
        except IndexError:
            fld_chainseq[fld_chain[i]].append(aa3to1[name])
            break
        if next_res != res:
            fld_chainseq[fld_chain[i]].append(aa3to1[name])
    for chain in fld_chainseq:
        fld_fasta = ''.join(fld_chainseq[chain]).upper()
        matches = []
        for seq in input_seq:
            in_fasta = input_seq[seq]
            matches.append(SequenceMatcher(None, fld_fasta, in_fasta).ratio())  # noqa: E501
        max_match = max(matches)
        match_index = matches.index(max_match)
        chain_lst = []
        fld_lst = fld_struc.tolist()
        for i, c in enumerate(fld_chain):
            if c == chain:
                chain_lst.append(fld_lst[i])
        chain_arr = np.array(chain_lst)
        fld_chainseq[chain] = (fld_fasta, match_index, chain_arr)
    
    return fld_chainseq


def group_chains(pdb_paths):
    """
    Find fragments of PDBs that belong in the same structure.
    
    Essentially grouping the PDBs in the "splitted" tarball by
    their PDB IDs so ones that belong in the same structure could
    be saved.
    
    Parameters
    ----------
    pdb_paths : list
        List of Path objects. Usually expecting "pdbs2operate".
    
    Returns
    -------
    pdb_groups : list
        List of list of paths for each group of fragments belonging
        to the same protein structure
    """
    assert type(pdb_paths) is list
    file_names = []
    folder_path = pdb_paths[0].parent
    for path in pdb_paths:
        file_names.append(path.stem)
    
    pdb_groups = []
    current_group = [Path(f"{folder_path}/{file_names[0]}.pdb")]
    
    for name in file_names[1:]:
        id = name[:4]
        first_name_id = current_group[0].stem[:4]
        
        if id == first_name_id and name[5] != current_group[0].stem[5]:
            current_group.append(Path(f"{folder_path}/{name}.pdb"))
        else:
            if len(current_group) > 1:
                pdb_groups.append(current_group)
            current_group = [Path(f"{folder_path}/{name}.pdb")]
    
    # Append the last group of files
    if len(current_group) > 1 and current_group[0].stem[5] != current_group[-1].stem[5]:  # noqa: E501
        pdb_groups.append(current_group)
    
    return pdb_groups


def calc_intrachain_ca_contacts(pdb, max_dist):
    """
    Find the residues that are in contact with each other.
    
    Parameters
    ----------
    pdb : Path
        Path to a PDB file of interest.
    
    max_dist : int or float
        Maximum distance allowed between CA to be considered a contact.
        Default is 6.
    
    Returns
    -------
    pdbid : str
        The ID file name
    
    contacts : dict
        Dictionary of pairs of sequences, residues, and CA distances
    
    counter : int
        Number of intermolecular contacts.
    """
    pdb_struc = Structure(pdb)
    pdb_struc.build()
    ca_arr = np.array(
        [pdb_struc.data_array[i]
            for i, data in enumerate(pdb_struc.data_array[:, col_name])
            if data == 'CA']
        )
    ca_coordinates = ca_arr[:, cols_coords].astype(float)

    with open(pdb) as f:
        pdb_raw = f.read()
    pdbid, fasta = get_fasta_from_PDB([pdb, pdb_raw])

    num_residues = len(ca_coordinates)
    contacts = []
    counter = 0
    
    for i in range(num_residues):
        for j in range(i + 1, num_residues):
            d_ca = np.linalg.norm(np.array(ca_coordinates[i]) - np.array(ca_coordinates[j]))  # noqa: E501
            ca_dists = []
            chain1_seq = ""
            chain2_seq = ""
            # Euclidian distance must be within range (default is 6 A)
            # residues must be at least 5 apart
            if d_ca <= max_dist and j > i + 4:
                c1_count = 0
                for k in range(i - 2, i + 3):
                    if 0 <= k < num_residues:
                        d_ca = np.linalg.norm(np.array(ca_coordinates[k]) - np.array(ca_coordinates[j]))  # noqa: E501
                        ca_dists.append(d_ca)
                        chain1_seq += f"{k},"
                        c1_count += 1
                half_res = c1_count / 2
                upper = ceil(half_res)
                lower = floor(half_res)
                for k in range(j - lower, j + upper):
                    if 0 <= k < num_residues:
                        chain2_seq += f"{k},"

                contacts.append({chain1_seq: [chain2_seq, ca_dists]})
                counter += 1
    
    if counter == 0:
        return False

    return pdbid, contacts, counter


def calc_interchain_ca_contacts(pdb_groups, max_dist):
    """
    Find CA contacts below a certain distance between different chains.
    
    Returns contact information based on specific chains.
    Has an additional datapoint compared to intrachain CA contacts.
    
    Format of chain_contacts:
    [
        {
            "SEQ_IDX":[
                "pdbid_B_seg0",
                "SEQ_IDX",
                [
                    distance1,
                    distance2,
                    distance3,
                    ...
                ]
            ]
        },
        ...
    ]
    
    Parameters
    ----------
    pdb_groups : list of lists of Path
        Groups of Path variables for multiple segments belonging to the
        same PDB structure.
    
    max_dist : float or int
        Maximum distance between CA to be considered a contact.
    
    Returns
    -------
    pdb_ids : list of str
        List of the different chains which the contacts at
        a specific list position belongs to.
    
    chain_contacts : list of dict
        List of different contacts between chains.
        Contacts are stored as dict as indicated above.
    
    counter : int
        Number of intermolecular contacts.
    """
    chain_contacts = {}
    counter = 0
    pdb_ids = []
    pdb_strucs = {}
    
    for path in pdb_groups:
        struc = Structure(path)
        struc.build()
        struc_arr = struc.data_array
        chain = path.stem[5]
        seg = path.stem[10]
        pdb_strucs[chain + seg] = (struc_arr, path)
    
    chains = list(pdb_strucs.keys())
    chain_combos = [combo for combo in combinations(chains, 2) if combo[0][0] != combo[1][0]]  # noqa: E501
    
    for combo in chain_combos:
        c1 = combo[0]
        c2 = combo[1]
        c1_arr = pdb_strucs[c1][0]
        c2_arr = pdb_strucs[c2][0]
        c1_path = pdb_strucs[c1][1]
        c2_path = pdb_strucs[c2][1]
        
        with open(c1_path) as f1:
            c1_raw = f1.read()
        c1_pdbid, _ = get_fasta_from_PDB([c1_path, c1_raw])
        with open(c2_path) as f2:
            c2_raw = f2.read()
        c2_pdbid, _ = get_fasta_from_PDB([c2_path, c2_raw])
        
        c1_ca_arr = np.array([c1_arr[i] for i, data in enumerate(c1_arr[:, col_name]) if data == 'CA'])  # noqa: E501
        c1_ca_coords = c1_ca_arr[:, cols_coords].astype(float)
        c2_ca_arr = np.array([c2_arr[i] for i, data in enumerate(c2_arr[:, col_name]) if data == 'CA'])  # noqa: E501
        c2_ca_coords = c2_ca_arr[:, cols_coords].astype(float)
        
        c1_res_tot = len(c1_ca_coords)
        c2_res_tot = len(c2_ca_coords)

        for i in range(c1_res_tot):
            for j in range(c2_res_tot):
                d_ca = np.linalg.norm(np.array(c1_ca_coords[i]) - np.array(c2_ca_coords[j]))  # noqa: E501
                ca_dists = []
                chain1_seq = ""
                chain2_seq = ""
                # Euclidian distance must be within range (default is 12 A)
                if d_ca <= max_dist:
                    c1_count = 0
                    for k in range(i - 2, i + 3):
                        if 0 <= k < c1_res_tot:
                            d_ca = np.linalg.norm(np.array(c1_ca_coords[k]) - np.array(c2_ca_coords[j]))  # noqa: E501
                            ca_dists.append(d_ca)
                            chain1_seq += f"{k},"
                            c1_count += 1
                    half_res = c1_count / 2
                    upper = ceil(half_res)
                    lower = floor(half_res)
                    for k in range(j - lower, j + upper):
                        if 0 <= k < c2_res_tot:
                            chain2_seq += f"{k},"
                    
                    if c1_pdbid not in chain_contacts:
                        chain_contacts[c1_pdbid] = []
                    chain_contacts[c1_pdbid].append({chain1_seq: [c2_pdbid, chain2_seq, ca_dists]})  # noqa: E501
                    pdb_ids.append(c1_pdbid)
                    counter += 1
    
    if counter == 0:
        return False
    
    return pdb_ids, chain_contacts, counter


def has_consecutive_match(s, q, consecutive_length=2):
    """Find matches of consecutive characters in s within q."""
    for i in range(len(q) - 1):
        consecutive_chars_q = q[i:i + consecutive_length]
        if consecutive_chars_q in s:
            return True
        
    return False


def contact_matrix(db, sequence):
    """
    Generate a matrix of the sequence against database of pairs of contacts.
    
    Parameters
    ----------
    db : dict of list of tuple pairs
        Minimized version of the database that only contains
        information of the sequence pairs. There will also be
        positional information so matches can be back-mapped
        to the full database.
        E.g. {"seg1": [("ABC", "DEF"),]}
    
    sequence : list or str
        Combination of input protein sequences.
        Can be single or dual sequences for intra- or
        intermolecular contacts respectfully.
    
    Returns
    -------
    matrix : np.ndarray
        Disitrubtion matrix of all the matches and weights

    positions : list
        List of positions corresponding to input database
    """
    segid = next(iter(db))
    
    if type(sequence) is list:
        seq1 = sequence[0]
        len1 = len(seq1)
        seq2 = sequence[1]
        len2 = len(seq2)
        
        h_mtx = np.zeros((len1, len2))
        l_mtx = np.empty((len1, len2), dtype=object)
        l_mtx.fill([])
        
        for pairs in db.values():
            for idx, pair in enumerate(pairs):
                p1, p2 = pair
                for i in range(len1):
                    for j in range(len2):
                        c1 = ""
                        c2 = ""
                        for k in range(i - 2, i + 3):
                            if 0 <= k < len1:
                                c1 += f"{seq1[k]}"
                        for k in range(j - 2, j + 3):
                            if 0 <= k < len2:
                                c2 += f"{seq2[k]}"
                        p1_c1 = has_consecutive_match(p1, c1)
                        p1_c2 = has_consecutive_match(p1, c2)
                        p2_c1 = has_consecutive_match(p2, c1)
                        p2_c2 = has_consecutive_match(p2, c2)
                        if ((p1_c1 and p2_c2) or (p1_c2 and p2_c1)):
                            h_mtx[i, j] += 1
                            l_mtx[i, j].append(idx)
        
        hit_matrix = np.flipud(h_mtx)
        loc_matrix = np.flipud(l_mtx)
    else:
        seq_len = len(sequence)
        hit_matrix = np.zeros((seq_len, seq_len))
        loc_matrix = np.empty((seq_len, seq_len), dtype=object)
        loc_matrix.fill([])
        
        for pairs in db.values():
            for idx, pair in enumerate(pairs):
                p1, p2 = pair
                for i in range(seq_len):
                    for j in range(i + 1, seq_len):
                        c1 = ""
                        c2 = ""
                        for k in range(i - 2, i + 3):
                            if 0 <= k < seq_len:
                                c1 += f"{sequence[k]}"
                        for k in range(j - 2, j + 3):
                            if 0 <= k < seq_len:
                                c2 += f"{sequence[k]}"
                        p1_c1 = has_consecutive_match(p1, c1)
                        p1_c2 = has_consecutive_match(p1, c2)
                        p2_c1 = has_consecutive_match(p2, c1)
                        p2_c2 = has_consecutive_match(p2, c2)
                        if ((p1_c1 and p2_c2) or (p1_c2 and p2_c1)) and j > i + 4:  # noqa: E501
                            hit_matrix[i, j] += 1
                            hit_matrix[j, i] += 1
                            loc_matrix[i, j].append(idx)

        hit_matrix = np.flipud(hit_matrix)
        loc_matrix = np.flipud(loc_matrix)
    return segid, hit_matrix, loc_matrix


def find_sequence_net_charge(seq, pH):
    """
    Find the net charge of a given sequence based on pK and pH.

    Parameters
    ----------
    seq : str
        String of 1 letter AAs.
    
    pH : int or float
        pH of interest.
    
    Returns
    -------
    net_charge : list of float
        Net charge per residue.
    """
    aa = [*seq]
    net_charge = []
    for a in aa:
        charge = 0
        # Accounts for phosphorylation, add a negative charge
        # Uses the pKa of phosphoate bioesters
        if a.islower():
            if pH > 1.5:
                charge -= 1
            elif pH == 1.5:
                charge -= 0.5
            if pH > 6.3:
                charge -= 1
            elif pH == 6.3:
                charge -= 0.5
            a = a.upper()
            
        pKa, pKb, pKx = pk_aa_dict[a]
        
        if pKa > pH:
            charge += 1
        elif pKa == pH:
            charge += 0.5
        if pKb < pH:
            charge -= 1
        elif pKb == pH:
            charge -= 0.5
        
        if a in ["R", "H", "K"]:
            if pKx > pH:
                charge += 1
            elif pKx == pH:
                charge += 0.5
        elif a in ["D", "C", "E", "Y"]:
            if pKx < pH:
                charge -= 1
            elif pKx == pH:
                charge -= 0.5
        net_charge.append(charge)
    
    return net_charge


def electropotential_matrix(sequences, pH=7.0):
    """
    Generate a matrix of contacts based on pH and residue electrostatics.
    
    Parameters
    ----------
    sequences : list or str
        Can accept single or dual sequences based on
        intra- or intermolecular contacts respectively.
    
    pH : int or float
        Desired pH value of interest.
        Defaults to 7.0.
    
    Returns
    -------
    matrix : np.ndarray
        Disitrubtion matrix of all the matches and weights
    """
    assert (type(sequences) is list) or (type(sequences) is str)

    if type(sequences) is list:
        seq1 = sequences[0]
        len1 = len(seq1)
        seq2 = sequences[1]
        len2 = len(seq2)
        charge1 = find_sequence_net_charge(seq1, pH)
        charge2 = find_sequence_net_charge(seq2, pH)
        
        charge_matrix = np.zeros((len1, len2))
        for i in range(len1):
            for j in range(len2):
                window1 = charge1[i:i + 5]
                window2 = charge2[j:j + 5]
                avg1 = sum(window1) / 5
                avg2 = sum(window2) / 5
                charge_matrix[i, j] = abs(avg1 - avg2)
        charge_matrix = np.flipud(charge_matrix)
    else:  # For single sequences only, intramolecular
        len_seq = len(sequences)
        net_charge = find_sequence_net_charge(sequences, pH)
        
        charge_matrix = np.zeros((len_seq, len_seq))
        for i in range(len_seq):
            for j in range(i + 5, len_seq):
                window1 = net_charge[i:i + 5]
                window2 = net_charge[j:j + 5]
                avg1 = sum(window1) / 5
                avg2 = sum(window2) / 5
                charge_matrix[i, j] = abs(avg1 - avg2)
                charge_matrix[j, i] = abs(avg1 - avg2)
        charge_matrix = np.flipud(charge_matrix)
    
    return charge_matrix


def extract_intrapairs_from_db(intra_seg):
    """
    Extract the sequence pairs of intra- contacts.

    Parameters
    ----------
    intra_seg : dict
        PDB segment from IDPConformerGenerator extended database
    
    Return
    ------
    contact_pairs : dict
        dict of list of tuple pairs of sequences in order
        Key value being the seg ID from database entry
    """
    seg_id = next(iter(intra_seg))
    seg_vals = list(intra_seg.values())[0]
    seq = seg_vals["fasta"]
    contact_pairs = {}
    if "intra" in seg_vals:
        contact_pairs[seg_id] = []
        intra = seg_vals["intra"]
        for contact in intra:
            s1 = next(iter(contact))
            s2 = contact[s1][0]
            f1 = s1.split(",")
            f1.pop(-1)
            f2 = s2.split(",")
            f2.pop(-1)
            r1 = ""
            r2 = ""
            for r in f1:
                r1 += seq[int(r)]
            for r in f2:
                r2 += seq[int(r)]
            contact_pairs[seg_id].append((r1, r2))
    else:
        return False

    return contact_pairs


def extract_interpairs_from_db(inter_segs):
    """
    Extract the sequence pairs of inter- contacts.

    Parameters
    ----------
    inter_seg : dict
        PDB segment from IDPConformerGenerator extended database
    
    Return
    ------
    contact_pairs : dict
        dict of list of pairs of sequences in order for specific segment
    """
    contact_pairs = {}
    for seg in inter_segs:
        seg_vals = inter_segs[seg]
        
        if "inter" in seg_vals:
            contact_pairs[seg] = []
            inter = seg_vals["inter"]
            s1_seq = seg_vals["fasta"]
            for contact in inter:
                s1 = next(iter(contact))
                s2_seg = contact[s1][0]
                try:
                    s2_seq = inter_segs[s2_seg]["fasta"]
                except KeyError:
                    return False
                s2_seq_idx = contact[s1][1]
                f1 = s1.split(",")
                f1.pop(-1)
                f2 = s2_seq_idx.split(",")
                f2.pop(-1)
                r1 = ""
                r2 = ""
                for r in f1:
                    r1 += s1_seq[int(r)]
                for r in f2:
                    r2 += s2_seq[int(r)]
                contact_pairs[seg].append((r1, r2))
    
    return contact_pairs


def pick_point_from_heatmap(prob_mtx, max_num=1):
    """
    Select one point from the probability heatmap.
    
    Parameters
    ----------
    prob_mtx : 2D np array heatmap
        Should be normalized ranging from 0 - 1.
    
    max_num : int
        Maximum number of contacts to select from matrix.
    
    Returns
    -------
    picked_points : list of tuple
        Coordinates of where the selected point is
    """
    flat_heatmap = prob_mtx.flatten()
    
    prob_distribution = flat_heatmap / np.sum(flat_heatmap)
    
    num_to_pick = np.random.randint(1, max_num + 1)
    
    picked_points = []
    for _ in range(num_to_pick):
        picked_index = np.random.choice(
            len(prob_distribution),
            p=prob_distribution
            )
        _, num_cols = prob_mtx.shape
        x = picked_index // num_cols
        y = picked_index % num_cols
        picked_points.append((x, y))
    
    return picked_points


def calculate_max_contacts(sequences):
    """
    Calculate the maximum number of possible contacts.
    
    Calculation is based on multiples of '5' since maximum of
    up to 5 residues is stored in the extended database.

    Parameters
    ----------
    sequences : dict or str
        Can accept multiple sequences based on `input_seq`
        variable from `cli_complex.py`.
    
    Returns
    -------
    num_contacts : dict
        Key is the chain value is the maximum number of contacts
        that chain can make.
    """
    assert (type(sequences) is dict) or (type(sequences) is str)
    num_contacts = {}
    if type(sequences) is dict:
        chain = list(sequences.keys())
        seq = list(sequences.values())
        len_seq = [len(s) for s in seq]
        for c, idx in enumerate(chain):
            if c == "":
                c = "A"
            num_contacts[c] = len_seq[idx] // 5
    else:
        len_seq = len(sequences)
        num_contacts['A'] = len_seq // 5
    
    return num_contacts


def find_sa_residues(
        structure_path,
        min_area=31.65,
        probe_radius=1.40,
        n_points=500
        ):
    """
    Use SASA module in BioPython to find surface accessible resiudes.
    
    Reference: https://biopython.org/docs/dev/api/Bio.PDB.SASA.html

    Parameters
    ----------
    structure_path : string
        Path to the .PDB structure of interest.
    
    min_area : float
        Minimum SASA to consider to be on the surface.
        Defaults to 31.65 Å^2, half the surface area of Glycine.
    
    probe_radius : float
        Radius of rolling-ball, defaults to 1.4 Å,
        radius of a water molecule.
    
    n_points : int
        Resolution of the surface of each atom.
        Defaults to 500.
    
    Returns
    -------
    residue_idx : dict
        Indices for the consecutive residues that are
        considered on the surface to make a contact.
        Keys are the chain ID of interest.
    """
    if structure_path.endswith('.pdb'):
        p = PDBParser(QUIET=1)
    elif structure_path.endswith('.cif'):
        p = MMCIFParser(QUIET=1)
    
    sr = ShrakeRupley(
        probe_radius=probe_radius,
        n_points=n_points,
        )
    struc = p.get_structure("structure", structure_path)
    sr.compute(struc, level="R")
    
    residue_idx = {}
    for chain in struc.get_chains():
        temp_res_idx = []
        for i, res in enumerate(list(chain.get_residues())):
            area = res.sasa
            if area > min_area:
                temp_res_idx.append(i)
        consecutive_idx = split_consecutive_groups(temp_res_idx)
        residue_idx[chain.id] = consecutive_idx

    # Check that the last set of indices does not go over
    fld_struc = Structure(Path(structure_path))
    fld_struc.build()
    for chain, sequence in fld_struc.fasta.items():
        last_idx = residue_idx[chain][-1][-1]
        if last_idx == len(sequence):
            residue_idx[chain][-1][-1] = residue_idx[chain][-1][-1][:-1]
    
    return residue_idx


def create_coordinate_combinations(data, modifier=0):
    """
    Create coordinates based on a tuple of two sets of residues of interest.
    
    Parameters
    ----------
    data : tuple, N = 2
        Tuple of lists ([[]], [[]])
    
    modifier : int
        Change the values of coordinates for indexing purposes
    
    Returns
    -------
    coordinates : list of tuple
        Combination of coordinates between each of the lists of lists
        within the tuple [(,)]
    """
    if len(data) != 2:
        raise ValueError("Input must be a tuple of two elements.")

    list1, list2 = data

    coordinates = []
    for sublist1 in list1:
        for item1 in sublist1:
            for sublist2 in list2:
                for item2 in sublist2:
                    coordinates.append((item1 + modifier, item2 + modifier))

    return coordinates


def process_custom_contacts(file, res_combos, chain_combos):
    """
    Process text (.txt) file containing inter- and/or intra- contacts.
    
    Follows the following format:
        Each new line is a new contact.
        Inter- and intra-contacts are denoted by a single slash (/).
        IDP sequences are given as string (multi-character).
        Chains are given as single upper case character.
        
    For example:
        seq1:1,2,3,4,5/A:33,34,35,36,37
        seq1:6,7,8,9/seq2:1,2,3,4,5,6,7
        seq2:20,21,27,35,39/B:1,10,12,13,14
        seq1:15,17,19,20,21/seq1:37,38,40,42,44

    Parameters
    ----------
    file : string or Path
        Path to the custom contacts file of interest.
    
    res_combos : list
        List of combinations of residues for each of the respected chains

    chain_combos : list
        List of combinations of different chains and sequences by their
        name.
    
    Returns
    -------
    custom_res_combo : list
        New residue combinations now with certain resiudes
    
    positions : list
        List of positions/coordinates for each of the combinations
        we know for sure to make a contact.
    """
    custom_res_combo = []
    custom_chain_combo = []
    positions = []
    
    with open(file) as cc_file:
        lines = cc_file.readlines()
        for line in lines:
            try:
                splitted = line.split("/")
                chain1 = splitted[0]
                chain2 = splitted[1]
                
                chain1_split = chain1.split(":")
                chain2_split = chain2.split(":")
                
                chain1ID = chain1_split[0]
                chain2ID = chain2_split[0]
                chain1Seq_grouped = \
                    consecutive_grouper([int(s) for s in chain1_split[1].split(',')])  # noqa: E501
                chain2Seq_grouped = \
                    consecutive_grouper([int(s) for s in chain2_split[1].split(',')])  # noqa: E501
                
                chain1Seq = []
                chain2Seq = []
                for g1 in chain1Seq_grouped:
                    chain1Seq.append([s for s in range(g1[0], g1[1])])
                for g2 in chain2Seq_grouped:
                    chain2Seq.append([s for s in range(g2[0], g2[1])])
                
                if len(chain1ID) > 1:
                    chain1ID = f">{chain1ID}"
                if len(chain2ID) > 1:
                    chain2ID = f">{chain2ID}"
                
                custom_chain_combo.append((chain1ID, chain2ID))
                custom_res_combo.append((chain1Seq, chain2Seq))
            except Exception:
                # We have to skip lines that don't follow the correct formatting
                continue
    
    for c, combo in enumerate(custom_chain_combo):
        to_add = -1
        c1, c2 = combo
        if len(c1) == 1 and len(c2) == 1:
            # Ignore improperly formatted files
            # i.e. cannot have both folded domains
            continue
        elif len(c1) > 1 and len(c2) > 1:
            # Get the positions straight up since they're both IDPs
            temp_coordinates = create_coordinate_combinations(
                data=custom_res_combo[c],
                modifier=-1
                )
            positions.append(temp_coordinates)
            # No need to update the residue combinations
            continue
        try:
            combo_rev = combo[::-1]
            if combo in chain_combos:
                idx = chain_combos.index(combo)
            elif combo_rev in chain_combos:
                custom_chain_combo[c] = combo_rev
                custom_res_combo[c] = custom_res_combo[c][::-1]
                idx = chain_combos.index(combo_rev)
            else:
                # If the chains do not match input files
                continue
        except ValueError:
            # Catches improper formatting
            continue

        # Make sure we're only considering the folded domain
        if c1.isupper() and len(c1) == 1:
            to_add = 0
        elif c2.isupper() and len(c2) == 1:
            to_add = 1
        
        custom_fld_res = custom_res_combo[c][to_add]  # noqa: F841
        original_fld_res = res_combos[idx]  # noqa: F841
        
    return custom_res_combo, positions
        

def reverse_position_lookup(coords, location_mtx, database):
    """
    Return database entry based on a point in the contacts frequency heatmap.

    Parameters
    ----------
    coords : tuple of int
        The location of the contact per the contacts frequency heatmap.
    
    location_mtx : dict of np.ndarray
        Keys are the segid as seen in the database. Values are an array
        of list of indicies of where the contact is.
    
    database : dict
        IDPConformerGenerator database
    
    Returns
    -------
    entry : dict
        Key-value pairs of the entry for a specific segid in the database.
        TODO: this return can be changed to only the torsion angle? or remain
        as we need information on the secondary structure as well
    """
    return
