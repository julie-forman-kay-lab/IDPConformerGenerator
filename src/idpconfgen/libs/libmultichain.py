"""Functions for recognizing and processing multiple protein chains."""
from difflib import SequenceMatcher
from itertools import combinations
from math import ceil, floor

import numpy as np

from idpconfgen import Path
from idpconfgen.core.definitions import aa3to1
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
    
    sequence : str
        Input protein sequence
    
    ncores : int
        Number of workers to execute search
    
    Returns
    -------
    matrix : np.ndarray
        Disitrubtion matrix of all the matches and weights

    positions : list
        List of positions corresponding to input database
    """
    seq_len = len(sequence)
    hit_matrix = np.zeros((seq_len, seq_len))

    for pairs in db.values():
        for pair in pairs:
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
    
    return hit_matrix


def extract_intrapairs_from_db(intra_seg):
    """
    Extract the sequence pairs of intra- contacts.

    Parameters
    ----------
    intra_seg : dict
        PDB segment from IDPConformerGenerator extended database
    
    Return
    ------
    contact_pairs : list
        List of tuple pairs of sequences in order
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
