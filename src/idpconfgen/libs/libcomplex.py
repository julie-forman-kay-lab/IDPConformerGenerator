"""Functions that are useful for generating dynamic complexes."""
from collections import Counter
from random import randint

import numpy as np

from idpconfgen import log
from idpconfgen.core.definitions import pk_aa_dict
from idpconfgen.libs.libparse import has_consecutive_match


contact_type = ["intra", "inter"]


def contact_matrix(database, sequence):
    """
    Generate a matrix of the sequence against database of pairs of contacts.
    
    Parameters
    ----------
    database : tuple of the following
        contact_db : dict of list of tuple pairs
            Minimized version of the database that only contains
            information of the sequence pairs. There will also be
            positional information so matches can be back-mapped
            to the full database.
            E.g. {"seg1": [("ABC", "DEF"),]}

        dist_db : dict of list of tuples
            Aligned with contact_db but for the distances associated
            with the pairs.
            E.g. {"seg1": [(1, 2, 3),]}
    
    sequence : list or str
        Combination of input protein sequences.
        Can be single or dual sequences for intra- or
        intermolecular contacts respectfully.
    
    Returns
    -------
    hit_matrix : np.ndarray
        Disitrubtion matrix of all the matches and weights

    dist_matrix : np.ndarray
        Matrix of lists corresponding to distances from database
    """
    contact_db = database[0]
    dist_db = database[1]
    
    segid = next(iter(contact_db))
    # Make sure databases are aligned
    assert segid == next(iter(dist_db))
    
    if type(sequence) is list:
        seq1 = sequence[0]
        len1 = len(seq1)
        seq2 = sequence[1]
        len2 = len(seq2)
        
        h_mtx = np.zeros((len1, len2))
        d_mtx = np.empty((len1, len2), dtype=object)
        d_mtx.fill([])
        
        for pairs in contact_db.values():
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
                            min_len = min(len(c1), len(c2))
                            sub_distances = list(dist_db.values())[0][idx][:min_len]  # noqa: E501
                            d_mtx[i, j].append(sub_distances)
        
        hit_matrix = np.flipud(h_mtx)
        dist_matrix = np.flipud(d_mtx)
    else:
        seq_len = len(sequence)
        h_mtx = np.zeros((seq_len, seq_len))
        d_mtx = np.empty((seq_len, seq_len), dtype=object)
        d_mtx.fill([])
        
        for pairs in contact_db.values():
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
                            h_mtx[i, j] += 1
                            h_mtx[j, i] += 1
                            min_len = min(len(c1), len(c2))
                            sub_distances = list(dist_db.values())[0][idx][:min_len]  # noqa: E501
                            d_mtx[i, j].append(sub_distances)
                            d_mtx[j, i].append(sub_distances)

        hit_matrix = np.flipud(h_mtx)
        dist_matrix = np.flipud(d_mtx)
    return hit_matrix, dist_matrix


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
    
    contact_dists : dict
        Dictionary of list of tuples of distances for each residue
        in the specfic pair. Key value being seg ID from database entry
    """
    seg_id = next(iter(intra_seg))
    seg_vals = list(intra_seg.values())[0]
    seq = seg_vals["fasta"]
    contact_pairs = {}
    contact_dists = {}
    if contact_type[0] in seg_vals:
        contact_pairs[seg_id] = []
        contact_dists[seg_id] = []
        intra = seg_vals[contact_type[0]]
        for contact in intra:
            s1 = next(iter(contact))
            s2 = contact[s1][0]
            dist = tuple(contact[s1][1])
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
            contact_dists[seg_id].append(dist)
    else:
        return False, False

    return contact_pairs, contact_dists


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
    
    contact_dists : dict
        Dictionary of list of distances between pairs of sequences
    """
    contact_pairs = {}
    contact_dists = {}
    for seg in inter_segs:
        seg_vals = inter_segs[seg]
        
        if contact_type[1] in seg_vals:
            contact_pairs[seg] = []
            contact_dists[seg] = []
            inter = seg_vals[contact_type[1]]
            s1_seq = seg_vals["fasta"]
            for contact in inter:
                s1 = next(iter(contact))
                s2_seg = contact[s1][0]
                dist = tuple(contact[s1][2])
                try:
                    s2_seq = inter_segs[s2_seg]["fasta"]
                except KeyError:
                    return False, False
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
                contact_dists[seg].append(dist)
    
    return contact_pairs, contact_dists


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


def process_custom_contacts(
        file,
        combo_chains,
        combo_res,
        input_seq,
        ignore_sasa=False
        ):
    """
    Process text (.txt) file containing inter- and/or intra- contacts.
    
    Lets the user know if a set of residues are not considered to be SASA
    by comparing to `combo_res`. For unspecified contacts in `combo_chains`,
    store as "none" in the final list of `custom_res`.
    
    Follows the following format:
        Each new line is a new contact.
        Inter- and intra-contacts are denoted by a single slash (/).
        IDP sequences are given as string (multi-character).
        Chains are given as single upper case character.
        
        Folded chains MUST be before IDP sequences.
        
    For example:
        A:33,34,35,36,37/seq1:1,2,3,4,5
        seq1:6,7,8,9/seq2:1,2,3,4,5,6,7
        B:1,10,12,13,14/seq2:20,21,27,35,39/
        seq1:15,17,19,20,21/seq1:37,38,40,42,44

    Parameters
    ----------
    file : string or Path
        Path to the custom contacts file of interest.
    
    combo_chains : list
        List of combinations of sequences or sequences and chains.

    combo_res : list
        List of combinations of different chains and sequences by their
        name.
    
    input_seq : dict
        Dictionary of input sequences for intramolecular contacts.
    
    ignore_sasa : Bool
        Whether or not to include custom sequences that are not found to be
        surface accessible.
    
    Returns
    -------
    cus_inter_res : list
        Combinations of interchain residues from custom contacts.
        Aligned to combo_chains.
    
    cus_intra_res : dict
        Custom contacts for intraIDP residues
    """
    custom_res = []
    custom_chains = []
    
    with open(file) as cc_file:
        lines = cc_file.readlines()
        for line in lines:
            try:
                splitted = line.split("/")
                chain1 = splitted[0]  # Must be chain of folded domain
                chain2 = splitted[1]
                
                chain1_split = chain1.split(":")
                chain2_split = chain2.split(":")
                
                chain1ID = chain1_split[0]
                chain2ID = chain2_split[0]
                chain1Seq = [[int(s) for s in chain1_split[1].split(',')]]  # noqa: E501
                chain2Seq = [[int(s) for s in chain2_split[1].split(',')]]  # noqa: E501

                if len(chain1ID) > 1:
                    chain1ID = f">{chain1ID}"
                if len(chain2ID) > 1:
                    chain2ID = f">{chain2ID}"
                
                custom_chains.append((chain1ID, chain2ID))
                custom_res.append((chain1Seq, chain2Seq))
            except Exception:
                # We have to skip lines that don't follow the correct formatting
                log.info(f"Incorrect formatting detected for line: {line}. Skipping.")  # noqa: E501
                continue

    og_combos = {}
    for i, chain in enumerate(combo_chains):
        og_combos[chain] = combo_res[i]
    
    custom_combos = {}
    ordered_combos = {}
    cus_intra_res = {}
        
    for i, ccombo in enumerate(custom_chains):
        if ccombo not in combo_chains:
            # Should be an intra-contact
            if ccombo in list(cus_intra_res.keys()):
                group = cus_intra_res[ccombo]
                group[0].append(custom_res[i][0][0])
                group[1].append(custom_res[i][1][0])
                cus_intra_res[ccombo] = group
            else:
                cus_intra_res[ccombo] = custom_res[i]
        else:
            # Should be an inter-contact
            if ccombo in list(custom_combos.keys()):
                group = custom_combos[ccombo]
                group[0].append(custom_res[i][0][0])
                group[1].append(custom_res[i][1][0])
                custom_combos[ccombo] = group
            else:
                custom_combos[ccombo] = custom_res[i]
    
    for key in list(input_seq.keys()):
        intra_combo = (key, key)
        if intra_combo not in list(cus_intra_res.keys()):
            cus_intra_res[intra_combo] = None
        
    for all_ccombo in combo_chains:
        ordered_combos[all_ccombo] = None
        if all_ccombo not in custom_chains:
            custom_combos[all_ccombo] = None
    
    sorted_chains = sorted(custom_combos, key=lambda x: list(ordered_combos.keys()).index(x))  # noqa: E501
    custom_combos = {key: custom_combos[key] for key in sorted_chains}
    
    for combo, res in og_combos.items():
        chain = combo[0]
        if len(chain) == 1:
            if custom_combos[combo] is None:
                continue
            else:
                custom_res = custom_combos[combo][0][0]
                chain_res = res[0]
                flat_res = [r for subr in chain_res for r in subr]
                non_sasa = set(custom_res) - set(flat_res)
                if not ignore_sasa:
                    log.info(
                        f"Please note that residue {non_sasa} in chain "
                        f"{chain} is not identified to be surface accessible. "
                        "Since --ignore-sasa has not been initialized, these "
                        "residues will be removed from consideration."
                        )
                    custom_combos[combo][0][0] = [val for val in custom_res if val not in non_sasa]  # noqa: E501
                else:
                    log.info(
                        f"Please note that residue {non_sasa} in chain "
                        f"{chain} is not identified to be surface accessible. "
                        "Since --ignore-sasa has been initialized, these "
                        "residues will not be removed from consideration."
                        )
        else:
            continue
    
    cus_inter_res = list(custom_combos.values())
    
    return cus_inter_res, cus_intra_res


def select_contacts(
        coords,
        min_x_spacing=3,
        max_num_points=10,
        ):
    """
    Select random number of contacts based on the contact heatmap.
    
    Parameters
    ----------
    coords : np.array
        Contact heatmap of choice.
    
    min_x_spacing : int
        Minimum number of residues to space out for in the X-axis.
        Defaults to 3.
    
    max_num_points : int
        Maximum number of points to select from contact heatmap.
        Defaults to 10.
    
    Returns
    -------
    x_coordinates : tuple
        Indices for the coordinates on the X-axis.
    
    y_coordinates : tuple
        Indices for the coordinates on the Y-axis.
    """
    y_size, x_size = coords.shape
    
    num_points = min(max_num_points, x_size // min_x_spacing)

    # Calculate probabilities for each point in the heatmap
    flat_heatmap = coords.flatten()
    probabilities = flat_heatmap / np.sum(flat_heatmap)
    
    # Choose random points based on probabilities
    chosen_indices = np.random.choice(
        np.arange(x_size * y_size),
        size=num_points,
        replace=False,
        p=probabilities
        )

    x_coordinates, y_coordinates = np.unravel_index(chosen_indices, (y_size, x_size))  # noqa: E501

    chosen_points = []
    while len(chosen_points) < num_points:
        index = np.random.choice(np.arange(x_size * y_size), p=probabilities)
        x, y = np.unravel_index(index, (y_size, x_size))
        
        if not chosen_points or abs(chosen_points[-1][0] - x) >= min_x_spacing:
            chosen_points.append((x, y))
            if len(chosen_points) == num_points:
                break

    # Sort the chosen points by x coordinate to ensure they are consecutive
    chosen_points.sort(key=lambda point: point[0])
    y_coordinates, x_coordinates = zip(*chosen_points)

    return x_coordinates, y_coordinates


def select_custom_contacts(contacts, idx, c_type, max_contacts):
    """
    Similar to select_contacts() but built for custom contacts.

    Parameters
    ----------
    contacts : list or dict
        Inter- and Intra-contacts respectively depending on the datatype.
        For the list case, it should be aligned by indices to the combination
        of chains determined from the main CLI.
    
    idx : int
        Index aligned to the combination of chains determined from the main CLI.
    
    c_type : string
        Type of inter- or intra-contact. Check against datatype of contacts.
    
    max_contacts : int
        Maximum number of contacts to have.
        Need to check we're not returning greater than the maximum.
    
    Returns
    -------
    list pair of residue contacts.
    
    False if no matches are found or invalid formatting.
    """
    if c_type is contact_type[0]:
        assert type(contacts) is dict
        custom_contacts = list(contacts.values())[idx]
        if custom_contacts is None:
            return False
        selected = randint(0, len(custom_contacts[0]) - 1)
        if max_contacts >= len(custom_contacts[0]):
            # Return all the contacts because it's within the max
            return custom_contacts[0], custom_contacts[1]
        else:
            # Return one contact
            return custom_contacts[0][selected], custom_contacts[1][selected]
    
    elif c_type is contact_type[1]:
        assert type(contacts) is list
        custom_contacts = contacts[idx]
        if custom_contacts is None:
            return False
        selected = randint(0, len(custom_contacts[0]) - 1)
        if max_contacts >= len(custom_contacts[0]):
            return custom_contacts[0], custom_contacts[1]
        else:
            # Return one contact
            return custom_contacts[0][selected], custom_contacts[1][selected]
    
    else:
        log.info("Incorrect contact-type. Must be inter or intra.")
        return False


def update_distance_distribution_matrix(dist, d_mtx):
    """
    Update the distance and distribution matrices.

    Parameters
    ----------
    dist : np.ndarray of list
        Array of list of all possible contacts for given database entry.
    
    d_mtx : np.ndarray of dict
        Array of dict of number of counts for each distance tuple in the array.

    Returns
    -------
    d_mtx : np.ndarray of dict
        Updated distance_mtx from information in d_mtx.
    """
    # They should have the same shape
    assert np.shape(dist) == np.shape(d_mtx)

    dist = dist.astype(list)
    for x, row in enumerate(dist):
        for y, item in enumerate(row):
            db_entry = item
            existing = d_mtx[x, y]
            db_count = Counter(db_entry)
            for dist, count in db_count.items():
                existing[dist] = count
            d_mtx[x, y] = existing
    return d_mtx


def reverse_position_lookup(coords, location_mtx, database):
    """
    Return database entry based on a point in the contacts frequency heatmap.

    NOTE this will be heavily remodified given the updated way of database
    processing. Essentially this will take in all 3 aligned matrices and
    return the distance that we want.

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
