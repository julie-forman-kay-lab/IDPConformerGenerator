"""Functions for recognizing and processing multiple protein chains."""
from difflib import SequenceMatcher

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
                for k in range(i - 2, i + 3):
                    if 0 <= k < num_residues:
                        d_ca = np.linalg.norm(np.array(ca_coordinates[k]) - np.array(ca_coordinates[j]))  # noqa: E501
                        ca_dists.append(d_ca)
                        chain1_seq += fasta[k]

                for k in range(j - 2, j + 3):
                    if 0 <= k < num_residues:
                        chain2_seq += fasta[k]

                contacts.append({chain2_seq: [chain1_seq, ca_dists]})
                counter += 1
    
    if counter == 0:
        return False

    return pdbid, contacts, counter


def calc_interchain_ca_contacts(pdb, max_dist):
    """
    Find CA contacts below a certain distance between different chains.
    
    Returns contact information based on specific chains.
    Has an additional datapoint compared to intrachain ca contacts.
    
    Format of return:
    [
        {
            "SEQUENCE":[
                "pdbid_B_seg0",
                "SEQUENCE",
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
    chain_contacts = []
    counter = 0
    pdb_ids = []
    
    # TODO: when making comparisons, don't make comparisons
    # between different segments of the same chain.
    # We need to keep things between unique chains.
    
    if counter == 0:
        return False
    
    return pdb_ids, chain_contacts, counter
