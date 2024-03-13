"""Functions for recognizing and processing multiple protein chains."""
from difflib import SequenceMatcher

import numpy as np

from idpconfgen import log
from idpconfgen.core.definitions import aa3to1
from idpconfgen.libs.libstructure import (
    col_chainID,
    col_resName,
    col_resSeq,
    col_segid,
    )
from idpconfgen.logger import S


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
    skipped_chains = []
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
        if max_match == 1.0:
            log.info(S(f"Identical sequence found, skipping chain {chain}."))
            fld_chainseq[chain] = fld_fasta
            skipped_chains += chain_lst
            continue
        chain_arr = np.array(chain_lst)
        fld_chainseq[chain] = (fld_fasta, match_index, chain_arr)
    
    return fld_chainseq, skipped_chains
