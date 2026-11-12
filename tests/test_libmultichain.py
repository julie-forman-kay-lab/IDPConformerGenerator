"""Test processing of structures with several chains."""
import numpy as np

from idpconfgen.libs.libmultichain import process_multichain_pdb
from idpconfgen.libs.libstructure import Structure, col_chainID

from .tcommons import data_folder


def test_process_multichain_pdb_keeps_chain_order():
    """Test chains come in the order of the structure."""
    structure = Structure(data_folder / 'EXPL_A.pdb')
    structure.build()
    chains = []
    for chain in 'QBA':
        rows = np.array(structure.data_array, copy=True)
        rows[:, col_chainID] = chain
        chains.append(rows)
    rows = np.concatenate(chains)
    sequence = 'AYIAKQRQISFVKSHF'
    matched, skipped = process_multichain_pdb(rows, {'s': 'MM' + sequence})
    assert list(matched) == ['Q', 'B', 'A']
    assert skipped == []
    for chain, (fasta, index, chain_rows) in matched.items():
        assert fasta == sequence
        assert index == 0
        assert set(chain_rows[:, col_chainID]) == {chain}
    # chains with the sequence given are skipped, in order as well
    matched, skipped = process_multichain_pdb(rows, {'s': sequence})
    assert list(matched) == ['Q', 'B', 'A']
    assert all(value == sequence for value in matched.values())
    assert len(skipped) == len(rows)
