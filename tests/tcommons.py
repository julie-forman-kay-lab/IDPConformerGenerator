"""Common funcs and variables for tests."""
import random

from idpconfgen import Path


tests_folder = Path(__file__).myparents()

project_folder = tests_folder.parents[0]

data_folder = Path(tests_folder, 'data')
pdb_example = Path(data_folder, 'pdb_example.pdb')
cif_example = Path(data_folder, 'cif_example.cif')
cif_noatomsite = Path(data_folder, 'cif_no_atom_site.cif')
cif_nohash = Path(data_folder, 'cif_nohash.cif')
cif_EOF = Path(data_folder, 'cif_EOF.cif')

pdb_saved = Path(data_folder, 'pdb_saved.pdb')

folder_output = Path(
    tests_folder,
    'data_out_scratch',
    )


def random_type():
    types = [
        1,
        1.0,
        [1,2],
        {'a': 1},
        None,
        {},
        [],
        set(),
        ]
    return random.choice(types)
