"""Common funcs and variables for tests."""
import random

from idpconfgen import Path


tests_folder = Path(__file__).myparents()

project_folder = tests_folder.parents[0]
iofiles_folder = Path(tests_folder, 'readfiles')
data_folder = Path(tests_folder, 'data')

pdb_example = Path(data_folder, 'pdb_example.pdb')
pdb_models = Path(data_folder, 'pdb_models.pdb')
pdb_models_output = Path(data_folder, 'pdb_models_output.pdb')
cif_example = Path(data_folder, 'cif_example.cif')
cif_example_auth = Path(data_folder, 'cif_example_auth.cif')
cif_example_noasymid = Path(data_folder, 'cif_example_noasymid.cif')
cif_noatomsite = Path(data_folder, 'cif_no_atom_site.cif')
cif_nohash = Path(data_folder, 'cif_nohash.cif')
cif_EOF = Path(data_folder, 'cif_EOF.cif')

cif_example_headers = [
    '_atom_site.group_PDB',
    '_atom_site.id',
    '_atom_site.type_symbol',
    '_atom_site.label_atom_id',
    '_atom_site.label_alt_id',
    '_atom_site.label_comp_id',
    '_atom_site.label_asym_id',
    '_atom_site.label_entity_id',
    '_atom_site.label_seq_id',
    '_atom_site.pdbx_PDB_ins_code',
    '_atom_site.Cartn_x',
    '_atom_site.Cartn_y',
    '_atom_site.Cartn_z',
    '_atom_site.occupancy',
    '_atom_site.B_iso_or_equiv',
    '_atom_site.pdbx_formal_charge',
    '_atom_site.auth_seq_id',
    '_atom_site.auth_comp_id',
    '_atom_site.auth_asym_id',
    '_atom_site.auth_atom_id',
    '_atom_site.pdbx_PDB_model_num',
    ]

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
