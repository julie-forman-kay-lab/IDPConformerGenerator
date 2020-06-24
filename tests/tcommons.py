"""Common funcs and variables for tests."""
import random
import shutil
from contextlib import suppress

from idpconfgen import Path



class DummyFolder:
    def __init__(self, folder):
        self.folder = Path(folder)

    def __enter__(self):
        self.folder.mkdir(parents=True, exist_ok=True)
        return self.folder

    def __exit__(self, *args, **kwargs):
        with suppress(FileNotFoundError):
            shutil.rmtree(self.folder)

class DummyFile:
    def __init__(self, fpath):
        self.fpath = Path(fpath)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        with suppress(FileNotFoundError):
            self.fpath.unlink()


# folders
tests_folder = Path(__file__).myparents()
project_folder = tests_folder.parent
iofiles_folder = Path(tests_folder, 'readfiles')
folder_output = Path(tests_folder, 'data_out_scratch')
data_folder = Path(tests_folder, 'data')
pdbids_folder = Path(data_folder, 'pdbids')

# files
PFA1 = Path(pdbids_folder, 'PFA1_A.pdb')
PFA2 = Path(pdbids_folder, 'PFA2_A.pdb')
PFA3 = Path(pdbids_folder, 'PFA3_A.pdb')
cif_EOF = Path(data_folder, 'cif_EOF.cif')
cif_example = Path(data_folder, 'cif_example.cif')
cif_example_auth = Path(data_folder, 'cif_example_auth.cif')
cif_example_noasymid = Path(data_folder, 'cif_example_noasymid.cif')
cif_noatomsite = Path(data_folder, 'cif_no_atom_site.cif')
cif_nohash = Path(data_folder, 'cif_nohash.cif')
cull = Path(data_folder, 'cull.list')
file_tar = Path(data_folder, 'files.tar')
pdb_example = Path(data_folder, 'pdb_example.pdb')
pdb_models = Path(data_folder, 'pdb_models.pdb')
pdb_models_output = Path(data_folder, 'pdb_models_output.pdb')
pdbs_fasta = Path(data_folder, 'pdbs.fasta')
dict1json = Path(data_folder, 'dict1.json')


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
