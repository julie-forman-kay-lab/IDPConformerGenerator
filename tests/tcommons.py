"""Common funcs and variables for tests."""

from idpconfgen import Path


tests_folder = Path(__file__).myparents()

project_folder = tests_folder.parents[0]

data_folder = Path(tests_folder, 'data')

folder_output = Path(
    tests_folder,
    'data_out_scratch',
    )
