"""Common funcs and variables for tests."""

from idpconfgen import Path


tests_folder = Path(__file__).myparents()

data_folder = Path(tests_folder, 'data')

folder_output = Path(
    tests_folder,
    'data_out_scratch',
    )
