"""Common funcs and variables for tests."""
import random

from idpconfgen import Path


tests_folder = Path(__file__).myparents()

project_folder = tests_folder.parents[0]

data_folder = Path(tests_folder, 'data')

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
