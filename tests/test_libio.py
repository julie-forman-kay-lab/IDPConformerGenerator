"""Test I/O lib."""
import pytest

from idpconfgen import Path
from idpconfgen.libs import libio

from . import tcommons


def test_concatenate_1():
    """Test concatenate entries."""
    user_input = [
        'IDDXCH',
        'somefile_that_does_not_exist.list',
        Path(tcommons.data_folder, 'pdblist.list').str(),
        Path(tcommons.data_folder, 'pdblist.list'),
        ]

    expected_output = [
        'IDDXCH',
        'somefile_that_does_not_exist.list',
        '123A\n',
        '123AB\n',
        '123ABC\n',
        '123A90\n',
        '# some comment line\n',
        '123A\n',
        '123AB\n',
        '123ABC\n',
        '123A90\n',
        '# some comment line\n',
        ]

    output = libio.concatenate_entries(user_input)

    assert expected_output == output


def test_concatenate_typerror_path():
    """Raise TypeError when Path is given."""
    with pytest.raises(TypeError):
        libio.concatenate_entries(Path('somepath'))


def test_concatenate_typerror_str():
    """Raise TypeError when string is given."""
    with pytest.raises(TypeError):
        libio.concatenate_entries('somestr')


def test_add_existent_files_1():
    """Add a file that exists."""
    l = []
    tpath = tcommons.tests_folder / '__init__.py'
    libio.add_existent_files(
        l,
        [tpath]
        )

    assert l == [tpath]


def test_add_existent_files_2():
    """Do not add a file that does not exist."""
    l = []
    tpath = tcommons.tests_folder / 'doesnotexist.py'
    libio.add_existent_files(
        l,
        [tpath]
        )

    assert not(l)


def test_has_suffix_1():
    """Test ext is None."""
    assert libio.has_suffix('file') is True


def test_has_suffix_2():
    """Test ext is equals file."""
    assert libio.has_suffix('file.pdb', ext='.pdb') is True


def test_has_suffix_3():
    """Test ext is differs file."""
    assert libio.has_suffix('file.csv', ext='.pdb') is False


def test_has_suffix_4():
    """Test ext is None when file has ext."""
    assert libio.has_suffix('file.csv') is True


def test_has_suffix_5():
    """Test ext differs without dot."""
    assert libio.has_suffix('file.csv', ext='pdb') is False


def test_has_suffix_6():
    """Test ext correct without dot."""
    assert libio.has_suffix('file.csv', ext='csv') is True


def test_list_recursively_1():
    files = libio.list_files_recursively(
        tcommons.data_folder,
        ext='.pdb',
        )

    assert files == [tcommons.data_folder / '1A12_A.pdb']


def test_list_files_recursively_1():
    """Test None extention."""
    files = libio.list_files_recursively(
        tcommons.data_folder,
        )
    
    expected = [tcommons.data_folder.joinpath(p) for p in [
            '1A12_A.dssp',
            '1A12_A.pdb',
            'cull.list',
            'pdblist.list',
            ]
            ]
    
    assert sorted(files) == sorted(expected)


def test_list_files_recursively_2():
    """Test pdb extention."""
    files = libio.list_files_recursively(
        tcommons.data_folder,
        ext='pdb',
        )
    
    expected = [tcommons.data_folder.joinpath(p) for p in [
            '1A12_A.pdb',
            ]
            ]
    
    assert sorted(files) == sorted(expected)


def test_list_files_recursively_3():
    """Test .pdb."""
    files = libio.list_files_recursively(
        tcommons.data_folder,
        ext='.pdb',
        )
    
    expected = [tcommons.data_folder.joinpath(p) for p in [
            '1A12_A.pdb',
            ]
            ]

    assert sorted(files) == sorted(expected)


def test_read_paths1():
    """Test single Path in list."""
    paths = libio.read_path_bundle([
        Path(tcommons.data_folder, '1A12_A.pdb'),
        ])

    assert paths == [Path(tcommons.data_folder, '1A12_A.pdb')]


def test_read_paths2():
    """Read PDBs in folder."""

    paths = libio.read_path_bundle([tcommons.data_folder], ext='.pdb')

    assert paths == [Path(tcommons.data_folder, '1A12_A.pdb')]
