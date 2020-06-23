"""Test I/O lib."""
import collections
import pytest
import types
from pprint import pprint
from idpconfgen import Path
from idpconfgen.libs import libio

from . import tcommons


def test_concatenate_0():
    """Test against cull.list."""
    result = list(libio.concatenate_entries([tcommons.cull]))
    expected = [
        '# 5XLI chains renamed to lowercase\n',
        '12E8H       221  XRAY        1.900    0.22    0.27\n',
        '16PKA       415  XRAY        1.600    0.19    0.23\n',
        '16VPA       366  XRAY        2.100    0.19    0.26\n',
        '1A04A       215  XRAY        2.200    0.21    0.27\n',
        '1A05A       358  XRAY        2.000    0.20    0.28\n',
        '1A0JA       223  XRAY        1.700    0.17    0.21\n',
        '1A12A       413  XRAY        1.700    0.19    0.22\n',
        '1A1XA       108  XRAY        2.000    0.21    0.25\n',
        ]


def test_concatenate_1():
    """Test concatenate entries."""
    user_input = [
        'ABC1D',
        'somefile_that_does_not_exist.list',
        str(Path(tcommons.iofiles_folder, 'file.list')),
        Path(tcommons.iofiles_folder, 'file.list'),
        ]

    expected_output = [
        'ABC1D',
        'somefile_that_does_not_exist.list',
        '# some comment line',
        'file6',
        'file7',
        'file8',
        '# some comment line',
        'file6',
        'file7',
        'file8',
        ]

    output = list(libio.concatenate_entries(user_input))

    assert expected_output == output


def test_paths_from_flist():
    result = libio.paths_from_flist(tcommons.iofiles_folder / 'file.list')
    assert isinstance(result, collections.Iterable)
    assert list(result) == [Path(f'file{i}') for i in range(6,9)]


@pytest.mark.parametrize('in1', ['string', 1, 1.0, Path('somepath')])
def test_concatenate_2(in1):
    """Test raises TypeError."""
    with pytest.raises(TypeError):
        list(libio.concatenate_entries(in1))


def test_check_file_exists_1():
    """Test with a file that exists."""
    result = libio.check_file_exist([tcommons.iofiles_folder / 'file1'])
    assert any(result)


def test_check_file_exists_2():
    """Test with a file that do not exist."""
    result = libio.check_file_exist([
        tcommons.iofiles_folder / 'file1',
        'donotexist',
        ])
    assert (False, [Path('donotexist')]) == result


@pytest.mark.parametrize('in1', [1, 1.0, 's', Path('p')])
def test_check_file_exists_typeerror(in1):
    """Test raises TypeError."""
    with pytest.raises(TypeError):
        libio.check_file_exist(in1)


@pytest.mark.parametrize(
    'input_',
    [
        tcommons.tests_folder / '__init__.py',
        (str(tcommons.tests_folder / '__init__.py')),
        ],
    )
def test_add_existent_files_1(input_):
    """Add a file that exists."""
    files = []
    libio.add_existent_files(files, [input_])
    assert files == [Path(input_)]


def test_add_existent_files_2():
    """Do not add a file that does not exist."""
    files = []
    tpath = Path(tcommons.tests_folder, 'doesnotexist.py')
    libio.add_existent_files(files, [tpath])
    assert not(files)


def test_has_suffix_1():
    """Test ext is None."""
    assert libio.has_suffix('file')


def test_has_suffix_2():
    """Test ext is equals file."""
    assert libio.has_suffix('file.pdb', ext='.pdb')


def test_has_suffix_3():
    """Test ext is differs file."""
    assert not libio.has_suffix('file.csv', ext='.pdb')


def test_has_suffix_4():
    """Test ext is None when file has ext."""
    assert libio.has_suffix('file.csv')


def test_has_suffix_5():
    """Test ext differs without dot."""
    assert not libio.has_suffix('file.csv', ext='pdb')


def test_has_suffix_6():
    """Test ext correct without dot."""
    assert libio.has_suffix('file.csv', ext='csv')


def test_list_recursively_1():
    """Test finds .ext files."""
    files = libio.list_files_recursively(
        tcommons.iofiles_folder,
        ext='.ext',
        )

    expected = [
        Path(tcommons.iofiles_folder, 'file.ext'),
        Path(tcommons.iofiles_folder, 'other_files', 'file2.ext'),
        ]

    assert list(files) == expected


def test_list_files_recursively_1():
    """Test None as  extention."""
    files = libio.list_files_recursively(tcommons.iofiles_folder)

    expected = [
        tcommons.iofiles_folder.joinpath(p) for p in [
            'file1',
            'file2',
            'file3',
            'file6',
            'file7',
            'file8',
            'file.ext',
            'file.list',
            'README',
            Path('other_files', 'file2.ext'),
            Path('other_files', 'file4'),
            Path('other_files', 'file5'),
            ]
        ]

    assert sorted(files) == sorted(expected)


@pytest.mark.parametrize(
    'ext',
    ['ext', '.ext'],
    )
def test_list_files_recursively_2(ext):
    """Test ext extention."""
    files = libio.list_files_recursively(
        tcommons.iofiles_folder,
        ext=ext,
        )

    expected = [
        tcommons.iofiles_folder.joinpath(p) for p in [
            Path('other_files', 'file2.ext'),
            'file.ext',
            ]
        ]

    assert sorted(files) == sorted(expected)


@pytest.mark.parametrize('in1', ['a', Path('a')])
def test_read_path_bundle_typeerror(in1):
    """Test read_path_bundle raises TypeError."""
    with pytest.raises(TypeError):
        libio.read_path_bundle(in1)


@pytest.mark.parametrize(
    'in1,ext,expected',
    [
        (
        [
            tcommons.iofiles_folder / 'file1',
            tcommons.iofiles_folder / 'other_files',
            tcommons.iofiles_folder / 'file.list',
            ],
        None,
        [
            tcommons.iofiles_folder / 'file1',
            Path('file6'),
            Path('file7'),
            Path('file8'),
            tcommons.iofiles_folder / 'other_files' / 'file2.ext',
            tcommons.iofiles_folder / 'other_files' / 'file4',
            tcommons.iofiles_folder / 'other_files' / 'file5',
            ]
        )]
    )
def test_read_bundle_inputs(in1, ext, expected):
    """Test read_bundle multiple inputs."""
    result = libio.read_path_bundle(in1, ext=ext, listext='.list')
    assert sorted(expected) == sorted(result)


@pytest.mark.parametrize(
    'in1,ext,expected',
    [
        (
            tcommons.iofiles_folder,
            '.ext',
            [tcommons.iofiles_folder / 'file.ext'],
            ),
        (
            tcommons.iofiles_folder / 'other_files',
            'ext',
            [tcommons.iofiles_folder / 'other_files' / 'file2.ext'],
            ),
        ]
    )
def test_glob_folder(in1, ext, expected):
    """Test glob folder function."""
    assert expected == libio.glob_folder(in1, ext)
