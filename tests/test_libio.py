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


@pytest.mark.parametrize('in1', ['string', 1, 1.0, Path('somepath')])
def test_concatenate_2(in1):
    """Test raises TypeError."""
    with pytest.raises(TypeError):
        libio.concatenate_entries(in1)


def test_check_file_exists_1():
    """Test with a file that exists."""
    result = libio.check_file_exist([tcommons.data_folder / '1A12_A.pdb'])
    assert (True, []) == result


def test_check_file_exists_2():
    """Test with a file that do not exist."""
    result = libio.check_file_exist([
        tcommons.data_folder / '1A12_A.pdb',
        'donotexist',
        ])
    assert (False, ['donotexist']) == result


@pytest.mark.parametrize('in1', [1, 1.0, 's', Path('p')])
def test_check_file_exists_typeerror(in1):
    """Test raises TypeError."""
    with pytest.raises(TypeError):
        libio.check_file_exist(in1)


@pytest.mark.parametrize(
    'input_',
    [
        tcommons.tests_folder / '__init__.py',
        (tcommons.tests_folder / '__init__.py').str(),
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
    """Test finds .pdb files."""
    files = libio.list_files_recursively(
        tcommons.data_folder,
        ext='.pdb',
        )

    assert files == [Path(tcommons.data_folder, '1A12_A.pdb')]


def test_list_files_recursively_1():
    """Test None as  extention."""
    files = libio.list_files_recursively(
        tcommons.data_folder,
        )
    
    expected = [
        tcommons.data_folder.joinpath(p) for p in [
            '1ABC_D.dssp',
            '1ABC_E.dssp',
            'wrong.dssp',
            'wrong2.dssp',
            '1A12_A.pdb',
            'cull.list',
            'pdblist.list',
            'path_bundle.flist',
            ]
        ]
    
    assert sorted(files) == sorted(expected)


def test_list_files_recursively_2():
    """Test pdb extention."""
    files = libio.list_files_recursively(
        tcommons.data_folder,
        ext='pdb',
        )
    
    expected = [tcommons.data_folder.joinpath(p) for p in ['1A12_A.pdb']]
    
    assert sorted(files) == sorted(expected)


def test_list_files_recursively_3():
    """Test .pdb."""
    files = libio.list_files_recursively(
        tcommons.data_folder,
        ext='.pdb',
        )
    expected = [tcommons.data_folder.joinpath(p) for p in ['pdb_example.pdb', 'pdb_saved.pdb']]
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
            [tcommons.data_folder],
            None,
            [
                Path(tcommons.data_folder, '1ABC_D.dssp'),
                Path(tcommons.data_folder, '1ABC_E.dssp'),
                Path(tcommons.data_folder, 'cull.list'),
                Path(tcommons.data_folder, 'pdblist.list'),
                Path(tcommons.data_folder, 'pdb_example.pdb'),
                Path(tcommons.data_folder, 'pdb_saved.pdb'),
                Path(tcommons.data_folder, 'wrong.dssp'),
                Path(tcommons.data_folder, 'wrong2.dssp'),
                ]
            ),
        (
            [
                tcommons.data_folder,
                Path(tcommons.data_folder, 'path_bundle.flist'),
                Path(tcommons.data_folder, 'noexist_bundle.flist'),
                ],
            None,
            [
                Path(tcommons.project_folder, 'setup.py'),
                Path(tcommons.data_folder, '1ABC_D.dssp'),
                Path(tcommons.data_folder, '1ABC_E.dssp'),
                Path(tcommons.data_folder, 'cull.list'),
                Path(tcommons.data_folder, 'pdblist.list'),
                Path(tcommons.data_folder, 'wrong.dssp'),
                Path(tcommons.data_folder, 'wrong2.dssp'),
                ]
            ),
        (
            [Path(tcommons.data_folder, 'path_bundle.flist')],
            None,
            [
                Path(tcommons.project_folder, 'setup.py'),
                ]
            ),
        (
            [tcommons.data_folder],
            '.pdb',
            [Path(tcommons.data_folder, '1A12_A.pdb')],
            ),
        (
            [tcommons.data_folder],
            'pdb',
            [Path(tcommons.data_folder, '1A12_A.pdb')],
            ),
        (
            [tcommons.data_folder],
            '.none',
            [],
            ),
        ]
    )
def test_read_bundle_inputs(in1, ext, expected):
    """Test read_bundle multiple inputs."""
    assert expected == libio.read_path_bundle(in1, ext=ext)


@pytest.mark.parametrize(
    'in1,ext,expected',
    [
        (
            tcommons.data_folder,
            '.pdb',
            [
                Path(tcommons.data_folder, 'pdb_example.pdb'),
                Path(tcommons.data_folder, 'pdb_saved.pdb')
                ],
            ),
        (
            tcommons.data_folder,
            'pdb',
            [
                Path(tcommons.data_folder, 'pdb_example.pdb'),
                Path(tcommons.data_folder, 'pdb_saved.pdb')
                ],
            ),
        (
            tcommons.data_folder,
            '.dssp',
            [
                Path(tcommons.data_folder, '1ABC_D.dssp'),
                Path(tcommons.data_folder, '1ABC_E.dssp'),
                Path(tcommons.data_folder, 'wrong.dssp'),
                Path(tcommons.data_folder, 'wrong2.dssp'),
                ],
            ),
        (tcommons.data_folder, '.none', []),
        ]
    )
def test_glob_folder(in1, ext, expected):
    """Test glob folder function."""
    assert expected == libio.glob_folder(in1, ext)
