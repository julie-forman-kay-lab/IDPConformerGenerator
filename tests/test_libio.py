"""Test I/O lib."""
import collections
import glob
import os
import tarfile
from multiprocessing import Manager

import pytest

from idpconfgen import Path
from idpconfgen.core.definitions import XmerProbs
from idpconfgen.libs import libio
from idpconfgen.libs.libpdb import PDBList

from . import tcommons


def test_concatenate_0():
    """Test against cull.list."""
    result = list(libio.concatenate_entries([tcommons.cull]))
    expected = [
        '# 5XLI chains renamed to lowercase',
        '12E8H       221  XRAY        1.900    0.22    0.27',
        '16PKA       415  XRAY        1.600    0.19    0.23',
        '16VPA       366  XRAY        2.100    0.19    0.26',
        '1A04A       215  XRAY        2.200    0.21    0.27',
        '6XY7AAA',
        ]
    assert result == expected


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


@pytest.mark.parametrize(
    'ext,num',
    [
        ('.txt', 2),
        ('.pdb', 1),
        (None, 3),
        ]
    )
def test_extract_from_tar(ext, num):
    """Test extract files from tarfile."""
    with tcommons.TmpFolder('tar_extracted') as df:
        result = libio.extract_from_tar(tcommons.file_tar, output=df, ext=ext)
        assert all(isinstance(i, Path) for i in result)
        suf = f'*{ext}' if ext else '*.*'
        num_files = glob.glob(str(Path(df, suf)))
        assert len(num_files) == num


def test_log_nonexistent_1():
    """Log non existent files."""
    libio.log_nonexistent([tcommons.file_tar])


def test_log_nonexistent_2():
    """Log non existent files."""
    libio.log_nonexistent([Path('doesnotexist')])


@pytest.mark.parametrize(
    'folder',
    [
        Path('this_folder'),
        Path('this', 'folder'),
        ]
    )
def test_make_destination_folder_1(folder):
    """Test make dest folder."""
    rfolder = libio.make_destination_folder(folder)
    assert rfolder.exists()
    assert rfolder == folder
    rfolder.rmdir()


def test_make_destination_folder_2():
    """Test make dest folder."""
    folder = libio.make_destination_folder(None)
    assert folder.exists()
    assert folder == Path.cwd()


def test_read_dictionary_from_disk():
    """Read dictionary from disk."""
    result = libio.read_dictionary_from_disk(tcommons.dict1json)
    assert result == {'somestring': 'some value'}


def test_read_dictionary_from_disk_tar():
    """Read dictionary from disk in tar."""
    result = libio.read_dictionary_from_disk(tcommons.dict1tar)
    assert result == {'somestring': 'some value'}


def test_read_dictionary_from_disk_tar2():
    """Read dictionary from disk in tar."""
    result = libio.read_dict_from_tar(tcommons.dict1tar)
    assert result == {'somestring': 'some value'}


def test_save_read_from_pickle():
    """Save and read to pickle."""
    d = {'somekey': 'some value'}
    o = 'delme.pickle'
    libio.save_dict_to_pickle(d, output=o)
    result = libio.read_dict_from_pickle(o)
    assert result == d
    Path(o).unlink()


def test_readPDBIDs_from_folder():
    """Read PDBIDs from folder."""
    result = libio.read_PDBID_from_folder(tcommons.pdbids_folder)
    assert isinstance(result, PDBList)


def test_read_PDBID_from_tar():
    """Read PDBIds from tar."""
    result = libio.read_PDBID_from_tar(tcommons.file_tar)
    expected = PDBList(['XXXX.pdb'])
    assert result == expected


def test_read_PDBID_from_tar_2():
    """Read from non-existent."""
    result = libio.read_PDBID_from_tar('nothing.tar')
    expected = PDBList([])
    assert expected == result


@pytest.mark.parametrize(
    'source',
    [
        tcommons.file_tar,
        tcommons.pdbids_folder,
        ]
    )
def test_read_PDBID_from_source(source):
    """Read PDBID from source."""
    r = libio.read_PDBID_from_source(source)
    assert isinstance(r, PDBList)


@pytest.mark.parametrize(
    'dest',
    [
        Path('mydisk.json'),
        Path('mydict.pickle'),
        ]
    )
def test_save_dictionary(dest):
    """Save dictionary to disk."""
    d = {'key': 'value'}
    libio.save_dictionary(d, output=dest)
    assert dest.exists()
    dest.unlink()


# capsys, crazy pytest undeclared variable! :-o
# https://docs.pytest.org/en/stable/reference.html?highlight=capsys#capsys
def test_save_dictionary_stdout(capsys):
    """Capture output."""
    d = {'key': 'value'}
    libio.save_dictionary(d, output=None)
    captured = capsys.readouterr()
    assert captured.out == f'{str(d)}\n'


def test_save_dict_to_json():
    """Test with TypeError."""
    m = Manager()
    d = m.dict()
    d['key'] = 'value'
    libio.save_dict_to_json(d)
    p = Path('mydict.json')
    assert p.exists()
    p.unlink()


def test_save_dict_to_pickle_with_manager():
    """Test with TypeError."""
    m = Manager()
    d = m.dict()
    d['key'] = 'value'
    libio.save_dict_to_pickle(d)
    p = Path('mydict.pickle')
    assert p.exists()
    p.unlink()


@pytest.mark.parametrize(
    'func',
    [
        tcommons.PFA1.read_text,
        tcommons.PFA1.read_bytes,
        ]
    )
def test_save_file_to_tar_1(func):
    """Test file is saved to tar correctly."""
    data = func()
    ftar = 'dummy.tar'
    with tarfile.open(ftar, 'w') as tarout:
        libio.save_file_to_tar(tarout, 'pfa1.pdb', data)

    with tarfile.open(ftar, 'r') as tarin:
        m = tarin.getmembers()
        names = tarin.getnames()
        d2 = tarin.extractfile(m[0]).read()

    assert len(names) == 1
    assert names[0] == 'pfa1.pdb'
    try:
        assert d2.decode() == data
    except AssertionError:
        assert d2 == data
    Path(ftar).unlink()


@pytest.mark.parametrize(
    'func',
    [
        libio.save_pairs_to_files,
        libio.save_pairs_to_disk,
        ]
    )
def test_save_pairs_to_disk(func):
    """Save pairs to disk."""
    df = 'save_pairs_to_disk_folder'
    pairs = dict([
        ('pair1.txt', 'some data'),
        ('pair2.txt', b'other data'),
        ])
    with tcommons.TmpFolder(df):
        func(pairs.items(), df)
        rfiles = [Path(p).name for p in glob.glob(str(Path(df, '*.txt')))]
        assert sorted(pairs.keys()) == sorted(rfiles)
        for fname, data in pairs.items():
            with open(Path(df, fname)) as fin:
                txt = fin.read()
                try:
                    assert data == txt
                except AssertionError:
                    assert data == txt.encode()


@pytest.mark.parametrize(
    'func',
    [
        libio.save_pairs_to_tar,
        libio.save_pairs_to_disk,
        ]
    )
def test_save_pairs_to_tar(func):
    """Save pairs to tar."""
    df = 'saved_pairs.tar'
    pairs = dict([
        ('pair1.txt', 'some data'),
        ('pair2.txt', b'other data'),
        ])
    with tcommons.TmpFile(df):
        func(pairs.items(), df)
        with tarfile.open(df, 'r') as tin:
            names = tin.getnames()
            for name, data in pairs.items():
                assert name in names
                txt = tin.extractfile(name).read()
                try:
                    assert data == txt
                except AssertionError:
                    assert data == txt.decode()


def test_paths_from_flist():
    """Test paths from list."""
    result = libio.paths_from_flist(tcommons.iofiles_folder / 'file.list')
    assert isinstance(result, collections.Iterable)
    assert list(result) == [Path(f'file{i}') for i in range(6, 9)]


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


@pytest.mark.parametrize(
    'in1,ext,expected',
    [(
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


@pytest.fixture
def cifs():
    """CIF files."""
    return [
        Path(tcommons.data_folder, p).resolve()
        for p in os.listdir(tcommons.data_folder)
        if Path(tcommons.data_folder, p).suffix == '.cif'
        ]


@pytest.mark.parametrize(
    'Iterator',
    [
        libio.FileIterator,
        libio.FileReaderIterator,
        ]
    )
def test_file_reader_iterator(Iterator, cifs):
    """FileReaderIterator."""
    assert len(cifs) > 0  # on purpose
    fi = Iterator(cifs, ext='.cif')
    assert len(fi) == len(cifs)
    for (n, f), c in zip(fi, cifs):
        assert n == c
        with open(c, 'rb') as ci:
            assert f == ci.read()


def test_FileIterator(cifs):
    """Test __getitem__."""
    lc = len(cifs)
    assert lc > 0
    fi = libio.FileIterator(cifs, ext='.cif')
    for i in range(0, len(cifs), 2):
        S = slice(i, i + 2)
        for (n, f), c in zip(fi[S], cifs[S]):
            assert n == c
            with open(c, 'rb') as ci:
                assert f == ci.read()


def test_FileIterator_all(cifs):
    """Test __getitem__."""
    lc = len(cifs)
    assert lc > 0
    fi = libio.FileIterator(cifs, ext='.cif')
    for (n, f), c in zip(fi[::], cifs):
        assert n == c
        with open(c, 'rb') as ci:
            assert f == ci.read()


@pytest.mark.parametrize(
    'Iterator',
    [
        libio.TarFileIterator,
        libio.FileReaderIterator,
        ]
    )
def test_file_reader_iterator_tar(Iterator):
    """FileReaderIterator."""
    pairs = dict([
        ('fileintar1.txt', 'data in tar\n'),
        ('fileintar2.txt', 'data file in 2\n'),
        ])
    fi = Iterator(str(tcommons.file_tar), ext='.txt')
    assert fi
    for n, data in fi:
        assert data.decode() == pairs[n]


def test_is_valid_fasta_file():
    """Test is valid fasta file."""
    assert libio.is_valid_fasta_file(tcommons.fasta1)
