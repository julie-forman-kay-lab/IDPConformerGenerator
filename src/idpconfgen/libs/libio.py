"""Functions and Classes regarding Input/Output."""

import glob
import itertools as it
import json
import os
import pickle
import sys
import tarfile
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcheck
from idpconfgen.libs.libpdb import PDBList
from idpconfgen.logger import S, T


@libcheck.argstype(list)
# NOT USED ANYWHERE
def add_existent_files(storage, source):
    """
    Add files that exist to a list.

    Given a list of `source` Paths, if Path exist adds it to `storage`.

    Adds Path instances.
    """
    for path in source:
        p = Path(path).resolve()
        if p.is_file():
            storage.append(p)
        else:
            log.error(S('file not found: {}', p.str()))


def check_file_exist(files_list):
    """
    Confirm all files in a list exist.

    Logs each entry for each file not found.

    Parameters
    ----------
    files_list : list
        The list of file Paths or strings.

    Returns
    -------
    tuple of two elements
        True
            If all files are found. False if there is at least one file
            missing.

        list
            A list of the missing files.
    """
    log.info(T('checking files exist'))
    # I have avoided to use it.tee on purpose
    paths1 = map(Path, files_list)
    paths2 = map(Path, files_list)
    log_nonexistent(paths1)
    # list() because of the return statement
    files_not_found = list(it.filterfalse(Path.exists, paths2))
    return not bool(files_not_found), files_not_found


@libcheck.argstype(((list, tuple),))
def concatenate_entries(entry_list):
    """
    Concatente entries.

    Entries can be given in a list of entries or file paths with
    entry lists. Single entries in the input list are used directly
    while files are read and their lines added one by one to the
    concatenated list.

    Notice:
        Does not descriminate between single entries and mispelled
        file paths. Every string that cannot be openned as path
        is considered an individual entry.

    Parameters
    ----------
    entry_list : lits
        List containing strings or file paths

    Returns
    -------
    list
        Concatenated strings plus lines in files.
    """
    for entry in entry_list:
        try:
            with Path(entry).open('r') as fh:
                for i in filter(bool, map(str.strip, fh)):
                    yield i
        except (FileNotFoundError, NotADirectoryError):
            yield entry


def extract_from_tar(tar_file, output='.'):
    """
    Extract files from tarfile.

    Parameters
    ----------
    tar_file : str
        The tarfile to extract.

    output : str, optional
        The folder to where extract tar files.
        Defaults to current working directory.

    Returns
    -------
    list
        A list of Path-objects pointing to the extracted files.
    """
    tar = tarfile.open(tar_file)
    names = tar.getnames()
    tar.extractall(path=output)
    tar.close()
    return [Path(output, i) for i in names]


def glob_folder(folder, ext):
    """
    List files with extention `ext` in `folder`.

    Does NOT perform recursive search.

    Parameters
    ----------
    folder : str
        The path to the folder to investigate.

    ext : str
        The file extention. Can be with or without the dot [.]
        preffix.

    Returns
    -------
    list of Path objects
        SORTED list of matching results.
    """
    ext = f'*{period_suffix(ext)}'
    files = glob.glob(Path(folder, ext).str())
    log.debug(f'folder {folder} read {len(files)} files with extension {ext}')
    return list(map(Path, files))


def has_suffix(path, ext=None):
    """
    Evaluate file suffix according to `ext` condition.

    Parameters
    ----------
    path : str of Path
        The file path.

    ext : str
        The file extension.
        Can be dotted '.' (.csv) or plain (csv).

    Return
    ------
    bool
        True
            Always if `ext` is None.
            If path suffix equals `ext`.

        False
            Otherwise
    """
    PS = period_suffix
    return ext is None or Path(path).suffix == PS(ext)


def list_files_recursively(folder, ext=None):
    """
    List files recursively from source folder.

    Parameters
    ----------
    folder : string or Path
        The folder from where to start searching.

    ext : string
        The file extension to consider.
        Files without the defined `ext` will be ignored.
        Defaults to ``None``, all files are considered.

    Returns
    -------
    unsorted list
        Of the file paths relative to the source `folder`.
    """
    phassuff = partial(has_suffix, ext=ext)
    for root, _subdirs, files in os.walk(folder):
        only_ext = filter(phassuff, files)
        for file_ in only_ext:
            yield Path(root, file_).resolve()


def log_nonexistent(files):
    """
    Log to ERROR files that do not exist.

    Parameters
    ----------
    files : iterable of Paths
    """
    log.info(T('Logging non-existent files'))
    _not_exist = it.filterfalse(Path.exists, files)
    strings = (S('file not found: {}', path) for path in _not_exist)
    try:
        log.error('\n'.join(strings))
    except TypeError:
        log.info(S('all files exist'))
    log.info(S('done'))


def make_destination_folder(dest):
    """
    Make a destination folder.

    Returns
    -------
    Path-object
        A path pointing to the folder created.
        If ``dest`` is ``None`` returns a Path poiting to the CWD.
    """
    try:
        dest_ = Path(dest)
    except TypeError:
        return Path.cwd()
    else:
        dest_.mkdir(parents=True, exist_ok=True)
        return dest_


@libcheck.argstype(((list, tuple),))
def paths_from_flist(path):
    """
    Read Paths from file listing paths.

    Returns
    -------
    Map generator
        Path representation of the entries in the file.
    """
    with open(path, 'r') as fin:
        lines = fin.readlines()
    lines = map(str.strip, lines)
    valid = filter(bool, lines)
    no_comments = it.filterfalse(lambda x: x.startswith('#'), valid)
    return list(map(Path, no_comments))


def period_suffix(ext):
    """
    Represent a suffix of a file.

    Example
    -------
    period_suffix('.pdf')
    '.pdf'

    period_suffix('pdf')
    '.pdf'

    Parameters
    ----------
    ext : str
        String to extract the suffix from.

    Returns
    -------
    str
        File extension with leading period.
    """
    return f'.{ext[ext.find(".") + 1:]}'


def read_path_bundle(path_bundle, ext=None, listext='.list'):
    """
    Read path bundle.

    Read paths encoded in strings, folders or files that are list of files.

    If a string or Path object points to an existing file,
        register that path.

    If a string points to a folder, registers all files in that folder
        that have extension `ext`, recursively.

    If a string points to a file with extension `listext`, registers
        all files referred in the `listext` file and which exist in disk.

    Non-existing files are log as error messages.

    Parameters
    ----------
    path_bundle : list-like
        A list containing strings or paths that point to files or folders.

    ext : string
        The file extension to consider. If ``None`` considers all files.
        Defaults to ``None``.

    listext : string
        The file extension to consider as a file listing other files.
        Defaults to ``.flist``.

    Returns
    -------
    generator
        A generator that complies with the specifications described
    """
    # defines functions and partials
    LFR = list_files_recursively
    CFI = it.chain.from_iterable
    p_lfr_ext = partial(LFR, ext=ext)
    p_suff_list = partial(has_suffix, ext=listext)

    pbundle = map(Path, path_bundle)
    pbundle2 = map(Path, path_bundle)

    # files that do not exist
    log_nonexistent(pbundle2)

    # files that exist
    p1 = filter(Path.exists, pbundle)
    p2 = filter(Path.exists, pbundle)
    p3 = filter(Path.exists, pbundle)

    # is dir filter
    f1 = filter(Path.is_dir, p1)

    # reads list files
    f2 = filter(p_suff_list, p2)

    # files with valid extension
    def _is_f3(x):
        return all((
            x.is_file(),
            not has_suffix(x, ext=listext),
            has_suffix(x, ext=ext),
            ))

    f3 = filter(_is_f3, p3)

    from_folders = CFI(map(p_lfr_ext, f1))
    from_lists = CFI(map(paths_from_flist, f2))

    return it.chain(from_folders, from_lists, f3)


def read_dictionary_from_disk(path):
    """
    Reads a dictionary from disk.

    Accepted formats:
        * pickle
        * json

    Returns
    -------
    dict
    """
    _path_suffix = Path(path).suffix
    options = {
        #'.pickle': read_dict_from_pickle,
        '.json': read_dict_from_json,
        }

    the_dict = options[_path_suffix](path)

    assert isinstance(the_dict)
    return the_dict


def read_dict_from_json(path):
    """Read dict from json."""
    return json.loads(path)


def read_dict_from_pickle(path):
    """Read dictionary from pickle."""
    return pickle.load(open( path, "rb" ))


def read_PDBID_from_folder(folder):
    """
    Read PDBIDs from folder.

    Parameters
    ----------
    folder : str
        The folder to read.

    Returns
    -------
    :class:`idpconfgen.libs.libpdb.PDBList`
    """
    log.info(T('reading destination folder'))
    log.info(S(f'from: {folder}'))
    pdblist = PDBList(glob_folder(folder, '*.pdb'))
    log.info(
        f"{S(f'found: {str(pdblist)}')}\n"
        f"{S('done')}\n"
        )
    return pdblist


def read_PDBID_from_source(source):
    """
    Read PDBIDs from destination.

    Accepted destinations:
        * folder
        * tarfile

    Returns
    -------
    :class:`idpconfgen.libs.libpdb.PDBList`.
    """
    options = {
        source.suffix == '.tar': read_PDBID_from_tar,
        source.is_dir(): read_PDBID_from_folder,
        }
    return options[True](source)


def read_PDBID_from_tar(tar_file):
    """
    Read PDBIDs from tarfile.

    .. note::

        Case-specific function, not an abstraction.

    Parameters
    ----------
    tar_file : :class:`idpconfgen.Path`
        The tarfile to read.

    Returns
    -------
    :class:`idpconfgen.libs.libpdb.PDBList`
        If file is not found returns an empty PDBList.
    """
    title = (
        f"{T('Reading PDB IDs in tar file.')}\n"
        f"{S(f'from: {tar_file}')}\n"
        )
    log.info(title)

    try:
        tar = tarfile.open(tar_file.str(), 'r:*')
        p = PDBList(tar.getnames())
        tar.close()

    except FileNotFoundError:
        log.info(f'tar_file {tar_file} not found. Returning empty PDBList.')
        p = PDBList([])

    end = (
        f"{S(f'found: {str(tar_file)}')}"
        f"{S('done')}\n"
        )
    log.info(end)
    return p


def save_dictionary(mydict, output='mydict.pickle'):
    """
    Save dictionary to disk.

    Accepted formats:
        * pickle
        * json

    Parameters
    ----------
    mydict : dict
        The dict to be saved.

    output : str or Path
        The output file. Format is deduced from file extension.

    Raises
    ------
    KeyError
        If extension is not a compatible format.
    """
    options = {
        '.pickle': save_dict_to_pickle,
        '.json': save_dict_to_json,
        }

    options[Path(output).suffix](mydict, output)


def save_dict_to_json(
        mydict,
        output='mydict.json',
        indent=True,
        sort_keys=True,
        ):
    """Save dictionary to JSON."""
    assert Path(output).suffix == '.json'
    with open(output, 'w') as fout:
        json.dump(mydict, fout, indent=indent, sort_keys=sort_keys)


def save_dict_to_pickle(mydict, output='mydict.pickle'):
    """Save dictionary to pickle file."""
    assert Path(output).suffix == '.pickle'
    with open(output, 'wb') as handle:
        pickle.dump(mydict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def write_text(text, output=None):
    """
    Write text to output.

    If output is ``None`` writes to stdout.
    Else, writes to file.
    If ``output`` points to a directory that does not exist, creates it.
    """
    try:
        opath = Path(output)
    except TypeError:
        sys.stdout.write(text)
    else:
        opath.myparents().mkdir(parents=True, exist_ok=True)
        opath.write_text(text)
