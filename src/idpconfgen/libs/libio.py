"""Functions and Classes regarding Input/Output."""
import itertools as it
import glob
import os
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcheck
from idpconfgen.logger import S, T


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
    concatenated = []
    for entry in entry_list:
        try:
            with Path(entry).open('r') as fh:
                concatenated.extend(fh.readlines())
        except (FileNotFoundError, NotADirectoryError):
            concatenated.append(entry)

    return concatenated


@libcheck.argstype(list)
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
    True
        If all files are found. False if there is at least one file
        missing.

    list
        A list of the missing files.
    """
    log.info(T('checking files exist'))

    # list() because of the return statement
    files_not_found = list(filter(
        lambda x: not Path(x).is_file(),
        files_list,
        ))

    for file_ in files_not_found:
        log.info(S('File NOT found: {}', file_))

    return not bool(files_not_found), files_not_found


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
    True
        Always if `ext` is None.
        If path suffix equals `ext`.

    False
        Otherwise

    """
    return ext is None or Path(path).suffix == '.{}'.format(ext.lstrip('.'))


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
    files_list = []
    for root, _subdirs, files in os.walk(folder):

        only_ext = filter(
            partial(has_suffix, ext=ext),
            files,
            )

        for file_ in only_ext:
            files_list.append(Path(root, file_).resolve())

    return files_list


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


def read_from_list(path):
    with open(path, 'r') as fin:
        lines = map(str.strip, fin.readlines())
        valid = filter(bool, lines)
        no_comments = filter(lambda x: not x.startswith('#'), valid)
        return list(map(Path, no_comments))


@libcheck.argstype((list, tuple),)
def read_path_bundle(path_bundle, ext=None, listext='.flist'):
    """
    Read path bundle.

    Read paths encoded in strings, folders or files that are list of files.

    If a string or Path object points to an existing file,
        register that path.

    If a string points to a folder, registers all files in that folder
        that have extension `ext`, recursively.

    If a string points to a file with extension `listext`, registers
        all files refered in the `listext` file that exist in disk.

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
    list
        A sorted list of all the files registered that exist in the disk.
    """
    files = []

    LFR = list_files_recursively
    CFI = it.chain.from_iterable

    pbundle, pbundle2 = it.tee(map(Path, path_bundle))

    _exist = filter(Path.exists, pbundle)
    _not_exist = filter(lambda x: not x.exists(), pbundle2)

    p1, p2, p3 = it.tee(_exist, 3)

    # is dir filter
    f1 = filter(Path.is_dir, p1)

    # reads list files
    partial2 = partial(has_suffix, ext=listext)
    f2 = filter(partial2, p2)

    # files with valid extension
    def _is_f3(x):
       return x.is_file() \
            and not has_suffix(x, ext=listext) \
            and has_suffix(x, ext=ext)

    f3 = filter(_is_f3, p3)

    # appends everything
    files.extend(CFI(map(partial(LFR, ext=ext), f1)))
    files.extend(CFI(map(read_from_list, f2)))
    files.extend(f3)

    [log.error(S('file not found: {}', path)) for path in _not_exist]

    return sorted(files)


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
    ext = f'*.{ext[ext.find(".") + 1:]}'
    files = sorted(glob.glob(Path(folder, ext).str()))
    log.debug(f'folder {folder} read {len(files)} files with extension {ext}')
    return [Path(p) for p in files]
