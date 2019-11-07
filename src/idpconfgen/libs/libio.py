"""Functions and Classes regarding Input/Output."""
import os
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcheck
from idpconfgen.logger import S


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
    pdblist : lits
        List containing PDBID:CHAIN identifiers of file list.

    Returns
    -------
    list
        Concatenated entries.
    """
    concatenated = []
    for entry in entry_list:
        try:
            with Path(entry).open('r') as fh:
                concatenated.extend(fh.readlines())
        except (FileNotFoundError, NotADirectoryError):
            concatenated.append(entry)

    return concatenated


def check_file_exist(files_list):
    """
    Confirms all files in a list exist.
    
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

    files_not_found = list(filter(
        lambda x: Path(x).exists(),
        files_list))
    
    for file_ in files_not_found:
        log.info(S('File NOT found: {}', file_))
    
    return not bool(files_not_found), files_not_found


def read_paths(path_bundle):
    """
    Reads paths from string bundle.

    String bundle can be list of strings or strings in files.

    Tandens the usage of :func:`concatenate_entries` and
    :func:`check_file_exists`.

    Returns the result from :func:`check_file_exists`.
    """
    list_of_paths = concatenate_entries(path_bundle)
    return check_file_exist(list_of_paths)


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
    try:
        for root, subdirs, files in os.walk(folder):
            
            only_ext = filter(
                partial(has_suffix, ext=ext),
                files,
                )
            
            for file_ in only_ext:
                files_list.append(Path(root, file_))
    except TypeError:
        pass
    return files_list


def add_existent_files(storage, source):
    """
    Add files that exist to a list.

    Given a list of `source` Paths, if Path exist adds it to `storage`.
    """
    for path in source:
        if path.is_file():
            storage.append(path)
        else:
            log.error(S('file not found: {}', path.str()))


def read_path_bundle(path_bundle, ext=None, listext='.list'):
    """
    Read path bundle.

    Read paths encoded in strings, folders or files that are list of files.

    If a string points to an existing file, register that path.

    If a string points to a folder, registers all files in that folder
        that have extension `ext`, recursively.
    
    If a string points to a file with extension `listext`, registers
        all files refered in the `listext` file that exist in disk.
    
    Parameters
    ----------
    path_bundle : list-like
        A list containing strings or paths that point to files or folders.
    
    ext : string
        The file extension to consider. If ``None`` considers any file.
        Defaults to ``None``.

    listext : string
        The file extension to consider as a file listing other files.
        Defaults to ``.list``.

    Returns
    -------
    list
        A list of all the files registered that exist in the disk.
    """
   
    func = add_existent_files
    hsuffix = has_suffix
    #path_bundle = list(map(Path, path_bundle))

    files = []
    
    folders = filter(
        lambda x: x.is_dir(),
        path_bundle,
        )

    extfiles = filter(
        partial(hsuffix, ext=ext),
        path_bundle,
        )

    listfiles = filter(
        partial(hsuffix, ext=listext),
        path_bundle,
        )
    
    for folder in folders:
        files.extend(list_files_recursively(folder, ext=ext))

    func(files, extfiles)
   
    for path in listfiles:
        try:
            with path.open('r') as fh:
                possible_files = fh.readlines()
        except FileNotFoundError:
            log.error(S('file not found: {}', path))
            continue
        
        files_w_ext = filter(
            partial(hsuffix, ext=ext),
            possible_files,
            )

        func(files, files_w_ext)
    
    return files


