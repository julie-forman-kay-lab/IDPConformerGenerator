"""Functions and Classes regarding Input/Output."""
import glob
import os
import sys
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

    files_not_found = list(filter(
        lambda x: not Path(x).exists(),
        files_list))
    
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
            files_list.append(Path(root, file_))
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
    func = add_existent_files
    hsuffix = has_suffix

    files = []
    
    folders = filter(
        lambda x: Path(x).is_dir(),
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
                possible_files = [l.strip() for l in fh.readlines()]
        except FileNotFoundError:
            log.error(S('file not found: {}', path))
            continue
        
        files_w_ext = filter(
            partial(hsuffix, ext=ext),
            possible_files,
            )

        func(files, files_w_ext)
    
    return sorted(p for p in files if p.suffix != listext)


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
    ext = f"*.{ext.strip().lstrip('*').lstrip('.')}"
    files = sorted(glob.glob(Path(folder, ext).str()))
    log.debug(f'folder {folder} read {len(files)} files with extension {ext}')
    return [Path(p) for p in files]


def write_text(text, output=None):
    """
    Writes text to output.

    If output is ``None`` writes to stdout.
    Else, writes to file.
    """
    try:
        opath = Path(output)
    except TypeError:
        sys.stdout.write(text)
    else:
        opath.myparents().mkdir(parents=True, exist_ok=True)
        opath.write_text(text)
