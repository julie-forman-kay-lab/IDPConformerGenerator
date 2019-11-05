"""Functions and Classes regarding Input/Output."""
from idpconfgen import Path, log
from idpconfgen.libs import libcheck


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


def read_paths(path_bundle, ext=None):
    
    files = []

    for path in path_bundle:
        
        p = Path(path)

        if p.is_dir():
            files.extend(read_files_recursively(p, ext=ext))
        
        elif p.suffix == ext:
            files.append(p)

        elif p.is_file():
            files.append(p)
