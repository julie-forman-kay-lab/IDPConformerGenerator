"""Functions and Classes regarding Input/Output."""
from idpconfgen import Path
from idpconfgen.libs import libcheck


@libcheck.checkargtype(((list, tuple),))
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
