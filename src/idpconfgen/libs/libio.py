"""Functions and Classes regarding Input/Output."""
import itertools as it
import glob
import json
import os
import sys
import pickle
from functools import partial
import tarfile


from idpconfgen import Path, log
from idpconfgen.libs import libcheck
from idpconfgen.libs.libpdb import PDBList
from idpconfgen.logger import S, T


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
                concatenated.extend(
                    filter(
                        bool,
                        map(str.strip, fh.readlines()),
                        ),
                    )
        except (FileNotFoundError, NotADirectoryError):
            concatenated.append(entry)

    return concatenated


@libcheck.argstype(((list, tuple),))
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

    pbundle, pbundle2 = it.tee(map(Path, path_bundle))

    # files that do not exist
    log_nonexistent(pbundle2)

    # files that exist
    _exist = filter(Path.exists, pbundle)
    p1, p2, p3 = it.tee(_exist, 3)

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


def log_nonexistent(files):
    """
    Log to ERROR files that do not exist.

    Parameters
    ----------
    files : iterable of Paths
    """
    log.info(T('Logging non-existent files'))
    _not_exist = it.filterfalse(Path.exists, files)
    for path in _not_exist: log.error(S('file not found: {}', path))
    log.info(S('done'))


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
    paths1, paths2 = it.tee(map(Path, files_list))
    log_nonexistent(paths1)
    files_not_found = list(it.filterfalse(Path.exists, paths2))
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


def paths_from_flist(path):
    """
    Read Paths from file listing paths.

    Returns
    -------
    Map generator
        Path representation of the entries in the file.
    """
    with open(path, 'r') as fin:
        lines = map(str.strip, fin.readlines())
        valid = filter(bool, lines)
        no_comments = it.filterfalse(lambda x: x.startswith('#'), valid)

    return map(Path, no_comments)


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



def extract_from_tar(pdbtar, output='.'):
    """
    """
    tar = tarfile.open(pdbtar)
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

def make_destination_folder(dest):
    try:
        dest_ = Path(dest)
    except TypeError:
        return Path.cwd()
    else:
        dest_.mkdir(parents=True, exist_ok=True)
        return dest_

def read_PDBID_from_tar(tfile):
    """
    """



def read_PDBID_from_folder(folder):
    """
    Function description here.

    Parameters
    ----------

    Returns
    -------
    """

    pdblist = PDBList(glob_folder(folder, '*.pdb'))
    log.info(T('reading destination folder'))
    log.info(S(f'from: {folder}'))
    log.info(S(f'{str(pdblist)}'))
    log.info(S('done\n'))

    return pdblist


def save_dictionary(mydict, output='mydict.pickle'):
    options = {
        '.pickle': save_pickle_dict,
        '.json': save_json_dict,
        }

    options[Path(output).suffix](mydict, output)


def save_json_dict(mydict, output='mydict.json'):
    """
    """
    with open(output, 'w') as fout:
        json.dump(mydict, fout, indent=True, sort_keys=True)


def save_pickle_dict(mydict, output='mydict.pickle'):
    assert Path(output).suffix == '.pickle'
    with open(output, 'wb') as handle:
        pickle.dump(mydict, handle, protocol=pickle.HIGHEST_PROTOCOL)
