import glob

from idpconfgen import log, Path


def glob_folder(folder, ext):
    """
    Lists files with extention `ext` in `folder`.

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
    list
        SORTED list of matching results
    """

    ext = f"*.{ext.strip().lstrip('*').lstrip('.')}"
    files = sorted(glob.glob(Path(folder, ext).str()))
    log.debug(f'folder {folder} read {len(files)} files with extension {ext}')
    return files
