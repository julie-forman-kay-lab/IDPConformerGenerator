import glob

from idpconfgen import log, Path

def glob_folder(folder, ext):
    ext = f"*.{ext.lstrip('*').lstrip('.')}"
    files = sorted(glob.glob(Path(folder, ext).str()))
    log.debug(f'folder {folder} read {len(files)} files with extension {ext}')
    return files
