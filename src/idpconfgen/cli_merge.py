r"""# noqa: D205, D400, E501
Collates and re-numbers conformers within subfolders of a given folder.
Has the option to remove previous and rename current conformer files after processing.

For ex., if there are 5 subfolders filled with 100 conformers each,
copy and re-number them to make 500 conformers in the working dir. or
destination dir.

USAGE:
    $ idpconfgen merge --source <PATH_TO_SUBFOLDERS>
    $ idpconfgen merge \
        --source <PATH_TO_SUBFOLDERS> \
        --destination <PATH_TO_OUTPUT> \
        --prefix <CUSTOM_NAME> \
        --delete
"""
import argparse
import glob
import shutil
from pathlib import Path

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import make_folder_or_cwd
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = '.idpconfgen_merge'
_name = 'merge'
_help = ('Collate and rename/renumber conformers '
         'from subfolders into current/new folder.'
         )

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-src',
    '--source',
    help=(
        'Directory containing subfolders of conformers to merge. '
        'Defaults to the current working directory.'
        ),
    type=Path,
    required=True,
    default=None,
    )

ap.add_argument(
    '-des',
    '--destination',
    help=(
        'Directory to collate conformers. Optional. '
        'Defaults to the current working directory.'
        ),
    type=Path,
    default=None,
    )

ap.add_argument(
    '-pre',
    '--prefix',
    help='Common prefix of the new files. Defaults to `conformer`.',
    type=str,
    default='conformer',
    )

ap.add_argument(
    '-del',
    '--delete',
    help='Deleted conformers and subfolders after merging. Defaults to False.',
    action='store_true',
    )


def main(
        source,
        destination=None,
        prefix='conformer',
        delete=False,
        **kwargs,
        ):
    """
    Perform main logic of the script.

    Parameters
    ----------
    source : str, required
        Path to the directory with all of the subfolders.

    destination : str, optional
        Path to the directory to save the collated conformers from the source.

    prefix : str, optional
        Replace the default `conformer` name for each file.

    delete : bool, optional
        Flag to turn on deletion of conformers and subfolders of interest.
        Defaults to True.
    """
    init_files(log, LOGFILESNAME)
    
    # initialize some variables
    count = 0

    if destination is None:
        destination = source

    source = make_folder_or_cwd(source)
    destination = make_folder_or_cwd(destination)

    log.info(T('Attempting to collate files'))
    log.info(S('Obtaining paths to all conformers...'))
    # get all the paths to the subfolders
    subpaths = list(glob.glob(f'{source}/*/'))
    pdbpaths = list(glob.glob(f'{source}/*/*.pdb', recursive=True))

    log.info(S('done'))

    log.info(S(f'Copying {len(pdbpaths)} files with prefix {prefix}...'))
    for path in pdbpaths:
        # using copy2 to preserve metadata and file permissions
        try:
            shutil.copy2(path, Path(destination, f'{prefix}_{count}.pdb'))
        except OSError as e:
            log.info(S(f'Error: {path} : {e.strerror}'))
        count += 1
    log.info(S('done'))

    if delete:
        log.info(S('Cleaning up subfolders and PDB files...'))
        for path in subpaths:
            try:
                shutil.rmtree(path)
            except OSError as e:
                log.info(S(f'Error: {path} : {e.strerror}'))
        log.info(S('done'))

    return
