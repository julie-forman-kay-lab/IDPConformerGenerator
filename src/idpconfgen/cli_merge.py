"""
Collates and re-numbers conformers within subfolders of a given folder.
Has the option to remove previous and rename current conformer files after processing.

For ex., if there are 5 subfolders filled with 100 conformers each,
copy and re-number them to make 500 conformers in the working dir. or
destination dir.

USAGE:
    $ idpconfgen merge --target <PATH_TO_SUBFOLDERS>
    $ idpconfgen merge \\
        --target <PATH_TO_SUBFOLDERS> \\
        --destination <PATH_TO_OUTPUT> \\
        --prefix <CUSTOM_NAME> \\
        --delete \\
"""

import argparse, shutil, glob

from idpconfgen.libs.libio import make_folder_or_cwd
from idpconfgen.libs import libcli
from idpconfgen import log
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = 'idpconfgen_merge'
_name = 'merge'
_help = 'Collate and rename/renumber conformers from subfolders into current/new folder.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-tgt',
    '--target',
    help=('Directory containing subfolders of conformers to merge. Required.'
          ),
    type=str,
    required = True,
    default='./',
)

ap.add_argument(
    '-des',
    '--destination',
    help=('Directory to collate conformers. Optional.'
          ),
    type=str,
    default=None,
)

ap.add_argument(
    '-pre',
    '--prefix',
    help=('Name to replace default `conformer` for each file. Optional.'
          ),
    type=str,
    default='conformer',
)

ap.add_argument(
    '-del',
    '--delete',
    help='Deleted conformers and subfolders after merging. Defaults to True.',
    action='store_true',
)

def main(
    target,
    destination=None,
    prefix='conformer',
    delete=False,
    **kwargs,
    ):
    """
    Perform main logic of the script.

    Parameters
    ----------
        target : str, required
            Path to the directory with all of the subfolders.
        
        destination : str, optional
            Path to the directory to save the collated conformers from the target.
        
        prefix : str, optional
            Replace the default `conformer` name for each file.
        
        delete : bool, optional
            Flag to turn on deletion of conformers and subfolders of interest.
            Defaults to True.
    """
    init_files(log, LOGFILESNAME)
    
    # initialize some variables
    count=0
    if destination == None: destination = target
    destination = str(make_folder_or_cwd(destination))
    
    log.info(T('Attempting to collate files'))
    log.info(S('Obtaining paths to all conformers...'))
    # get all the paths to the subfolders
    pdbpaths=[]
    subpaths=[]
    for spath in glob.glob(f'{target}/*/'): subpaths.append(spath)
    for fpath in glob.glob(f'{target}/*/*.pdb', recursive=True): pdbpaths.append(fpath)
    log.info(S('done'))
    
    log.info(S(f'Copying {len(pdbpaths)} files with prefix {prefix}...'))
    for path in pdbpaths:
        # using copy2 to preserve metadata and file permissions
        try:
            shutil.copy2(path, destination+f'/{prefix}_{count}.pdb')
        except OSError as e:
            log.info(S(f'Error: {path} : {e.strerror}'))
        count+=1
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
