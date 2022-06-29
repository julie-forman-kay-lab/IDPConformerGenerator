"""# noqa: D205, D210, D400, E501
Interface for X-EISD logic to score and reweight conformer ensembles.

Has to ability to store and optimize conformer ensembles towards
back-calculated and experimental data.

Inspired and imported from: https://github.com/THGLab/X-EISD/

NOTES
-----
For `--filepath`, please name all your experimental files
ending with the extension(s) corresponding to the experiment
of interest and starting with `exp_`.
    
    For example:
    exp_*.saxs, exp_*.cs, exp_*.fret, exp_*.jc,
    exp_*.noe, exp_*.pre, exp_*.rdc, exp_*.rh
    
For back-calculated files, please start them with the `back_`
prefix, ending with the same file extensions as seen above.
    
    For example:
    back_*.saxs, back_*.cs, back_*.fret, etc.
    
Please note that for every `exp_*` file, there must be a corresponding
`back_*` file. EXCEPT for jc, noe, pre as those back-calculations are
handled internally.

USAGE:
    $ idpconfgen eisd [PATH-TO-PDBS] --filepath
"""
import argparse
import os
import shutil
from functools import partial
from turtle import back

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_path_bundle,
    FileReaderIterator,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libparse import pop_difference_with_log
from idpconfgen.logger import S, T, init_files, report_on_crash

from idpconfgen.components.eisd import (
    eisd_modes,
    eisd_modules,
    eisd_optimization_types,
    default_type,
    default_mode,
    )

LOGFILESNAME = '.idpconfgen_eisd'
TMPDIR = '__tmpeisd__'

_name = 'eisd'
_help = 'Score and reweight generated ensembles.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-fpath',
    '--filepath',
    help=('Path to the folder containing experimental '
          'and back-calculation data files.'),
    type=str,
    required=True,
    default=None,
    )

ap.add_argument(
    '-eps',
    '--epochs',
    help='Number of times to run main optimization, defaults to 250.',
    type=int,
    default=250,
    required=False,
    )

ap.add_argument(
    '-opt',
    '--optimization',
    help='If optimization is required, defaults to True.',
    action='store_true',
    required=False,
    )

ap.add_argument(
    '-m',
    '--mode',
    help='Optimization mode on all, single, or dual experimental observables.',
    choices=eisd_modes,
    default=default_mode,
    required=False,
    )

ap.add_argument(
    '--optimization-type',
    help='Type of optimzation done. Defaults to max.',
    choices=eisd_optimization_types,
    default=default_type,
    required=False,
    )

ap.add_argument(
    '--beta',
    help=('Hyperparameter for "mc" optimization type'
          ' (Metropolis Monte Carlo). '
          'Defaults to 0.1.'         
          ),
    default=0.1,
    type=float,
    required=False,
    )

libcli.add_argument_output(ap)

ap.add_argument(
    '--tmpdir',
    help=(
    'Temporary directory to store data during calculation '
    'if needed.'
    ),
    type=Path,
    default=TMPDIR,
    )

libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        filepath=None,
        optimization=False,
        optimization_type=default_type,
        mode=default_mode,
        epochs=250,
        beta=0.1,
        output="./eisd_output",
        tmpdir=TMPDIR,
        ncores=1,
        func=None,
        ):
    """
    Main logic of processing, scoring, and optimizing ensembles.
    
    Parameters
    ----------
    pdb_files : str or Path, required
        Path to PDB files on the disk. Accepts TAR file.
    
    filepath : str or Path, required
        Path to the folder containing experimental and
        back-calculated data files.
    
    optimization : Bool, optional
        Flag to turn on or off process of optimizing an ensemble using
        experimental and back-calculated data.
        Defaults to False.
    
    optimization_type : str, optional
        Type of optimization to run, requires `optimization`
        flag to be true.
        Defaults to `max`.
    
    mode : str, optional
        EISD scoring mode.
        Defaults to `all`.
        
    epochs : int, optional
        Number of times to run main optimization.
        Defaults to 250.
        
    beta : float, optional
        Hyperparameter for "mc" type optimization.
        Defaults to 0.1.
    
    output : str or Path, optional
        Path to the folder to store any files or conformers.
        Defaults to working directory.
        
    ncores : int, optional
        Number of workers to use for multiprocessing.
        Defaults to 1.   
    """
    init_files(log, LOGFILESNAME)
    
    log.info(T('Reading conformer ensemble paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, FileNotFoundError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))
    
    log.info(T('Checking experimental data files'))
    
    exp_paths=[]
    back_paths=[]
    valid_exp_modules=[]
    valid_back_modules=[]
    all_files = [f for f in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, f))]  # noqa: E501
    for f in all_files:
        if f.startswith('exp_'):
            if f.endswith(eisd_modules):
                exp_paths.append(os.path.join(filepath, f))
                _ext = f[f.rindex('.')+1:]
                valid_exp_modules.append(f'.{_ext}')
        elif f.startswith('back_'):
            if f.endswith(eisd_modules):
                back_paths.append(os.path.join(filepath, f))
                _ext = f[f.rindex('.')+1:]
                valid_back_modules.append(f'.{_ext}')
    
    valid_exp_modules.sort()
    valid_back_modules.sort()
    
    if valid_exp_modules == []:
        log.info(S('WARNING: no valid experimental files found.'
                   ' Please refer to the help documentation for'
                   ' this module.'
                   ))
        return
    else:
        if valid_exp_modules != valid_back_modules:
            _diff = tuple(set(valid_exp_modules) ^ (set(valid_back_modules)))
            exp_paths = [exp for exp in exp_paths if not exp.endswith(_diff)]
            back_paths = [bck for bck in back_paths if not bck.endswith(_diff)]
            log.info(S('Note: found inconsistent experimental and back-calculation'
                       ' data pairs. Keeping only paths of matching pairs of data.'
                       ))
    
    log.info(S('done'))
    
    
    
    
    if _istarfile:
        shutil.rmtree(tmpdir)
    
    return



if __name__ == '__main__':
    libcli.maincli(ap, main)