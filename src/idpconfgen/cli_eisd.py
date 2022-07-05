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

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libfunc import consume
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_path_bundle,
    FileReaderIterator,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import consume_iterable_in_list, pool_function, starunpack
from idpconfgen.libs.libparse import pop_difference_with_log
from idpconfgen.logger import S, T, init_files, report_on_crash

from idpconfgen.components.eisd import (
    eisd_run_all,
    eisd_run_pairs,
    eisd_modes,
    eisd_modules,
    eisd_optimization_types,
    default_type,
    default_mode,
    parse_mode_exp,
    parse_mode_back,
    make_pairs,
    meta_data,
    modes,
    )
from idpconfgen.components.eisd.optimizer import core_eisd
from idpconfgen.components.eisd.parser import parse_data

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
        output="./",
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
        Path to the folder to store eisd outputs.
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
    filenames, errs = meta_data(filepath)
    if errs:
        for e in errs:
            log.info(S(e))
        if filenames == {}:
            log.info(S('done'))
            return
    log.info(S('done'))
    
    exp_paths = parse_data(filenames[parse_mode_exp], mode=parse_mode_exp)
    back_paths = parse_data(filenames[parse_mode_back], mode=parse_mode_back)
    
    _ens_size = len(pdbs2operate)
    _output = []
    if mode == eisd_run_all:
        _output.append(os.path.join(output, mode))
    elif mode == eisd_run_pairs:
        pairs = make_pairs()
        for pair in pairs:
            _output.append(os.path.join(output, "%s_%s_%s"%(mode, pair[0], pair[1])))
    else:
        for module in eisd_modules:
            _output.append(os.path.join(output, "%s_%s"%(mode, module)))
    
    execute = partial (
        report_on_crash,
        core_eisd,
        exp_data=exp_paths,
        bc_data=back_paths,
        ens_size=_ens_size,
        epochs=epochs,
        mode=mode,
        beta=beta,
        opt_type=optimization_type,
        )
    
    execute_pool = pool_function(execute, structure=pdbs2operate, ncores=ncores)
    
    
    
    if _istarfile:
        shutil.rmtree(tmpdir)
        
    log.info(S('done'))
    
    return



if __name__ == '__main__':
    libcli.maincli(ap, main)
