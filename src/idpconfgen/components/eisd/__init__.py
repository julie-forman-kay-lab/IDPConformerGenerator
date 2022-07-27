"""
Module contains variable definitions and
miscellaneous utility functions to help the `eisd` subclient.

Functions and logic inspired/imported from https://github.com/THGLab/X-EISD/
"""
import os

exp_idx = 'index'
exp_resnum = 'resnum'
exp_val = 'value'
exp_max = 'upper'
exp_atmID = 'atomname'
exp_min = 'lower'
exp_err = 'error'

eisd_run_all = 'all'
eisd_run_single = 'single'
eisd_run_pairs = 'pairs'
eisd_modes = (
    eisd_run_all,
    eisd_run_single,
    eisd_run_pairs
    )
default_mode = eisd_run_all

opt_max = 'max'
opt_mc = 'mc'
eisd_optimization_types = (
    opt_max,
    opt_mc,
    None,
    )
default_type = opt_max

parse_mode_exp = 'exp'
parse_mode_back = 'bc'
parse_modes = (
    parse_mode_exp,
    parse_mode_back,
    )

saxs_name = 'saxs'
cs_name = 'cs'
fret_name = 'fret'
jc_name = 'jc'
noe_name = 'noe'
pre_name = 'pre'
rdc_name = 'rdc'
rh_name = 'rh'
eisd_modules = (
    saxs_name,
    cs_name,
    fret_name,
    jc_name,
    noe_name,
    pre_name,
    rdc_name,
    rh_name,
    )

# define back calculation uncertainties
# refer to Lincoff et al. 2020 for details
default_bc_errors = {
    pre_name: 0.0001,
    noe_name: 0.0001,
    saxs_name: 0.006,
    fret_name: 0.0074,
    rh_name: 0.812,
    rdc_name: 0.88,
    # cs error reported from UCBShift
    cs_name: {'C': 1.31, 'CA': 0.97, 'CB': 1.29, 'H': 0.38, 'HA': 0.29} 
    # J-coupling errors set by default
    }


# The following two functions have been imported from:
# https://github.com/THGLab/X-EISD/blob/master/eisd/utils/miscell.py
# https://github.com/Oufan75/X-EISD/blob/master/eisd/utils.py
def make_pairs(all):
    pairs=[]
    for i in range(len(all)):
        for j in range(i + 1, len(all)):
            pairs.append([all[i], all[j]])
    return pairs


def modes(mode, all):
    flags = {}
    for prop in all:
        flags[prop] = False

    if mode is eisd_run_all:
        return {flag:True for flag in flags}
    elif type(mode) is list:
        for flag in mode:
            flags[flag] = True
        return flags
    elif type(mode) is str:
        flags[mode] = True
        return flags

    else:
        raise ValueError("The mode in the main function is not recognized.")


def meta_data(fpath):
    """
    Function filters through experimental and back-calculated
    data files to return only the paths that exist in both cases.
    
    Automatically removes paths to datafiles if e.g. there are more
    back-calculated data than experimental.
    
    Parameters
    ----------
    fpath : str or Path
        Path to the directory containing both experimental
        and 

    Returns
    -------
    meta : dict
        First layer keys of `exp` and `bc` having values of dictionaries
        with keys being the module and the value being the path to the data.
    
    errlog : list
        List of errors to relay to the user if there are any.
    """
    exp_paths = []
    back_paths = []
    valid_exp_modules = []
    valid_back_modules = []
    
    meta = {}
    errlog = []
    
    all_files = [f for f in os.listdir(fpath) if os.path.isfile(os.path.join(fpath, f))]  # noqa: E501
    for f in all_files:
        if f.startswith('exp_'):
            if f.endswith(eisd_modules):
                exp_paths.append(os.path.join(fpath, f))
                _ext = f[f.rindex('.') + 1:]
                valid_exp_modules.append(f'.{_ext}')
        elif f.startswith('back_'):
            if f.endswith(eisd_modules):
                back_paths.append(os.path.join(fpath, f))
                _ext = f[f.rindex('.') + 1:]
                valid_back_modules.append(f'.{_ext}')
    
    valid_exp_modules.sort()
    valid_back_modules.sort()
    
    if valid_exp_modules == []:
        errlog.append(
            'WARNING: no valid experimental files found.'
            ' Please refer to the help documentation for'
            ' this module.'
            )
        return meta, errlog
    else:
        differences = tuple(set(valid_exp_modules) ^ (set(valid_back_modules)))
        
        if differences:
            exp_paths = [exp for exp in exp_paths if not exp.endswith(differences)]
            back_paths = [bck for bck in back_paths if not bck.endswith(differences)]
            errlog.append(
                'Note: found inconsistent experimental and back-calculation'
                ' data pairs. Keeping only paths of matching pairs of data.'
                )
            errlog.append(f'Excluding: {differences}...')
    
    EXP_DATA_FILENAMES = {}
    BACK_DATA_FILENAMES = {}
    
    for module in eisd_modules:
        for exp in exp_paths:
            if module in exp:
                EXP_DATA_FILENAMES[module] = exp
        for bc in back_paths:
            if module in bc:
                BACK_DATA_FILENAMES[module] = exp

    meta = {
        parse_mode_exp: EXP_DATA_FILENAMES,
        parse_mode_back: BACK_DATA_FILENAMES,
        }
    
    return meta, errlog
