"""
Module contains variable definitions and
miscellaneous utility functions to help the `eisd` subclient.

Functions and logic inspired/imported from https://github.com/THGLab/X-EISD/
"""
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


# The following two functions have been imported from:
# https://github.com/THGLab/X-EISD/blob/master/eisd/utils/miscell.py
def make_pairs():
    all = ['saxs', 'cs', 'fret', 'jc', 'noe', 'pre', 'rdc', 'rh']
    pairs=[]
    for i in range(len(all)):
        for j in range(i+1, len(all)):
            pairs.append([all[i],all[j]])
    return pairs


def modes(mode):
    flags = {'saxs': False, 'cs': False, 'fret':False, 'jc': False,
             'noe': False,  'pre': False, 'rdc':False, 'rh': False
             }

    if mode is 'all':
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