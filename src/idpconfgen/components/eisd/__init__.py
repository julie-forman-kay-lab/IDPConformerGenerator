"""
Module contains variable definitions and
miscellaneous utility functions to help the `eisd` subclient.

Functions and logic inspired/imported from https://github.com/THGLab/X-EISD/
"""
eisd_modes = (
    'all',
    'single',
    'pairs',
    )

eisd_optimization_types = (
    'max',
    'mc',
    None,
    )

eisd_modules = (
    'saxs',
    'cs',
    'fret',
    'jc',
    'noe',
    'pre',
    'rdc',
    'rh',
    )

default_mode = 'all'
default_type = 'max'


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