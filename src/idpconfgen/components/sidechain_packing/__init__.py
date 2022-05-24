"""
Modules to compute sidechains.

This package contains all methods integrated in idpconfgen to compute
side chains. It implements a strategy pattern where depending on the
sidechain option (str) the building workflow will use one function or
another.

For this reason, all functions computing sidechains should output the
same types and have a compatible input api.

All functions should return an array mask to be used in
`all_atom_coords` and return the coordinate compatible with such a mask.

For that, in general, the ``init_`` functions receive the template and
all_atom masks named tuples (see libbuild).
"""
from idpconfgen.components.sidechain_packing.faspr import init_faspr_sidechains
from idpconfgen.components.sidechain_packing.mcsce import (
    add_mcsce_subparser,
    init_mcsce_sidechains,
    )

sidechain_packing_methods = {
    'faspr': init_faspr_sidechains,
    'mcsce': init_mcsce_sidechains,
    }
"""Sidechain packing algorithms."""
_spm = tuple(sidechain_packing_methods.keys())


DEFAULT_SDM = 'faspr'


def add_sidechain_method(parser):
    """Add option to choose from sidechain packing algorithms."""
    parser.add_argument(
        '-scm',
        '--sidechainmethod',
        dest='sidechain_method',
        default=DEFAULT_SDM,
        choices=_spm,
        )

def get_sidechain_packing_parameters(parameters, mode):
    """Extract sidechain packaging related parameters."""
    assert mode in _spm
    spm_params = {}
    for param, value in parameters.items():
        if param.startswith(mode):
            parts = param.split("_")
            # we want to discard the first part of the parameter name
            # which corresponds to the `mode`
            spm_params['_'.join(parts[1:])] = value
    return spm_params
