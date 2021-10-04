"""Modules to compute sidechains."""
from idpconfgen.components.sidechain_packing.faspr import init_faspr_sidechains
from idpconfgen.components.sidechain_packing.mcsce import init_mcsce_sidechains


sidechain_packing_methods = {
    'faspr': init_faspr_sidechains,
    'mcsce': init_mcsce_sidechains,
    }


def add_sidechain_method(parser):
    """Add option to choose from sidechain packing algorithms."""
    parser.add_argument(
        '-scp',
        '--sidechainmethod',
        dest='sidechain_method',
        default='faspr',
        choices=list(sidechain_packing_methods.keys()),
        )
