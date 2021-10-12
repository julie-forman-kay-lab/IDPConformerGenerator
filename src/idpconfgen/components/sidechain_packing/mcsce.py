"""
Implement MCSCE sidechain packing algorithm logic.

MCSCE repository at: https://github.com/THGLab/MCSCE
"""
from functools import partial


mcsce_defaults = {
    'efunc_terms': ('lj', 'clash'),
    'n_trials': 200,
    'mode': 'simple',
    'temperature': 300,
    'parallel_worker': 1,
    }


def init_mcsce_sidechains(input_seq, **kwargs):
    """."""
    from mcsce.libs.libstructure import Structure
    from mcsce.libs.libenergy import prepare_energy_function
    from mcsce.core import build_definitions
    from mcsce.core.side_chain_builder import create_side_chain

    params = {**mcsce_defaults, **kwargs}

    # initiates only the backbone atoms
    s = Structure(fasta=input_seq)
    s.build()

    ff = build_definitions.forcefields["Amberff14SB"]
    ff_obj = ff(add_OXT=True, add_Nterminal_H=True)

    _mode = params.pop('mode')
    mcsce_sampling_options = {'simple': True, 'exhaustive': False}

    try:
        return_first_valid = mcsce_sampling_options[_mode]
    except KeyError as err:
        _msg = "Mode has to be either 'simple' or 'exhaustive'"
        raise ValueError(_msg) from err

    efunc_partial = partial(
        prepare_energy_function,
        forcefield=ff_obj,
        terms=params.pop('efunc_terms'),
        )

    params['efunc_creator'] = efunc_partial
    params['return_first_valid'] = return_first_valid

    def calc(coords):

        s.coords = coords

        final_structure = create_side_chain(s, **params)

        if final_structure is None:
            return None
        else:
            return final_structure.coords

    return calc
