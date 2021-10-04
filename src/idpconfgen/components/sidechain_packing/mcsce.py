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
    from mcsce.libs.libstructure import Stucture
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
        raise ValueError("Mode has to be either simple or exhaustive.") from err

    efunc_partial = partial(
        prepare_energy_function,
        forcefield=ff_obj,
        terms=efunc_terms,
        )

    def calc(coords):

        s.coords = coords

        final_structure = create_side_chain(
            s,
            n_trials,
            efunc_partial,
            return_first_valid=return_first_valid,
            **params)

        if final_structure is None:
            return None
        else:
            return final_structure.coords

    return calc
