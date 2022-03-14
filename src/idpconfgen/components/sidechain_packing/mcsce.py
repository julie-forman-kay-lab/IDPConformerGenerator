"""
Implement MCSCE sidechain packing algorithm logic.

MCSCE repository at: https://github.com/THGLab/MCSCE
"""
from functools import partial

import numpy as np


mcsce_defaults = {
    'efunc_terms': ('lj', 'clash'),
    'n_trials': 200,
    'batch_size': 16,
    'mode': 'simple',
    'temperature': 300,
    'parallel_worker': 1,
    }

need_H_mask = True
only_H_mask = None


def init_mcsce_sidechains(input_seq, template_masks, all_atom_masks, **kwargs):
    """."""
    from mcsce.libs.libstructure import Structure
    from mcsce.libs.libenergy import prepare_energy_function
    from mcsce.core import build_definitions
    from mcsce.core.side_chain_builder import create_side_chain, initialize_func_calc

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

    initialize_func_calc(
        partial(
            prepare_energy_function,
            batch_size=params.pop('batch_size'),
            forcefield=ff_obj,
            terms=params.pop("efunc_terms"),
            ),
        structure=s,
        )

    params['return_first_valid'] = return_first_valid

    def calc(coords):

        # yes, I could have made a class, store these attributes, bla bla
        # ... but I didn't want to :-)
        global need_H_mask
        global only_H_mask

        s.coords = coords[template_masks.non_sidechains]

        final_structure = create_side_chain(s, **params)

        if final_structure is None:
            return None

        if need_H_mask:
            need_H_mask = False
            only_H_mask = np.where(
                np.in1d(
                    final_structure.data_array[:, 2],
                    ('N', 'CA', 'C', 'O', 'H1', 'H2', 'H3', 'OXT', 'H'),
                    invert=True,
                    ),
                )

        return all_atom_masks.all_sidechains, final_structure.coords[only_H_mask]

    return calc
