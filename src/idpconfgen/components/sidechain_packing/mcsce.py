"""
Implement MCSCE sidechain packing algorithm logic.

MCSCE repository at: https://github.com/THGLab/MCSCE
"""
from functools import partial
from idpconfgen.core.build_definitions import convert_one2three


mcsce_defaults = {
    'efunc_terms': ('lj', 'clash'),
    'n_trials': 200,
    'batch_size': 16,
    'mode': 'simple',
    'temperature': 300,
    'parallel_worker': 1,
    }


def init_mcsce_sidechains(input_seq, template_masks, all_atom_masks, all_atom_labels, **kwargs):
    """."""
    from mcsce.libs.libstructure import Structure
    from mcsce.libs.libenergy import prepare_energy_function
    from mcsce.core import build_definitions
    from mcsce.core.side_chain_builder import create_side_chain, initialize_func_calc

    params = {**mcsce_defaults, **kwargs}

    # initiates only the backbone atoms
    s = Structure(fasta=input_seq)
    s.build()
    s = s.remove_side_chains()

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

    # params['efunc_creator'] = efunc_partial
    input_seq_3 = convert_one2three(input_seq)
    initialize_func_calc(efunc_partial, input_seq_3, s)
    params['return_first_valid'] = return_first_valid

    def calc(coords):

        s.coords = coords[template_masks.non_sidechains]

        final_structure = create_side_chain(s, **params)

        final_structure.reorder_by_atom_labels(all_atom_labels)

        if final_structure is None:
            return None
        else:
            return True, final_structure.coords

    return calc
