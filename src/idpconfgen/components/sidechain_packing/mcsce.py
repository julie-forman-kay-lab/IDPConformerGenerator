"""
Implement MCSCE sidechain packing algorithm logic.

MCSCE repository at: https://github.com/THGLab/MCSCE
"""
from functools import partial

import numpy as np


mcsce_defaults = {
    'efunc_terms': ('lj', 'clash'),
    'mode': 'simple',
    'n_trials': 10,
    'batch_size': 16,
    'temperature': 300,
    'parallel_worker': 1,
    }

need_H_mask = True
only_H_mask = None


def add_mcsce_subparser(ap):
    """Add MCSCE related parameters to client."""
    group = ap.add_argument_group(
        title="MCSCE related parameters",
        description=(
            "Parameters configuring MCSCE sidechain sampling. "
            "Used only if `-scm mcsce` is selected."
            ),
        )
    group.add_argument(
        '--mcsce-mode',
        help=f'Sidechain sampling method. Defaults to {mcsce_defaults["mode"]!r}.',
        choices=('simple', 'exhaustive'),
        default=mcsce_defaults['mode'],
        )
    group.add_argument(
        '--mcsce-n_trials',
        help=f'Sampling trials in exhaustive mode. Defaults to {mcsce_defaults["n_trials"]}.',
        type=int,
        default=mcsce_defaults['n_trials'],
        )
    group.add_argument(
        '--mcsce-batch_size',
        help=f'The MCSCE batch size. Defaults to {mcsce_defaults["batch_size"]}.',
        type=int,
        default=mcsce_defaults['batch_size'],
        )
    group.add_argument(
        '--mcsce-temperature',
        help=f'Sampling temperature. Defaults to {mcsce_defaults["temperature"]}.',
        type=int,
        default=mcsce_defaults['temperature'],
        )
    return


def init_mcsce_sidechains(
        input_seq,
        template_masks,
        all_atom_masks,
        user_parameters=None,
        **ignore,
        ):
    """."""
    from mcsce.libs.libstructure import Structure
    from mcsce.libs.libenergy import prepare_energy_function
    from mcsce.core import build_definitions
    from mcsce.core.side_chain_builder import create_side_chain, initialize_func_calc

    _up = user_parameters or {}
    params = {**mcsce_defaults, **_up}

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
            return None, None

        if need_H_mask:
            # the atoms in final_structure are not ordered the same as the
            # atoms in cli_build. Instead of reordering the atoms in
            # final_structure, we extract a mask to take all sidechain atoms
            # and pass those to cli_build. The globals are used to perform
            # this operation only once.
            need_H_mask = False
            only_H_mask = np.where(
                np.in1d(
                    final_structure.data_array[:, 2],
                    ('N', 'CA', 'C', 'O', 'H1', 'H2', 'H3', 'OXT', 'H'),
                    invert=True,
                    ),
                )

        return all_atom_masks.all_sidechains, final_structure.coords[only_H_mask]  # noqa: E501

    return calc
