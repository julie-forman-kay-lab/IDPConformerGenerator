"""Components controlling the energy threshold evaluation."""
from libfuncpy import vartial

from idpconfgen.libs.libcli import add_general_arg
from idpconfgen.libs.libenergyij import post_calc_options


et_type_help = \
"""How to calculate the energy for 'ij' atom pairs. `whole` calculates the
energy as a single float for the whole conformer (performs sum of pair
contributions. `pairs` calculates energy for individual pairs. In the current
implementation, the energy threshold (-etbb and -etss) will be compared to the
`whole` calculation or to the individual `pairs`. In `pairs`, if a single pair
has energy above threshold the conformer is discarded. If you use `pairs`, the
minimum values for `-etbb` and `-etss` should not be bellow 2 or 3 because pairs
at 1-4 (three bonds apart) will always have a small positive energy
contribution.  Hence, setting `-etbb` or `-etss` to 0 with `pairs` will result
in infinite calculations. If you choose `whole`, you can even define negative
numbers for `-etbb` and `-etss`."""

et_type_args = ['-et', '--et-type']
et_type_kwargs = {
    'help': et_type_help,
    'choices': post_calc_options,
    'dest': 'energy_type_ij',
    }

add_et_type_arg = vartial(add_general_arg, *et_type_args, **et_type_kwargs)
