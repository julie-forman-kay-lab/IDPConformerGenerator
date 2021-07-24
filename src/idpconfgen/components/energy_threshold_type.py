"""Components controlling the energy threshold evaluation."""
from libfuncpy import vartial

from idpconfgen.libs.libcli import add_general_arg
from idpconfgen.libs.libenergyij import post_calc_options


et_type_help = \
"""How to evaluate the energy threshold. `whole` evaluates threshold against the
sum of all contributions. `pairs` evaluates threshold against the pair of
highest energy. If you use `pairs` define a minimum -etbb and -etss of 2 or 3,
depending on the system. This is because pairs at 1-4 (three bonds apart) will
always have a small positive energy contribution. If you choose `whole`, you can
even defined negative numbers for -etbb and -etss.
"""

et_type_args = ['-et', '--et-type']
et_type_kwargs = {
    'help': et_type_help,
    'choices': post_calc_options,
    'dest': 'energy_type_ij',
    }

add_et_type_arg = vartial(add_general_arg, *et_type_args, **et_type_kwargs)
