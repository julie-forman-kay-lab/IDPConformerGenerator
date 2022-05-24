"""Bond geomtry strategies."""
from idpconfgen.components.bgeo_strategies.int2cart.bgeo_int2cart import \
    name as bgeo_int2cart_name
    
from idpconfgen.components.bgeo_strategies.sampling import \
    name as bgeo_sampling_name


bgeo_strategies_default = bgeo_sampling_name
"""The default bond geometry sampling strategy."""


bgeo_strategies = (
    bgeo_sampling_name,
    bgeo_int2cart_name,
    )
"""Available bond geometry sampling strategies."""
