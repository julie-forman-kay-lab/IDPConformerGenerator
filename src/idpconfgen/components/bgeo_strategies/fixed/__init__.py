"""Module to build backbone geometry angles with fixed values."""
from itertools import cycle

from idpconfgen.core.build_definitions import (
    build_bend_angles_CA_C_Np1,
    build_bend_angles_Cm1_N_CA,
    build_bend_angles_N_CA_C,
    )


name = "fixed"


def get_cycle_bend_angles():
    """
    Return an infinite iterator of the bend angles.
    """
    return cycle((
        build_bend_angles_Cm1_N_CA,  # used for OMEGA
        build_bend_angles_N_CA_C,  # used for PHI
        build_bend_angles_CA_C_Np1,  # used for PSI
        ))
