"""Store functions for using folded domains as templates."""
from numba import njit
#def calc_energy(conf_coords):
#    for coord in conf_coords:
#        dist = np.sum(template_coords - coord)**2, axis=1)
#        if np.any(dist < 12.96):
#            return True
#    return False


def calc_dist_threshold_against_template(folded_coords):
    """
    Assert no atom from target is closer than a threshold to the template.
    """
    @njit
    def calc_dist_threshold(conf_coords):
        for coord in conf_coords:
            for tcoord in folded_coords:
                x = coord[0] - tcoord[0]
                y = coord[1] - tcoord[1]
                z = coord[2] - tcoord[2]
                dist = (x * x + y * y + z * z)
                if dist < 12.96:
                    return True
        return False
    return calc_dist_threshold
