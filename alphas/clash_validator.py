import numpy as np
from numpy import array

class ClashValidator:

    def __init__(self):
        self.clashes = 0

    def clash_found(self, current_conformer, last_conformer):
        # VDW = {"CA":1.7,"C":1.7,"N":1.55,"O":1.52}
        VDW = np.array([3.2, 3.2, 3.1, 3.04])
        
        current_chain = current_conformer.coords
        last_chain = last_conformer.coords
        # compare the new loop with everything we had so far and make sure
        # no clashes exist
        for seen_index in range(0, len(last_chain)-1):
            
            previous_coords = np.array(list(last_chain[seen_index].values()))

            for new_index in range(len(last_chain)+2, len(current_chain)):
                new_coords = np.array(list(current_chain[new_index].values()))

                try:
                    distances = np.linalg.norm(new_coords-previous_coords, axis=1)
                    clashes = any(distances < VDW)
                except ValueError:
                    # we are at the last new_index, therefore its made up of only 3 atoms
                    try:
                        distances = np.linalg.norm(new_coords-previous_coords[:-1], axis=1)
                    except ValueError:
                        # COO
                        continue
                    clashes = any(distances < VDW[:-1])
                if clashes:
                    return True

        return False