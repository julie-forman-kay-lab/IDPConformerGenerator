import numpy as np
from scipy.spatial import distance
from numpy import array

class ClashValidator:

    def __init__(self):
        self.clashes = 0

    def clash_found(self, current_chain, last_chain):
        last_chain_vector = np.dstack([np.array(list(atom.values())) for atom in last_chain][:-2]) #ORIGINAL
        # last_chain_vector = np.dstack([np.array(list(atom.values())) for atom in last_chain])
        l_axis1, l_axis2, l_axis3 = last_chain_vector.shape
        # stacking atoms on top of each other by reshaping it to (nx4, 3)
        last_chain_vector = last_chain_vector.reshape(l_axis3*l_axis1, l_axis2)

        current_chain_vector = np.dstack([np.array(list(atom.values())) for atom in current_chain][len(last_chain):-1]) #ORIGINAL
        # current_chain_vector = np.dstack([np.array(list(atom.values())) for atom in current_chain])
        c_axis1, c_axis2, c_axis3 = current_chain_vector.shape
        # stacking atoms on top of each other by reshaping it to (nx4, 3)
        current_chain_vector = current_chain_vector.reshape(c_axis3*c_axis1, c_axis2)

        distances = distance.cdist(last_chain_vector, current_chain_vector, 'euclidean')

        # VDW = {"N":1.55, "CA":1.7, "C":1.7, "O":1.52}
        # VDW = np.array([1.55, 1.7, 1.7, 1.52])
        VDW = np.array([1.7, 1.7, 1.7, 1.7])
        allowed_distances = np.add.outer(VDW,VDW)
        x, y = distances.shape #will be multiples of 4
        allowed_distances = np.tile(allowed_distances, (int(x/4), int(y/4)))
        
        clashes = (distances < allowed_distances)
        for row in clashes:
            if any(row):
                return True

        return False

        # compare the new loop with everything we had so far and make sure
        # no clashes exist
        # for seen_index in range(0, len(last_chain)-2):
            
        #     previous_coords = np.array(list(last_chain[seen_index].values()))

        #     for new_index in range(len(last_chain)+2, len(current_chain)):
        #         new_coords = np.array(list(current_chain[new_index].values()))

        #         try:
        #             distances = np.linalg.norm(new_coords-previous_coords, axis=1)
        #             clashes = any(distances < VDW)
        #         except ValueError:
        #             # we are at the last new_index, therefore its made up of only 3 atoms
        #             try:
        #                 distances = np.linalg.norm(new_coords-previous_coords[:-1], axis=1)
        #             except ValueError:
        #                 # COO
        #                 continue
        #             clashes = any(distances < VDW[:-1])
        #         if clashes:
        #             return True

        # return False