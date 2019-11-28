import numpy as np
from scipy.spatial import distance
from numpy import array
import math

class ClashValidator:

    def __init__(self):
        self.clashes = 0

    def clash_found_vectorized(self, new_chain, last_chain):
        """
        This method calculates every possible distance between 
        each atom within new_chain and last_chain using
        numpy.
        """

        # concatenate every residue containing four atoms into one numpy array
        last_chain_vector = np.array([atoms_array for residue in last_chain[:-2] for atoms_array in residue.values()])
        new_chain_vector = np.array([atoms_array for residue in new_chain for atoms_array in residue.values()])

        # get the euclidean distance between each atom
        distances = distance.cdist(last_chain_vector, new_chain_vector, 'euclidean')

        # VDW = {"N":1.55, "CA":1.7, "C":1.7, "O":1.52}
        VDW = np.array([1.55, 1.7, 1.7, 1.52])
        allowed_distances = np.add.outer(VDW,VDW)
        x, y = distances.shape
        #the extra slicing is to allow us to include the last 3 atoms
        allowed_distances = np.tile(allowed_distances, (math.ceil(x/4), math.ceil(y/4)))[:,:y] 

        clashes = (distances < allowed_distances)
        if np.any(clashes):
                self.clashes += 1
                return True
        return False

    def clash_found(self, new_chain, last_chain):
        """
        This method is the same as clash_found_vectorized but manually
        does the work. Ie iterates once over last_chain and for every
        value, it iterates once over new_chain and calculates
        the necessary distances.

        """
        VDW = np.array([1.55, 1.7, 1.7, 1.52])
        allowed_distances = np.add.outer(VDW,VDW)

        # compare the new loop with everything we had so far and make sure
        # no clashes exist
        for seen_index in range(0, len(last_chain)-2):
            
            previous_coords = np.array(list(last_chain[seen_index].values()))

            for new_index in range(0, len(new_chain)):
                new_coords = np.array(list(new_chain[new_index].values()))

                try:
                    distances = distance.cdist(previous_coords, new_coords, 'euclidean')
                    clashes = distances < allowed_distances
                    if np.any(clashes):
                        self.clashes += 1
                        return True
                except ValueError:
                        # we are at the last new_index, therefore its made up of only 3 atoms
                        distances = distance.cdist(previous_coords, new_coords, 'euclidean')
                        clashes = distances < allowed_distances[:-1].T
                        if np.any(clashes):
                            self.clashes += 1
                            return True
        return False