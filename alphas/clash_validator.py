import numpy as np
from scipy.spatial import distance
from scipy.spatial import KDTree
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

    def clash_found_kdtree(self, new_chain,last_chain ):
        """
        Abstract and pretend all molecules are one. This will make things a lot faster
        """
        # concatenate every residue containing four atoms into one numpy array
        last_chain_vector = np.array([atoms_array for residue in last_chain[:-2] for atoms_array in residue.values()])
        new_chain_vector = np.array([atoms_array for residue in new_chain for atoms_array in residue.values()])

        last_chain_kdtree = KDTree(last_chain_vector)
        new_chain_kdtree = KDTree(new_chain_vector)

        # VDW = {"N":1.55, "CA":1.7, "C":1.7, "O":1.52}
        VDW = np.array([1.55, 1.7, 1.7, 1.52])
        allowed_distances = np.add.outer(VDW,VDW)

        # VDW = {"N":1.55, "CA":1.7, "C":1.7, "O":1.52}
        clashes = last_chain_kdtree.sparse_distance_matrix(new_chain_kdtree, 3.4) #maximum value to consider
        for indices, distance in zip(clashes.keys(),clashes.values()):
            if distance <= 3.04:
                # we don't care what the atoms are, this is a clash
                return True

            row = indices[0] % 4
            column = indices[1] % 4

            if distance <= allowed_distances[row, column]:
                return True
        
        return False


        print(clashes)
