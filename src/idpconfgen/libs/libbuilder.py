"""Library for the Conformer builder."""
import numpy as np

from idpconfgen.core import definitions as DEFS


class ConformerTemplate:
    """
    A conformer template for building.

    Used by the :class:`ConformerBuilder` as a molecule coordinate
        system.
    """
    def __init__(self, seq):
        
        self._seq = self._parse_seq(seq)
        
        # 1 is for the additional oxygen atom in the C-term carboxyl
        num_of_atoms = len(self.seq) * DEFS.num_bb_atoms + 1

        coords = np.full(
            3 * num_of_atoms,
            np.nan,
            )

        self._coords = coords.reshape(num_of_atoms, 3)
        
        self._atomnames = \
            np.array(
                DEFS.backbone_atoms * len(self.seq) + (DEFS.COO_atom,),
                dtype='<U1',
                )

    @property
    def seq(self):
        """
        The aminoacid sequence that represents the conformer.

        Is a tuple of aminoacids represented by 3-letter codes.
        """
        return self._seq

    @property
    def coords(self):
        """
        The XYZ coordinates that define the conformer.

        Is an Numpy array of length N * 4 + 1, where N is the length
        of the :attr:`seq` of the conformer, 4 is the number of backbone
        atoms (N, CA, C, O) and 1 accounts for the terminal oxygen
        of the carboxyl group.
        """
        return self._coords
    
    @property
    def atomnames(self):
        """
        An array of the backbone atom names.

        Does not include the terminal oxygen atom from the carboxyl group,
        therefore:

        :attr:`atomnames`.size == :attr:`coords`.shape[0] - 1
        """
        return self._atomnames

    def get_coord(self, residue_pos, atomname):
        """
        Return the XYZ coordinates for atom in residue position.

        Parameters
        ----------
        residue_pos : int
            The residue position with the :attr:`seq`.
            Position is 0-indexed.

        atomname : str
            The atom name. For example ``N``, ``C``.
        """
        return

    def add_atom_coords(self, residue_index, atom_name, coords):
        """
        Parameters
        ----------
        residue_index : int
            The residue 0-index according to :attr:`seq` length.

        atom_name : str
            The atom name. For example ``N``, ``CA``.

        coord : numpy.array of shape (3,) dtype=np.float32
            A Numpy array with the XYZ space coordinates of the atom.
        """

        # -1 because we are 0-indexed
        index = residue_index * DEFS.num_bb_atoms
        try:
            index += DEFS.backbone_atoms.index(atom_name)
        except ValueError:  # we are adding carboxyl oxygen
            index = -1

        self.coords[index,:] = coords
    
    def is_complete(self):
        return not np.any(np.isnan(self.coords))

    @staticmethod
    def _parse_seq(seq):
        """
        seq (list or str) -> tuple of 3 char strings
        """
        try:
            sequence = tuple(DEFS.aa1to3[c] for c in seq)
        except KeyError:
            # the sequece is not a sequence of 1 letter codes
            if all(c in DEFS.aa3to1.keys() for c in seq):
                sequence = tuple(seq)
            else:
                raise ValueError('seq aminoacids codes are not valid')
        return sequence
