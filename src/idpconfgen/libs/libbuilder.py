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
                DEFS.backbone_atoms * len(self.seq) + (DEFS.COO_name,),
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
        return self.coords[self._get_index(residue_pos, atomname)]

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
        self.coords[self._get_index(residue_index, atom_name),:] = coords

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

    @staticmethod
    def _get_index(residue_index, atom_name):
        try:
            index = DEFS.backbone_atoms.index(atom_name)
        except ValueError:  # we are adding carboxyl oxygen, are we?
            if atom_name == DEFS.COO_name:
                return -1
            else:
                raise ValueError('atom_name not valid: {!r}'.format(atom_name))
        else:
            index += residue_index * DEFS.num_bb_atoms
            return index


class ConformerBuilderNeRF:
    """
    Conformer builder.
    """
    def __init__(
            self,
            conformer,
            angledb,
            frag_size=None,
            ):
        self._conformer = conformer
        self._angledb = angledb
        self.frag_size = frag_size
    
    @property
    def conformer(self):
        return self._conformer

    @property
    def angledb(self):
        return self._angledb

    def _make_coord(self, atomname):
        
        parent_coord = self.conformer.get_coord(
            residue_index + atom.poff,
            atom.name,
            )

        xaxis_coord = self.conformer.get_coord(
            residue_index + atom.xoff,
            atom.name,
            )

        yaxis_coord = self.conformer.get_coord(
            residue_index + atom.yoff,
            atom.name,
            )



    def build(self):
        resindx = 0  # starts at 1 because the 0 is made by the seed coords

        while not self.conformer.is_complete():
            
            fragment_angles = self.angledb.get_fragment(size=self.frag_size)
            
            for residue_angles in fragment_angles:
                
                resindx += 1
                
                for atom in DEFS.NeRF_building_order:
                
                    coords = self.make_coord(
                        residue_angles,
                        atom,
                        resindx,
                        )

                    self.conformer.add_atom_coords(
                        residue_index=resindx,
                        atom_name=atom,
                        coords=coords,
                        )













