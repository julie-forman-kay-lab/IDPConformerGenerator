"""Library for the Conformer builder."""
from abc import ABC, abstractmethod
import math
from collections import defaultdict, deque, namedtuple
from typing import NamedTuple

import numpy as np

from idpconfgen import Path
from idpconfgen.core import definitions as DEFS
from idpconfgen.core import interfaces as ITF
from idpconfgen.libs import libfragment as LF
from idpconfgen.libs import libutil as UTIL


class ConformerTemplate:
    """
    A conformer template for building.

    Implements the basic logic of a conformer data structure to store
    and access XYZ coordinates.

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
        
        # this will need to be updated if we incorporated sidechains
        self._atomnames = \
            np.array(
                DEFS.backbone_atoms * len(self.seq) + (DEFS.COO_name,),
                dtype='<U1',
                )
         
        # backbone mask does not account for O from (COO)
        self.bb_mask = np.isin(self.atomnames, DEFS.backbone_atoms)

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
        """
        Evaluate if conformer is complete.

        True if all atoms in conformer coords array have assigned
        coordinates.
        """
        return not np.any(np.isnan(self.coords))

    def is_bb_complete(self):
        """
        True if all backbone atoms are assigned coordinates.

        .. note::
            
            Additional oxygen from terminal carboxyl group does not
            count as backbone atom.

        """
        return not np.any(np.isnan(self.coords[self.bb_mask,:]))
    
    def clean_coords(
            self,
            *,
            from_residue=None,
            from_atom=None,
            to_residue=None,
            to_atom=None,
            ):
        """
        Restore atom coordinates to np.nan.

        Parameters
        ----------
        """
        from_index = self._get_index(from_residue, from_atom)
        to_index = self._get_index(to_residue, to_atom)
        # remember index from None are None
        self.coords[from_index:to_index,:] = np.nan

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
        except ValueError:
            # ValueError occurs when adding the oxygen for the last
            # carboxyl group or when atom_name is None
            if atom_name == DEFS.COO_name:
                return -1
            elif atom_name is None:
                atom_name = 0
            else:
                raise ValueError('atom_name not valid: {!r}'.format(atom_name))
        else:
            try:
                index += residue_index * DEFS.num_bb_atoms
            except TypeError:  # residue_index is None
                return None
            else:
                return index


class AtomNeRF(NamedTuple):
    """
    Defines an atom building block for NeRF algorithm.

    Atoms building blocks are used during the confomer construction.
    
    The offset attributes (`_off`) refer to the residues from which
    the frame atoms will be read according to the Natural Extension
    Reference Frame building protocol. 0 means the reference atom is taken
    from the same residue where the Atom is going to the inserted.
    Likewise, -1 means from the previous residue and 1 from the next
    residue. Same applies for the `resindx` parameter, where indexes
    refer to the building progress: 0 the current residue being built,
    -1 the previously built residue and 1 the following residue.

    There is no differentiation from N to C building or C to N building,
    that is done at the level of the ConformerNeRF and BuilderNeRF.

    See :class:`ConformerBuilder`.

    Attributes
    ----------
    name
        The name of the atom when built in the conformer.

    rosetta_name
        The name this atom is represented in the Rosetta atom
            parameter DB. This relates to the building order.
            For example, N is UPPER in NtoC building and C is LOWER
            in CtoN building.

    poff
        The residue index where the parent atom will be taken.

    xoff
        The residue index where the xaxis atom wil be taken.
    
    yoff
        The residue index where the yaxis atom wil be taken.

    resoff
        The residue index offset relative to the building process
        where to add the atom. See :class:`ConformerNerf`.
    """
    name: str
    rosetta_name: str
    poff: int
    xoff: int
    yoff: int
    resoff: int

# to build forward
N_atom_NeRF = AtomNeRF(DEFS.N_name, 'UPPER', -1, -1, -1, 0)
CA_atom_NeRF = AtomNeRF(DEFS.CA_name, 'CA', 0, -1, -1, 0)
C_atom_NeRF = AtomNeRF(DEFS.C_name, 'C', 0, 0, -1, 0)
O_atom_NeRF = AtomNeRF(DEFS.O_name, 'O', 0, 0, 1, -1)
NeRF_building_order_NtoC = [
    N_atom_NeRF,
    O_atom_NeRF,
    CA_atom_NeRF,
    C_atom_NeRF,
    ]

# to build backwards
N_atom_NeRF_b = AtomNeRF(DEFS.N_name, 'N', 1, 1, 1, 0)
CA_atom_NeRF_b = AtomNeRF(DEFS.CA_name, 'CA', 0, 1, 1, 0)
C_atom_NeRF_b = AtomNeRF(DEFS.C_name, 'LOWER', 0, 0, 1, 0)
O_atom_NeRF_b = AtomNeRF(DEFS.O_name, 'O', 0, 0, -1, 0)
NeRF_building_order_CtoN = [
    C_atom_NeRF_b,
    CA_atom_NeRF_b,
    O_atom_NeRF_b,
    N_atom_NeRF_b,
    ]


class FragDBNeRFFactory:
    """
    Factory for FragmentAngleDBNeRF.

    Generates equal and independent copies of FragmentAngleDBNeRF.
    """
    def __init__(self):
        self._fragdb = None
    
    def read_fragdb(self, fname):
        """
        Reads an angle database from a file.
        """
        self._fragdb = LF.FragmentAngleDB.from_file(fname)
    
    def get_db_for_nerf(self):
        """
        Create a :class:`FragmentAngleDBNeRF`.

        DB is created from a previously loaded fragment angle file.

        Returns the new DB object.
        """
        fragnerf = FragmentAngleDBNeRF(self._fragdb)
        return fragnerf


class FragmentAngleDBNeRF:
    """Data base of angles from loop fragments."""     
    
    def __init__(self, fragmentdb):
        self._fragdb = fragmentdb

    def __eq__(self, other):
        return self._fragdb == other._fragdb

    def get_angle_fragment(self, fragsize=None):
        """
        Return a random fragment from database.
        
        Parameters
        ----------
        fragsize : int or None, optional
            If give, returns a random slice of size `fragsize` from
            the randomly selected fragment.

        Returns
        -------
        dict
            A dictionary containing the residue angle requirements
            to build a full residue with the scope of
            :class:`ConformerBuilderNeRF`.

            .. seealso::
                
                :method:`transform_frag2dict`
        """

        frag = self._fragdb.get_angle_fragment(fragsize)
        return self.transform_frag2dict(frag)

    def transform_frag2dict(self, fragment):
        """
        Transform a fragment form the fragment database (:attr:`self.fragdb`)
        to a res:angle dictionary.
        
        For the N-term residue the first PHI angle is discarded.
        The last PSI and OMEGA angles, are stored and not returned.
        This is such because these angles will only be used if
        an additional fragment is required in the building process.
        Therefore, on the second fragment the first PHI residue is
        combined with the last OMEGA and PSI from the previous fragment.

        This strategy abstracts the building process and avoids the need
        to build the first fragment separately and to rebuild the whole
        conformer from scratch at each fragment addition.
        
        Parameters
        ----------
        fragment
            An angle fragment list given by
            :method:`self.fragdb.get_pure_fragment`.

        Returns
        -------
        dict
            Contains PHI, PSI and OMEGA angles per key that are required
            to build an amino-acid. The returned dictionary has
            one key, starting at 0, for each residue to be built.
            
            {0: {
                'PHI': float,
                'PSI': float,
                'OMEGA': float,
                },
            1: {
                'PHI': float,
                'PSI': float,
                'OMEGA': float,
                },
            }
        """
        build_angles = defaultdict(dict)
        
        try:
            build_angles[-1] = self.last
            is_first = False
        except AttributeError:
            is_first = True

        for i, residue in enumerate(fragment):
            build_angles[i - 1]['PHI'] = residue.phi
            build_angles[i]['PSI'] = residue.psi
            build_angles[i]['OMEGA'] = residue.omega
        else:
            self.last = build_angles.pop(i)

        if is_first:
            # EYES OPEN: this pops key -1, NOT the last position
            build_angles.pop(-1)

        return build_angles


class ConformerBuilderNeRF:
    """
    Conformer builder.
    """
    def __init__(
            self,
            conformer,
            angledb,
            rosettadb,
            frag_size=None,
            reverse_build=False,
            start_residue_index=0,
            ):
        #
        self._conformer = conformer
        self._angledb = angledb
        self._rosettadb = rosettadb
        self.frag_size = frag_size
        self._previous_fragment_index = deque()
        # to clean: self.reverse_build = reverse_build
        self._current_residue_index = start_residue_index

        if reverse_build:  # build from C to N
            self.increment_conformer_residue = self.decrement_residue
            self.building_order = NeRF_building_order_CtoN
            self.get_fragment = self.angledb.get_fragment_reversed
            #
            self._clean_from_residue = None
            self._clean_from_atom = None
            self._clean_to_residue = self.current_residue_index
            self._clean_to_atom = DEFS.N_name

        else:  # build from N to C
            self.increment_conformer_residue = self.increment_residue
            self.building_order = NeRF_building_order_NtoC
            self.get_fragment = self.angledb.get_angle_fragment
            #
            self._clean_from_residue = self.current_residue_index
            self._clean_from_atom = DEFS.O_name
            self._clean_to_residue = None
            self._clean_to_atom = None

    @property
    def current_residue_index(self):
        """The current conformer residue at which the builder is operating."""
        return self._current_residue_index

    @property
    def conformer(self):
        """The conformer upon which the builder is operating."""
        return self._conformer

    @property
    def angledb(self):
        """
        The data base of angle fragments.
        
        :class:`FragmengDB`.
        """
        return self._angledb

    @property
    def rosettadb(self):
        """
        The Rosetta atom building parameters.

        .. seealso::

            :func:`read_rosetta_db`
        """
        return self._rosettadb
    
    @property
    def frag_size(self):
        return self._frag_size

    @frag_size.setter
    def frag_size(self, value):
        # +1 is intrinsic to this protocol
        # see tests/test_libbuilder.TestFragmentAngleDB.test_transform_frag2dict
        # therefore, to build a fragment of 5 residues we need to
        # request a fragment of 6 to the FragmentAngleDB.
        if value > 0:
            self._frag_size = int(value) + 1
        else:
            raise ValueError(f'value must be a positive number: {value}')

    def increment_residue(self):
        self._current_residue_index += 1
    
    def decrement_residue(self):
        self._current_residue_index -= 1
   
    def register_previous_index(self):
        self._previous_fragment_index.append(self._current_residue_index)

    def restore_previous_index(self):
        self._current_residue_index = self._previous_fragment_index.pop()
    
    def restore_two_states(self):
        # because we are looking for speed enhancement I avoid here
        # using condictionals for workflow and instead define a new
        # method for the case to rever back two states
        self._previous_fragment_index.pop()
        self.restore_previous_index()

    def build_backbone(self):
        """Build backbone of conformer."""
        while True: #  not self.conformer.is_bb_complete():

            self.register_previous_index()

            try:
                self._build_backbone_fragment()
            except StopIteration:
                # do the logic on building COO if that is the case
                pass
            
            try:
                self.validate_new_fragment()

            except ClashError:
                self.restore_previous_index()
                self.conformer.clean_coords(
                    from_residue=self.clean_from_residue,
                    from_atom=self.clean_from_atom,
                    to_residue=self.clean_to_residue,
                    to_atom=self.clean_to_atom,
                    )
            except RepetitiveClashError:
                self.restore_two_states()
                self.conformer.clean_coords(
                    from_residue=self.clean_from_residue,
                    from_atom=self.clean_from_atom,
                    to_residue=self.clean_to_residue,
                    to_atom=self.clean_to_atom,
                    )


    def _build_backbone_fragment(self):

        fragment_angles = self.get_fragment(size=self.frag_size)
         
        for residue_angles in fragment_angles:

            # this can be incremented directly in the first run
            # because the first residue of the conformer template
            # is already created by the seed coordinates
            self.increment_conformer_residue()

            for atom_nerf in self.building_order:
                
                if self._valid_to_build(atom_nerf.name):

                    XYZ_coords = self.make_coord(
                        residue_angles,
                        atom_nerf,
                        )

                    self.conformer.add_atom_coords(
                        residue_index=\
                            self.current_residue_index + atom_nerf.resoff,
                        atom=atom_nerf.name,
                        coords=XYZ_coords,
                        )
                else:
                    raise StopIteration

    def _valid_to_build(self, index, atom_name):
        """
        Returns `True` if there are no XYZ coordinates for residue
        `index` and atom `atom_name`, ie, if XYZ are np.nan.
        Returns `False` otherwise.
        """
        return np.all(
            np.isnan(
                self.conformer.get_coord(
                    self.current_residue_index,
                    atom_name,
                    )
                )
            )

    def _make_coord(self, residue_angles, atom):
         
        # example:    self.rosettadb['ALA']['UPPER']
        # returns a RosettaAtomData type
        residue_type = self.conformer.seq[self.current_residue_index]
        #atom_name = self.db_atom_map.get(atom.name, atom.name)

        # rosetta_atom isinstance RosettaAtomData
        rosetta_atom = self.rosettadb[residue_type][atom.rosetta_name]
        
        theta = residue_angles.get(
            rosetta_atom.polar_theta,
            rosetta_atom.polar_theta,
            )

        parent_coord, xaxis_coord, yaxis_coord = \
            self._get_conformer_coords(atom, rosetta_atom)

        coord = LIBCALC.makecoord(
            theta,
            rosetta_atom.polar_phi,
            rosetta_atom.polar_r,
            parent_coord,
            xaxis_coord,
            yaxis_coord,
            )

        return coord

    def _get_coords_from_conformer(self, atom_to_build, rosetta_atom):
        """
        :class:`AtomNeRF`, :class:`RosettaAtomData` -> tuple
            
        Returns tuple of np.array of shape (3,)
        """

        parent_coord = self.conformer.get_coord(
            self.current_residue_index + atom_to_build.poff,
            rosetta_atom.name,
            )

        xaxis_coord = self.conformer.get_coord(
            self.current_residue_index + atom_to_build.xoff,
            rosetta_atom.name,
            )

        yaxis_coord = self.conformer.get_coord(
            self.current_residue_index + atom_to_build.yoff,
            rosetta_atom.name,
            )

        return parent_coord, xaxis_coord, yaxis_coord


class RosettaAtomData(ITF.ReprClean):
    # RosettaAtomData does not inhering from typing.NamedTuple
    # because currently .polar_theta can be either float or str.
    # .. TODO: change this.
    """
    Rosetta atom data for building.

    Required atom information to perform the angle to coordinate
    transformation using the NeRF algorithm. Data was extracted from
    the Rosetta Commons suite.
    """
    def __init__(
            self,
            polar_theta=None,
            polar_phi=None,
            polar_r=None,
            parent_atom=None,
            xaxis_atom=None,
            yaxis_atom=None,
            ):
        self.polar_theta = polar_theta
        self.polar_phi = polar_phi
        self.polar_r = polar_r
        self.parent_atom = parent_atom
        self.xaxis_atom = xaxis_atom
        self.yaxis_atom = yaxis_atom
    
    def __str__(self):
            return repr(self)


def read_rosetta_db(folder, ext='.params'):
    """
    Read ROSETTA DB parameters need for conformer building.

    folder path -> dict
    """
    rosetta_db = {}
    for aa in DEFS.aa3to1.keys():

        afile = Path(DEFS.data_folder, aa).with_suffix(ext).open().readlines()
        
        rosetta_db[aa] = {}

        icoor_filter = filter(
            lambda x: x.startswith('ICOOR_INTERNAL'),
            afile,
            )

        for a in icoor_filter:
            
            current = RosettaAtomData()
            l = a.split()

            atom = l[1]
            current.polar_theta = math.radians(float(l[2]))
            current.polar_phi = math.radians(float(l[3]))
            current.polar_r = float(l[4])

            current.parent_atom = l[5]
            current.xaxis_atom = l[6]
            current.yaxis_atom = l[7]

            # Handles special cases
            if atom == "CA":
                # this replaces the 180.000 in the .params file
                current.polar_phi = math.radians(58.300)
                current.polar_theta = "OMEGA"

                current.parent_atom = "N"
                current.xaxis_atom = "C"
                current.yaxis_atom = "CA"
            
            elif atom == "UPPER":
                current.polar_theta = "PSI"
                
            elif atom == "C":
                current.polar_theta = "PHI"

            rosetta_db[aa][atom] = current

    return rosetta_db

