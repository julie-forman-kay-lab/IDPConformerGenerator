"""Library for the Conformer builder."""
from abc import ABC, abstractmethod
import math
from collections import namedtuple

import numpy as np

from idpconfgen.core import definitions as DEFS
from idpconfgen.core import interfaces as ITF
from idpconfgen.libs import libutil as UTIL


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
        
        # this will need to be updated if we incorporated sidechains
        self._atomnames = \
            np.array(
                DEFS.backbone_atoms * len(self.seq) + (DEFS.COO_name,),
                dtype='<U1',
                )
         
        # backbone mask does not account for O from (COO)
        self.bb_mask = np.isin(self.atom_names, DEFS.backbone_atoms)

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


class ConformerNeRF(ConformerTemplate):
    """Conformer template specific for NeRF building algorithm."""
    
    def __init__(self, *args, init_residue_index=0, **kwargs):
        self._current_residue_index = init_residue_index
        super().__init__(*args, **kwargs)
    
    @property
    def current_residue_index(self):
        return self._current_residue_index

    def increment_residue(self):
        self._current_residue_index += 1
    
    def decrement_residue(self):
        self._current_residue_index -= 1
    

    def add_coord_NtoC(self, atom, coord):
        """
        Add atom coord considering N-term to C-term building.

        Parameters
        ----------
        atom : :class:`AtomNeRF`
            The atom to add to the conformer.

        coords : Numpy array of shape (3,)
            The XYZ coordinates
        """
        self.add_atom_coords(
            self._current_residue_index + atom.resindx,
            atom.name,
            coord,
            )

    def add_coord_CtoN(self, atom, coord):
        """
        Add atom coord considering C-term to N-term building.

        Parameters
        ----------
        atom : :class:`AtomNeRF`
            The atom to add to the conformer.

        coords : Numpy array of shape (3,)
            The XYZ coordinates
        """
        self.add_atom_coords(
            self._current_residue_index - atom.resindx,
            atom.name,
            coord,
            )


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
            ):
        #
        self._conformer = conformer
        self._angledb = angledb
        self._rosettadb = rosettadb
        self.frag_size = frag_size
        self.reverse_build = reverse_build

        if reverse_build:  # build from C to N
            self.increment_conformer_residue = self._conformer.decrement_residue
            self.building_order = DEFS.NeRF_building_order_CtoN
            self.get_fragment = self.angledb.get_fragment_reversed
        else:  # build from N to C
            self.increment_conformer_residue = self._conformer.increment_residue
            self.building_order = DEFS.NeRF_building_order
            self.get_fragment = self.angledb.get_fragment

        self.db_atom_map = {
            'UPPER': DEFS.N_name,
            'LOWER': DEFS.C_name,
            }


    @property
    def conformer(self):
        return self._conformer

    @property
    def angledb(self):
        return self._angledb
    
    @property
    def rosettadb(self):
        return self._rosettadb

    def build_backbone(self):
        """Build backbone of conformer."""
        while True: #  not self.conformer.is_bb_complete():
            try:
                self._build_backbone()
            except StopIteration:
                pass
    
    def _build_backbone(self):

        fragment_angles = self.get_fragment(size=self.frag_size)
         
        for residue_angles in fragment_angles:
            
            #if self.conformer.is_bb_complete():
            #    raise StopIteration

            # this can be incremented directly in the first run
            # because the first residue of the conformer template
            # is already created by the seed coordinates
            self.increment_conformer_residue()
            
            for atom_nerf in self.building_order:
                
                if self._valid_to_build(
                        self.conformer.current_residue_index_resindx,
                        atom_nerf.name
                        ):

                    coords = self.make_coord(
                        residue_angles,
                        atom_nerf,
                        )

                    self.conformer.add_atom_coords(
                        residue_index= (
                            self.conformer.current_residue_index_resindx
                            + atom_nerg.resindx
                            ),
                        atom=atom_nerf.name,
                        coords=coords,
                        )
                else:
                    raise StopIteration

    def _valid_to_build(self, index, atom_name):
        return np.all(np.isnan(self.conformer.get_coord(index, atom_name)))

    def _make_coord(self, residue_angles, atom):
         
        # example:    self.rosettadb['ALA']['UPPER']
        # returns a RosettaAtomData type
        rosetta_atom = \
            self.rosettadb[self.conformer.current_residue_type][atom.rname]
        
        theta = residue_angles.get(
            rosetta_atom.polar_theta,
            rosetta_atom.polar_theta,
            )

        parent_coord, xaxis_coord, yaxis_coord = self._get_coords(atom)

        coord = LIBCALC.makecoord

    def _make_coord_CtoN(self, residue_angles, atom):
        
        # example:    self.rosettadb['ALA']['UPPER']
        # returns a RosettaAtomData type
        #rosetta_atom = \
        #    self.rosettadb[self.conformer.current_residue_type][atom.rname]
        
        #theta = residue_angles.get(
        #    rosetta_atom.polar_theta,
        #    rosetta_atom.polar_theta,
        #    )

        parent_coord, xaxis_coord, yaxis_coord = \
            self._get_coords(atom, reverse=-1)

        #coord = LIBCALC.makecoord

    def _get_coords(self, atom):
        """atom is a :class:`Atom_NeRF`."""

        parent_coord = self.conformer.get_coord(
            self.conformer.current_residue_index + atom.poff,
            atom.name,
            )

        xaxis_coord = self.conformer.get_coord(
            self.conformer.current_residue_index + atom.xoff,
            atom.name,
            )

        yaxis_coord = self.conformer.get_coord(
            self.conformer.current_residue_index + atom.yoff,
            atom.name,
            )


class FragmentDBABC(ABC):
    
    @abstractmethod
    def get_pure_fragment(self):
        return


class FragmentAngleDBNeRF(FragmentDBABC, ITF.Prototype):
    """
    Database for fragment angles.

    Loads angles from a file (text or pickle) of the following format::

        12asA  P L   101  229   98  -79.590   -2.188  179.809
        12asA  D L   102  228   99 -103.843   11.688  162.288
        12asA  E L   103  227  100  -58.624  134.132 -167.362
        12asA  D L   104  226  101  -85.815  -40.242  172.548

        16pkA  Y L   249  166  249  -83.645  148.597  176.997
        16pkA  S L   250  165  250  -81.829  132.302 -177.161
        16pkA  I L   251  164  251 -116.621    4.225  178.907
        16pkA  G L   252  163  252   57.598 -128.309 -178.929
        16pkA  K L   253  162  253  -97.343   16.681 -177.532
    
    Read from files with :attr:`from_file`.

    Attributes
    ----------
    db : database
    """
    @property
    def db(self):
        """
        Database attribute.
        """
        return self._db

    def get_pure_fragment(self):
        """
        Retrieve a random fragment from the fragment database.
        
        The fragment is in its pure form.

        .. seealso:: :attr:`get_angle_fragment`  
    
        """
        return random.sample(self.db, 1)[0]

    def get_angle_fragment(self, fragsize=None):
        """
        Select a random element from population.

        In our case selects a random loop from loop DB.
        """
        sliceObj = UTIL.random_fragment(fragsize)
        frag = self.get_random_fragment()
        return self._transform_frag2dict(frag[sliceObj])
    
    def _transform_frag2dict(self, fragment):
        """
        Transform a fragment form the fragment database (:attr:`self.db`)
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
            A :attr:`self.db` fragment given by :method:`get_pure_fragment`.

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

        for residue in enumerate(fragment):
            build_angles[i - 1]['PHI'] = residue.phi
            build_angles[i]['PSI'] = residue.psi
            build_angles[i]['OMEGA'] = residue.omega
        else:
            self.last = build_angles.pop(i)

        if is_first:
            # EYES OPEN: this pops key -1, NOT the last position
            build_angles.pop(-1)

        return build_angles
    
    @classmethod
    def from_file(cls, fname):
        try:
            data = cls.read_text_file(fname)
        except UnicodeDecodeError:
            data = cls.from_pickle(fname)

        parsed = cls._parse_raw_data(data)
        c = cls()
        c._db = parsed
        return c
    
    @staticmethod
    def read_pickle(fname):
        """
        Loads loop database from pickle file.
        
        Required format:
            [
                [
                    '12asA  P L   101  229   98  -79.590   -2.188  179.809'
                    '12asA  D L   102  228   99 -103.843   11.688  162.288'
                    '12asA  E L   103  227  100  -58.624  134.132 -167.362'
                    '12asA  D L   104  226  101  -85.815  -40.242  172.548'
                    ],
                (...)
                ]
             
        Parameters
        ----------
        database_file : pickle
            Pickle file
        """
        with open(database_file, 'rb') as fh:
             data = pickle.load(fh)
        return data

    @staticmethod
    def read_text_file(fname):
        """
        Loads database from text file.
        
        Required Format:
            
            12asA  P L   101  229   98  -79.590   -2.188  179.809
            12asA  D L   102  228   99 -103.843   11.688  162.288
            12asA  E L   103  227  100  -58.624  134.132 -167.362
            12asA  D L   104  226  101  -85.815  -40.242  172.548


            12asA  G L   156  174  153   68.074   16.392 -175.411
            12asA  L L   157  173  154  -81.735  119.822 -176.938
            12asA  A L   158  172  155  -71.433  133.117 -172.298
            12asA  I L   165  165  162 -101.494  150.487 -174.989


            16pkA  Y L   249  166  249  -83.645  148.597  176.997
            16pkA  S L   250  165  250  -81.829  132.302 -177.161
            16pkA  I L   251  164  251 -116.621    4.225  178.907
            16pkA  G L   252  163  252   57.598 -128.309 -178.929
            16pkA  K L   253  162  253  -97.343   16.681 -177.532
            16pkA  S L   254  161  254  -64.590  147.146  179.504

        Parameters
        ----------
        fname : str or Path
            Path to text file angle DB.
        """
        with open(fname, 'r') as fh:
            blocks = fh.read().strip().split('\n\n\n')

        data = []
        for block in blocks:
            data.append(block.split('\n'))
       
        return data

    @staticmethod
    def _parse_raw_data(data):
        parsed_data = []
        for block in data:
            parsed_block = []
            for line in block:
                ls = line.split()
                parsed_block.append(
                    ResidueAngle(
                        pdbid=ls[0],
                        letter=ls[1],
                        dssp=ls[2],
                        phi=math.radians(float(ls[6])),
                        psi=math.radians(float(ls[7])),
                        omega=math.radians(float(ls[8])),
                        )
                    )
            parsed_data.append(parsed_block)
        return parsed_data


ResidueAngle = namedtuple(
    'ResidueAngle',
    [
        'pdbid',
        'letter',
        'dssp',
        'phi',
        'psi',
        'omega',
        ]
    )


class RosettaAtomData:
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
    
    def __repr__(self):
        kwargs = ', '.join(
            f'{key}={val!r}' for key, val in self.__dict__.items()
            )
        rpr = '{}({})'.format(
            __class__.__name__,
            kwargs,
            )
        return rpr

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

