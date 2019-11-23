from collections import defaultdict, deque
import copy
import logging
import math
import os
from pathlib import Path
import pickle
import pprint
import random
import sys
import numpy as np
from clash_validator import ClashValidator


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_ch = logging.StreamHandler()
_ch.setLevel(logging.DEBUG)
_ch.setFormatter(logging.Formatter('%(message)s'))

log.addHandler(_ch)

aaletter_3to1 = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'}

aaletter_1to3 = {}
for K in aaletter_3to1.keys():
    aaletter_1to3[aaletter_3to1[K]] = K


class CustomException(Exception):
    emsg = 'Custom Exception Error!'
    def __init__(self, *args):
        self.args = args

    def __str__(self):
        return self.emsg.format(*args)


class SequenceIndexError(CustomException):
    emsg = 'Conformer sequence does not contain index {} (0-indexed)'


def convert_one2three(input_sequence):
    """
    Convert 1 letter code sequence (FASTA) in a 3 letter code list.
    
    AAA -> ['ALA', 'ALA', 'ALA']
    
    Parameters
    ----------
    input_sequence : str
        Input sequence in FASTA format.

    Returns
    -------
    List
        List of 3 letter code version of the input sequence.
    """
    return [aaletter_1to3[amino_acid] for amino_acid in input_sequence]


# ALAA: this will be the entry point of the database we are building
# now it just loads the Robert's loop database. We can make it load
# whatever.
class LoopDataBase:
    """
    Load loop database.

    Loop database is that of Robert's format.
    
    Parameters
    ----------
    database_file : pickle
    
    Attributes
    ----------
    db : database
    """
    def __init__(self, database_file):
        self.load_database(database_file)
    
    def load_database(self, database_file):
        """
        Loads loop database.

        Database is in the format:

        12asA  P L   101  229   98  -79.590   -2.188  179.809
        12asA  D L   102  228   99 -103.843   11.688  162.288
        12asA  E L   103  227  100  -58.624  134.132 -167.362
        12asA  D L   104  226  101  -85.815  -40.242  172.548
        
        Parameters
        ----------
        database_file : pickle
            Pickle file
        """
        self.db = pickle.load(open(database_file, 'rb'))
    
    def get_pure_fragment(self):
        """Retrieve a random fragment from the fragment database."""
        #return self.db[0]
        return random.sample(self.db, 1)[0]

    def get_fragment(self):
        """
        Select a random element from population.

        In our case selects a random loop from loop DB.
        """
        frag = self.get_pure_fragment()
        return self.transform_frag2dict(frag)
    
    def transform_frag2dict(self, fragment):
        """
        Transform a fragment form the fragment database (:attr:`self.db`)
        to a res:angle dictionary.
        
        Input format:

        12asA  P L   101  229   98  -79.590   -2.188  179.809
        12asA  D L   102  228   99 -103.843   11.688  162.288
        12asA  E L   103  227  100  -58.624  134.132 -167.362
        12asA  D L   104  226  101  -85.815  -40.242  172.548
        
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

        for i, line in enumerate(fragment):
            ls = line.split()
            build_angles[i - 1]['PHI'] = math.radians(float(ls[6]))
            build_angles[i]['PSI'] = math.radians(float(ls[7]))
            build_angles[i]['OMEGA'] = math.radians(float(ls[8]))
        else:
            self.last = build_angles.pop(i)

        if is_first:
            # EYES OPEN: this pops key -1, NOT the last position
            build_angles.pop(-1)

        return build_angles


class Atom:
    """
    Defines an atom building block.

    Atoms building blocks are used during the confomer construction.
    
    The offset attributes (`_off`) refer to the residues from which
    the frame atoms will be read according to the Natural Extension
    Reference Frame building protocol. 0 means the reference atom is taken
    from the same residue where the Atom is going to the inserted.
    Likewise, -1 means from the previous residue and 1 from the next
    residue.

    See :class:`ConformerBuilder`

    Attributes
    ----------
    db_name
        The name of the atom in the Rosetta Database.
        This only deppends on the input.

    name
        The name of the atom when built in the conformer.

    parentoff
        The residue index where the parent atom will be taken.

    xoff
        The residue index where the xaxis atom wil be taken.
    
    yoff
        The residue index where the yaxis atom wil be taken.
    """
    def __init__(self, db_name, name, parentoff, xoff, yoff):
        self.db_name = db_name
        self.name = name
        self.poff = parentoff
        self.xoff = xoff
        self.yoff = yoff
    
    def __repr__(self):
        kwargs = ', '.join(f'{key}={val!r}' for key, val in self.__dict__.items())
        rpr = '{}({})'.format(
            __class__.__name__,
            kwargs,
            )
        return rpr

# Backbone atoms are defined here.
Natom = Atom('UPPER', 'N', -1, -1, -1)
Oatom = Atom('O', 'O', 0, 0, 1)
CAatom = Atom('CA', 'CA', 0, -1, -1)
Catom = Atom('C', 'C', 0, 0, -1)

backbone_seed = [
    {
        'N': np.array([0.000, 0.000, 0.000]),
        'CA': np.array([1.458, 0.000, 0.000]),
        'C': np.array([2.009, 1.420, 0.000]),
        },
    ]


class StateHolder:

    def __init__(self):
        self.state = deque()
    
    def save(self, state):
        """Save something to the register."""
        self.state.append(copy.deepcopy(state))

    def load(self):
        """Returns the last saved item."""
        return self.state.pop()
    
    def clean(self):
        """
        Clean the whole register.

        ATTENTION: All states will be erased and there is no way
        back!!
        """
        self.state = deque()


class ConformerBuilder:
    """
    Conformer builder.
    
    Builds atoms on a coordinate system.

    Parameters
    ----------
    conformer : :class:`Conformer`
        The initial conformer upon which to build new atoms.
        The initial conformer must have at least the three main
        back bone atoms (N, CA, C) to seed the building.i

    sequence : str or list
        The sequence of the complete form of the conformer to build.
    
    basedb : dict
        Base angle database containing the theoretical values for
        positions not given during the building step.
        In our workcases this data base was extracted from Rosetta.
    """
    _bbatoms = set(('N', 'CA', 'O', 'C'))
    
    def __init__(self, conformer, sequence, basedb, angledb):
        
        self.conformer = conformer
        self.seq = sequence
        self.register = StateHolder()
        self.basedb = basedb
        self.angledb = angledb
    
    @property
    def seq(self):
        """The full sequence of the conformer being built."""
        return self._seq
    
    @seq.setter
    def seq(self, sequence):
        if isinstance(sequence, str):
            self._seq = convert_one2three(sequence)
        elif isinstance(sequence, list) and all(len(c) == 3 for c in sequence):
            self._seq = sequence
        else:
            raise ValueError(f'sequence not valid: {sequence}')
        log.info(f'Builder seq length: {len(self._seq)}')
    
    def is_bb_complete(self):
        """
        True if:
            - there are as many residues in the coordinate system
                as residues in the complete sequence
            - and, each coordinate system has at least the backbone
                N, CA, C, O atoms.
        False otherwise
        """
        bbatoms = self._bbatoms

        is_bb_complete = [
            len(conformer) == len(self.seq),
            all(
                bbatoms.issubset(res_coords.keys())
                    for res_coords in self.conformer.coords
                ),
            ]
        return all(is_bb_complete)
    
    def get_residue_type(self, pos):
        """
        Return the residue type for a given pos.

        If pos is -1 returns the residue type of the last residue
        registered in the :attr:`coords`. Otherwise return the type
        of that position if that position already in coords list.

        Raises
        ------
        ValueError
            If pos is higher than length of built coords list.
        """
        if pos < 0:
            # the `+ pos` is equivalent to do `- abs(pos)`
            seq_index = len(self.conformer.coords) + pos
            try:
                return self.seq[seq_index]
            except IndexError:
                raise SequenceIndexError(seq_index)
        
        elif pos <= len(self.coords):
            return self.seq[pos]

        else:
            emsg = \
                'pos must be negative or not higher than the Conformer length'
            raise ValueError(emsg)
    
    def add_coord(self, angles, atom, pos=-1, new_res=False):
        """
        Adds an atom to the system.

        Builds XYZ coordinates for that atom for the `pos` residue.
        
        Parameters
        ----------
        angles : dict
            A dictionary containing the angle information to build upon.
            Needs to match atom. Normally a dictionary containing
            'PHI', 'PSI' and 'OMEGA' is given.

        atom : instance of :class:`Atom`
            The atom to insert 

        pos : int, optional.
            Sequence position where to insert atom.
            If position already exists, overwrites it.
            Defaults to -1, that is, the last position.

        new_res : bool
            Whether the atom is the first of a new residue.
        """
        #self.register.save(self.conformer)
        #log.info(f'building atom: {atom.name}')
        #log.info(f'angles: {angles}')
        #log.info(f'atoms: {atom}')
        
        if new_res:
            self.conformer.add_coord_system({})
        try:
            residue = self.get_residue_type(pos)
        except SequenceIndexError as e:
            self.conformer.pop()
            raise e
        else:
            self.conformer.add_residue(residue)

        coordata = self.basedb[residue][atom.db_name]
        theta = angles.get(coordata.polar_theta, coordata.polar_theta)

        parent_coord = self.conformer.get_coord(pos + atom.poff, coordata.parent_atom)
        xaxis_coord = self.conformer.get_coord(pos + atom.xoff, coordata.xaxis_atom)
        yaxis_coord = self.conformer.get_coord(pos + atom.yoff, coordata.yaxis_atom)
 
        coord = makecoord(
            theta,
            coordata.polar_phi,
            coordata.polar_r,
            parent_coord,
            xaxis_coord,
            yaxis_coord,
            )
       

        self.conformer.add(pos, atom.name, coord)
        return  # None
    
    def add_COO(self):
        log.info('adding COO')
        coordata = self.basedb[self.seq[-1]][Oatom.name]

        # the strings '' and '1' are given to identify
        # the COO C-term atoms. '' is required to match the
        # is_bb_complete method so that at least one O has the same
        # nomenclature of CO carbonyl.
        for i, phi in zip(
                ['', '1'],
                [-120.0, 120.0],
                ):

            parent_coord = self.conformer.get_coord(-1 + Oatom.poff, coordata.parent_atom)
            xaxis_coord = self.conformer.get_coord(-1 + Oatom.xoff, coordata.xaxis_atom)
            yaxis_coord = self.conformer.get_coord(-1 + Oatom.yoff, coordata.yaxis_atom)
            
            coord = makecoord(
                0,
                phi,
                coordata.polar_r,
                parent_coord,
                xaxis_coord,
                yaxis_coord,
                )

            self.conformer.add(-1, Oatom.name + i, coord)

        return  # None
    
    def build_bb(self):
        n = 0
        while not self.is_bb_complete() and n < 100_000:
            n += 1


            build_angles = self.angledb.get_fragment()
            # Python now gives you ordered values =)
            # order of addition is remembered
            for residue_angles in build_angles.values():
                try:
                    self.add_coord(
                        residue_angles,
                        Natom,
                        new_res=True,
                        )
                except SequenceIndexError:
                    log.info('Reached the end of the sequence.')
                    log.info(f'Is sequence complete?: {self.is_bb_complete()}')
                    self.add_COO()
                    return

                self.add_coord(
                    residue_angles,
                    Oatom,
                    pos=-2,
                    )

                self.add_coord(
                    residue_angles,
                    #build_angles[i + 1],
                    CAatom,
                    )
                
                self.add_coord(
                    residue_angles,
                    Catom,
                    )

            try:
                # last state
                last_conformer = self.register.load()
            except IndexError:
                # register is empty
                self.register.save(self.conformer)
                continue

            validator = ClashValidator()
            clash_found = validator.clash_found(self.conformer.coords, last_conformer.coords)
            if not clash_found:
                # no clashes, save this loop
                self.register.save(self.conformer)
            else:
                self.register.save(last_conformer)
                print("here", len(last_conformer.coords), len(self.conformer.coords))
                self.conformer = last_conformer

        else:
            if n > 100_000:
                raise StopIteration('run above 100k cycles!!!!')
    
    def save(self, filename='conformer_gen.pdb'):
        self.conformer.save(filename)


class Conformer:
    """
    A molecule coordinate system.
    
    Parameters
    ----------
    coords_seed : dict
        The coordinates to initiate the conformer with.
    """
    def __init__(self, coords_seed, seq=('MET',)):
        
        self.coords = copy.deepcopy(coords_seed)
        self.seq = list(seq)
        return
    
    def __len__(self):
        return len(self.coords)
   
    @property
    def n_atoms(self):
        """ The total number of atoms in the conformer.""" 
        return sum(len(residue_coords) for residue_coords in self.coords)

    def get_coord(self, pos, atomname):
        """
        Return the coordinates of the `atomname` at position `pos`.
        """
        if atomname in self.coords[pos]:
            return self.coords[pos][atomname]

        elif atomname == 'UPPER':
            return self.coords[pos]['N']
        
        # in our current routines this line is never called
        elif atomname == 'LOWER':
            return self.coords[pos]['C']
        
        else:
            raise ValueError(f'Cant retrieve atom for {atomname}')
    
    def add(self, index, atom_name, coord):
        """
        Add an atom coordinate at given index.
        """
        self.coords[index][atom_name] = coord
    
    def add_coord_system(self, coords=None):
        """
        Appends a coordinate system to the conformer.

        If `None` given, adds an empty system (dictionary)
        """
        if coords:
            self.coords.append(coords)
        else:
            self.coords.append({})
    
    def add_residue(self, res):
        # here we can add all kind of checks
        self.seq.append(res)

    def pop(self, index=-1):
        """
        Removes and returns the coords system at index position.
        """
        return self.coords.pop(index)
    
    def save(self, filename='conformer.pdb'):
        """
        Saves the conformer to file.
        
        File format deduced from extension.

        Parameters
        ----------
        filename : str
            File path and name to write the conformer to.
            Current available formats: pdb.

            .. todo: mmCIF format.
        """
        s = StructureIO(filename, self)
        s.save()

class StructureIO:
    
    def __new__(cls, filename, conformer):
        
        suf = Path(filename).suffix
        
        if suf == '.pdb':
            return PDBIO(filename, conformer)
        elif suf == '.cif':
            return CIFIO(filename, conformer)

class PDBIO:
    # This is a very simple PDB output parser, definitively is not complete
    # here in can extend our preferences as much as required.
    # the important thing is that it is encapsulated alread.
    template_line = "ATOM {:>6}  {:2s}  {:>3} {:>1}{:>4}     {:7.3f} {:7.3f} {:7.3f}                           "  # noqa
    
    def __init__(self, filename, conformer):
        self.filename = filename
        self.conformer = conformer

    def save(self):
    
        ostring = self.gen_string()
        with open(self.filename, 'w') as fh:
            fh.write('\n'.join(ostring) + '\n')
        log.info(f'saved conformer: {self.filename}')
    
    def gen_string(self):
        ostring = []
        atom_count = 1

        for resnum, residue in enumerate(self.conformer.coords):
            for atom, xyz in residue.items():

                ostring.append(self.template_line.format(
                    atom_count,
                    atom,
                    self.conformer.seq[resnum],
                    "A",  # chain ID
                    resnum + 1,
                    xyz[0],
                    xyz[1],
                    xyz[2]
                    ))
                atom_count += 1

        return ostring

class CIFIO:
    pass
#########################
#
# Exciting Math Functions NEVER USED NEVER USED NEVER USED
#
########################## 


# NEVER USED 
def torsion_angle(A, B, C, D):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = A
    p1 = B
    p2 = C
    p3 = D
    
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


# NEVER USED
def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

# NEVER USED
def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.
    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.
    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    http://en.wikipedia.org/wiki/Kabsch_algorithm
    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix
    Returns:
    U -- Rotation matrix
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


# NEVER USED
def unit_vector(vector):
    #Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


# NEVER USED
def xyz_byplane( PARENT, XAXIS, PLANE, CHILD ):
    
    C = PDICT[K][1]["N"] # Planar
    B = PDICT[K][1]["CA"] # X-Axis
    A = PDICT[K][1]["C"] # PARENT
    D = PDICT[K+1][1]["N"] # Child

    R, T = RT_to_plane(PARENT, XAXIS, PLANE)

    return np.dot((CHILD+T), R)

#########################
#
# Key functions for algorithm
#
#########################
def make_axis_vectors(A, B, C): # For np.array([x,y,z]) of connected atoms A->B->C
    
    Avec = A-B # vector along A->B bond

    r1cen = A - B
    r2cen = C - B

    Nvec = np.cross(r1cen, r2cen) # vector parralel to plane

    r1cen = A - B
    r2cen = Nvec

    Cvec = np.cross(r1cen, r2cen) # 90 degree angle to Avec, in plane with A,B,C

    
    Alen = np.sqrt(Avec[0]**2+Avec[1]**2+Avec[2]**2)
    if Alen > 0.0:
        Avec /= Alen

    Nlen = np.sqrt(Nvec[0]**2+Nvec[1]**2+Nvec[2]**2)
    if Nlen > 0.0:
        Nvec /= Nlen

    Clen = np.sqrt(Cvec[0]**2+Cvec[1]**2+Cvec[2]**2)
    if Clen > 0.0:
        Cvec /= Clen

    # There might be some kind of, uh... gimbal lock event with rare combinations of A,B & C ?
    # It is... very rare?
    #
    #if Alen == 0.0 or Nlen == 0.0 or Clen == 0.0:
    #    print "ERROR"
    #    exit()
    
    return Nvec, Avec, Cvec


def RT_to_plane( A, B, C ):# PARENT, XAXIS, PLANE
    """
    This function defines the rotation and translation required to put 4 atoms into a planar orientation where the x-axis
    is determined by the vector A->B and the y axis is the normal vector to the plane defined by A-B-C

    0,0,0 is defined by atom A


        R or Rota  == the rotation matrix required to align the plane with the cartesian X/Y plane
                      (such that the Z axis aligns with the normal vector of the plane)
        T or Trans == the translation vector required to move the plane such that the central point
                      is at 0,0,0 in cartesian space
    """

    NVEC, VEC1, VEC2 = make_axis_vectors(A, B, C)

    #b = np.array([ VEC2, -VEC1, NVEC ] )


    b = np.array([ VEC1, -VEC2, NVEC ] )
    v = np.array([ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ])

    #Rota = np.matrix( [[np.dot( v[0], b[0]), np.dot( v[0], b[1]), np.dot( v[0], b[2])],\
    #                  [np.dot( v[1], b[0]), np.dot( v[1], b[1]), np.dot( v[1], b[2])],\
    #                   [np.dot( v[2], b[0]), np.dot( v[2], b[1]), np.dot( v[2], b[2])]])

    Rota = np.array( [[np.dot( v[0], b[0]), np.dot( v[0], b[1]), np.dot( v[0], b[2])],\
                      [np.dot( v[1], b[0]), np.dot( v[1], b[1]), np.dot( v[1], b[2])],\
                      [np.dot( v[2], b[0]), np.dot( v[2], b[1]), np.dot( v[2], b[2])]])


    #Trans = A

    return Rota, A  # Trans

def makecoord( theta, phi, rad, PARENT, XAXIS, YAXIS ):
    NEW = np.array([0.0, 0.0, 0.0])    
    NEW[0] = rad*math.cos(phi)
    NEW[1] = rad*math.sin(phi)*math.cos(theta)
    NEW[2] = rad*math.sin(phi)*math.sin(theta)

    R, T = RT_to_plane(PARENT, XAXIS, YAXIS)
    
    NEW = np.dot(NEW,R.T)

    NEW = NEW + T

    #print( "THERE", NEW )

    return NEW


class CoordData:
    """
    A coordinate data point.

    Used to store theoretical information from database (aka Rosetta).
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
    """
    rosetta_db = {}
    for aa in aaletter_3to1.keys():

        afile = Path(folder, aa).with_suffix(ext).open().readlines()
        
        rosetta_db[aa] = {}

        icoor_filter = filter(
            lambda x: x.startswith('ICOOR_INTERNAL'),
            afile,
            )

        for a in icoor_filter:
            
            current = CoordData()
            l = a.split()

            atom = l[1]
            current.polar_theta = math.radians(float(l[2]))
            current.polar_phi   = math.radians(float(l[3]))
            current.polar_r     = float(l[4])

            current.parent_atom = l[5]
            current.xaxis_atom  = l[6]
            current.yaxis_atom  = l[7]

            # Handles special cases
            if atom == "CA":
                # this replaces the 180.000 in the .params file
                current.polar_phi = math.radians( 58.300 )
                current.polar_theta = "OMEGA"

                current.parent_atom = "N"
                current.xaxis_atom  = "C"
                current.yaxis_atom  = "CA"
            
            elif atom == "UPPER":
                current.polar_theta = "PSI"
                
            elif atom == "C":
                current.polar_theta = "PHI"

            rosetta_db[aa][atom] = current

    return rosetta_db


# NEVER USED
def get_kdtree( UDICT, start, stop ):

    coords = {"C":[],"N":[],"O":[]}

    for i in UDICT.keys():
        if i >= start and i <= stop:

            for ATOM in UDICT[i]:
                if ATOM == "CA":
                    coords["C"].append(UDICT[i][ATOM])
                if ATOM == "C":
                    coords["C"].append(UDICT[i][ATOM])
                if ATOM == "N":
                    coords["N"].append(UDICT[i][ATOM])
                if ATOM == "O":
                    coords["O"].append(UDICT[i][ATOM])

    kdtrees = {}
    for C in coords.keys():
        kdtrees[C] = scipy.spatial.cKDTree( coords[C], leafsize=2 )

    return kdtrees


if __name__ == '__main__':

    loop_pickle = LoopDataBase('/Users/alaashamandy/Desktop/UNIWORK/CSC495/IDPCalcPDBDownloader/alphas/LOOPS.pickle4')
    rosetta_db = read_rosetta_db('/Users/alaashamandy/Desktop/UNIWORK/CSC495/IDPCalcPDBDownloader/alphas/l-caa')
    
    input_sequence = "MDEYSPKRHDVAQLKFLCESLYDEGIATLGDSHHGWVNDPTSAVNLQLNDLIEHIASFVMSFKIKYPDDGDLSELVEEYLDDTYTLFSSYGINDPELQRWQKTKERLFRLFSGEYISTLMKT"
    # input_sequence = "PDEDRLSPLHSVAAA"
    #input_sequence = "PDEDRLSPLHSV"
    #input_sequence = "D"
    
    conformer = Conformer(backbone_seed)
   
    builder = ConformerBuilder(
        conformer,
        input_sequence,
        rosetta_db,
        loop_pickle,
        )
    
    print(len(builder.conformer.coords))
    builder.build_bb()
    print(len(builder.conformer.coords))
    builder.save('conformer_gen.pdb')

