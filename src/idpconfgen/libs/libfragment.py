import math
import pickle
import random
from abc import ABC, abstractmethod
from collections import namedtuple

from idpconfgen.libs import libutil as UTIL

class ResidueAngle:
    def __init__(self,
            pdbid=None,
            residue=None,
            pdb_res_1=None,
            pdb_res_2=None,
            pdb_res_3=None,
            dssp=None,
            phi=None,
            psi=None,
            omega=None,
            ):
        #
        self.pdbid = pdbid
        self.residue = residue
        self.dssp = dssp
        self.pdb_res_1 = pdb_res_1
        self.pdb_res_2 = pdb_res_2
        self.pdb_res_3 = pdb_res_3
        self.phi = phi
        self.psi = psi
        self.omega = omega
    
    def __eq__(self, other):
        return all(vs == vo
            for vs, vo in zip(self.__dict__.values(), other.__dict__.values()))

    def __str__(self):
        return '{:>5}  {} {} {:>4} {:>4} {:>4} {:>8.3f} {:>8.3f} {:>8.3f}'.format(
            self.pdbid,
            self.residue,
            self.dssp,
            self.pdb_res_1,
            self.pdb_res_2,
            self.pdb_res_3,
            self.phi * 180 / math.pi,
            self.psi * 180 / math.pi,
            self.omega * 180 / math.pi,
            )
    
    def __repr__(self):
        return str(self)


#ResidueAngle = namedtuple(
#    'ResidueAngle',
#    [
#        'pdbid',
#        'letter',
#        'dssp',
#        'phi',
#        'psi',
#        'omega',
#        ]
#    )


class FragmentDBABC(ABC):
    """Provide abstract interface for FragmengDB."""

    @abstractmethod
    def get_pure_fragment(self):
        return


class FragmentAngleDB(FragmentDBABC):
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
    
    Read from files with :method:`from_file`.

    Attributes
    ----------
    db : database
    """
    
    def __len__(self):
        return len(self.db)

    def __eq__(self, other):
        is_eq = all((
            len(self) == len(other),
            self.total_len() == other.total_len(),
            all(sres == ores for sfrags, ofrags in zip(self.db, other.db)\
                for sres, ores in zip(sfrags, ofrags)),
            ))
        return is_eq
    
    def total_len(self):
        """
        The total length of the fragment db.

        Total length counts the sum of all residue entries from all the
        fragments in the db.
        """
        return sum(len(fragment) for fragment in self.db)

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
        sliceObj = UTIL.random_fragment(self.db, fragsize)
        frag = self.get_pure_fragment()
        return frag[sliceObj]
    
    def pickle_db(self, fname='fragment_angle_db.pickle'):
        """
        Save current db to a pickle file.

        Saved pickle file can later be read with .from_file() method.
        """
        with open(fname, 'wb') as fh:
            pickle.dump(self.db, fh)

    @classmethod
    def from_file(cls, fname):
        try:
            data = cls.read_text_file(fname)
            parsed = cls._parse_raw_data(data)
        except UnicodeDecodeError:
            parsed = cls.read_pickle(fname)

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
        with open(fname, 'rb') as fh:
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
                        residue=ls[1],
                        dssp=ls[2],
                        pdb_res_1=ls[3],
                        pdb_res_2=ls[4],
                        pdb_res_3=ls[5],
                        phi=math.radians(float(ls[6])),
                        psi=math.radians(float(ls[7])),
                        omega=math.radians(float(ls[8])),
                        )
                    )
            parsed_data.append(parsed_block)
        return parsed_data

