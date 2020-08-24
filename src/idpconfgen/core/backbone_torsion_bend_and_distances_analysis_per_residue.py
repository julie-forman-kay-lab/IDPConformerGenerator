from pathlib import Path
import sys
from statistics import mean, stdev
from functools import partial
import numpy as np
from collections import defaultdict

from idpconfgen.libs.libstructure import Structure, col_name, cols_coords, col_resName
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libio import read_path_bundle

from idpconfgen.libs.libtimer import ProgressCounter
from idpconfgen.core.definitions import aa3to1


NORM = partial(np.linalg.norm, axis=1)
ANGLES = defaultdict(dict)
DISTANCES = defaultdict(dict)


def get_angle(v1, v2):
    cosang = np.diagonal(np.dot(v1, v2.T))
    sinang = np.linalg.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def get_angle2(v1, v2):
    cosang = np.diagonal(np.dot(v1, v2.T)) / (NORM(v1) * NORM(v2))
    return np.arccos(cosang)

def calculate_angle(pdb):
    bbatoms = {'N', 'CA', 'C', 'O'}
    s = Structure(pdb)
    s.build()
    s.add_filter(lambda x: x[col_name] in bbatoms)  # not needed
    fa = s.filtered_atoms

    N = fa[fa[:, col_name] == 'N'][:, cols_coords].astype(float)
    CA = fa[fa[:, col_name] == 'CA'][:, cols_coords].astype(float)
    C = fa[fa[:, col_name] == 'C'][:, cols_coords].astype(float)
    O = fa[fa[:, col_name] == 'O'][:, cols_coords].astype(float)

    CA_N = N - CA
    CA_C = C - CA
    C_CA = CA - C
    C_Np1 = N[1:] - C[:-1]
    C_O = O[:-1] - C[:-1]  # ignores last O which can be a carboxyl
    N_Cm1 = C[:-1] - N[1:]
    N_CA = CA - N

    a1 = get_angle2(CA_N, CA_C)
    a2 = get_angle2(C_CA[:-1], C_Np1)
    a3 = get_angle2(C_CA[:-1], C_O)  # last not used not to confound with carboxyl
    a4 = get_angle2(C_Np1, C_O)
    a5 = get_angle2(N_Cm1, N_CA[1:])

    d1 = NORM(CA_N)
    d2 = NORM(CA_C)
    d3 = NORM(C_Np1)
    d4 = NORM(C_O)

    for residue in aa3to1.keys():

        N_mask = fa[:, col_name] == 'N'
        mask = fa[N_mask, col_resName] == residue

        ANGLES['N_CA_C'].setdefault(residue, []).extend(a1[mask])
        ANGLES['CA_C_Np1'].setdefault(residue, []).extend(a2[mask[:-1]])
        ANGLES['CA_C_O'].setdefault(residue, []).extend(a3[mask[:-1]])
        ANGLES['Np1_C_O'].setdefault(residue, []).extend(a4[mask[:-1]])
        # the last angle blongs to the left residue (from first to ignore minus 1)
        # because it is placed on omega torsion which belongs to thefirst residue
        ANGLES['Cm1_N_CA'].setdefault(residue, []).extend(a5[mask[:-1]])

        DISTANCES['N_CA'].setdefault(residue, []).extend(d1[mask])
        DISTANCES['CA_C'].setdefault(residue, []).extend(d2[mask])
        DISTANCES['C_Np1'].setdefault(residue, []).extend(d3[mask[:-1]])
        DISTANCES['C_O'].setdefault(residue, []).extend(d4[mask[:-1]])

    return #a1, a2, a3, a4, a5, d1, d2, d3, d4


#a1, a2, a3, a4, a5 = [], [], [], [], []
#d1, d2, d3, d4 = [], [], [], []

pdbs = sys.argv[1] #Path('ssplit')
bbatoms = {'N', 'CA', 'C', 'O'}
pdb_files = read_path_bundle([pdbs]) #2* [Path(sys.argv[1])] #list(pdbs.glob('*.pdb'))
#print(pdb_files)


#pool = pool_function(calculate_angle, pdb_files, ncores=15)
#pool = [calculate_angle(pdb) for pdb in pdb_files]
for pdb_file in pdb_files:
    calculate_angle(pdb_file)

#for n1, n2, n3, n4, n5, m1, m2, m3, m4 in pool:
#    a1.extend(n1)
#    a2.extend(n2)
#    a3.extend(n3)
#    a4.extend(n4)
#    a5.extend(n5)
#    d1.extend(m1)
#    d2.extend(m2)
#    d3.extend(m3)
#    d4.extend(m4)
#
#
#print(a1, a2, a3, a4, a5, d1, d2, d3, d4)
#a1 = np.array(a1)
#a2 = np.array(a2)
#a3 = np.array(a3)
#a4 = np.array(a4)
#a5 = np.array(a5)
#d1 = np.array(d1)
#d2 = np.array(d2)
#d3 = np.array(d3)
#d4 = np.array(d4)
#
#print(a1, a2, a3, a4, a5, d1, d2, d3, d4)


_ = '+-'


for key, resdict in ANGLES.items():
    print('*********************')
    print(key)
    for residue, values in resdict.items():
        calc_mean = np.mean(values)
        calc_std = np.std(values)
        print(
            len(values),
            residue,
            np.round(calc_mean, decimals=3),
            np.round(calc_std, decimals=3),
            "|||||",
            calc_mean.as_integer_ratio(),
            "|||||",
            calc_std.as_integer_ratio(),
            )

for key, resdict in DISTANCES.items():
    print('****DISTANCES****')
    print(key)
    for residue, values in resdict.items():
        calc_mean = np.mean(values)
        calc_std = np.std(values)
        print(
            len(values),
            residue,
            np.round(calc_mean, decimals=3),
            np.round(calc_std, decimals=3),
            "|||||",
            calc_mean.as_integer_ratio(),
            "|||||",
            calc_std.as_integer_ratio(),
            )


#print('######################')
#for key, resdict in ANGLES.items():
#    print('*********************')
#    print(key)
#    for residue, values in resdict.items():
#        print(
#            residue, 'average',
#            np.mean(values).as_integer_ratio(),
#            np.std(values).as_integer_ratio(),
#            )
#
#for key, resdict in DISTANCES.items():
#    print('****DISTANCES****')
#    print(key)
#    for residue, values in resdict.items():
#        print(residue, 'average',
#        np.mean(values).as_integer_ratio(), 
#        np.std(values).as_integer_ratio())
#
#
#print('N_CA_C    average',  np.mean(a1).as_integer_ratio(),  _,  np.std(a1))
#print('CA_C_Np1  average',  np.mean(a2).as_integer_ratio(),  _,  np.std(a2))
#print('CA_C_O    average',  np.mean(a3).as_integer_ratio(),  _,  np.std(a3))
#print('Np1_C_O   average',  np.mean(a4).as_integer_ratio(),  _,  np.std(a4))
#print('Cm1_N_CA  average',  np.mean(a5).as_integer_ratio(),  _,  np.std(a5))
#print('disntace  N_CA',     np.mean(d1).as_integer_ratio(),  _,  np.std(d1))
#print('disntace  CA_C',     np.mean(d2).as_integer_ratio(),  _,  np.std(d2))
#print('disntace  C_Np1',    np.mean(d3).as_integer_ratio(),  _,  np.std(d3))
#print('distance  C_O',      np.mean(d4).as_integer_ratio(),  _,  np.std(d4))
