"""
Create a database of covalent bond geometries for the backbone.

Given a PDB file:

1) reads each of its trimers, and for the middle residue:
2) Calculates phi/psi and rounds them to the closest 10 degree bin
3) assigns the planar angles found for that residue to the
    trimer/torsion key.
4) updates that information in dictionary (library)

You can provide a folder with several PDB/mmCIF files.

Created key:values have the following form (example) in the output JSON.

{
    'AAA:10,-30': {
        'Cm1_N_Ca': [],
        'N_Ca_C': [],
        'Ca_C_Np1': [],
        'Ca_C_O': [],
        }
    }

Where `AAA` is the trimer peptide in the structure, `10:-30` are the
PHI and PSI torsion angles found for the middle A, and the others are the
planar bond angles around the referred atoms. Cm1 is the carbonyl atom
for the first A, Np1 is the amide nitrogen of the third residue.

Angles are stored in radians and in sync with the expectations of the
`$ idpconfgen build` command. That is, angles are in the form of:

    (PI - angle) / 2


**NOTE:** the `-c` flag will convert the above representation into the one
discribed below. The actual form required by `$ idpconfgen build` is the one
below.

Main keys for the bond types are exact.

{
    'Cm1_N_Ca': {               <- bond type
        'A': {                  <- central residue
            'DE': {             <- previous and next residues (DAE)
                '10,-40': [     <- PHI and PSI angles rounded to 10º bins
                    ...         <- observed bend angles in (pi - radians) / 2
                    ],
                },
            },
        },
    'N_Ca_C': {                 <- other examples
        'G': {
            'PL': {
                '50,-110'[],
                },
            },
        },
    'Ca_C_Np1': {...},
    'Ca_C_O': {...},
        }
    }

USAGE:
    $ idpconfgen bgeo [PDBS]
    $ idpconfgen bgeo [PDBS] -c
    $ idpconfgen bgeo [PDBS] --convert
"""
import argparse
from collections import defaultdict

from idpconfgen import log
from idpconfgen.core.exceptions import PDBFormatError
from idpconfgen.libs import libcli
from idpconfgen.libs.libhigherlevel import (
    convert_bond_geo_lib,
    read_trimer_torsion_planar_angles,
    )
from idpconfgen.libs.libio import FileReaderIterator, save_dict_to_json
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = 'idpconfgen_bgeo'

_name = 'bgeo'
_help = 'Create a database of covalent bond geometry.'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-c',
    '--convert',
    help='Convert Bond Geometry DB to bond type hierarchy.',
    action='store_true',
    )


def main(pdb_files, convert=False, func=None):
    """Perform main logic."""
    log.info(T('Creating bond geometry database.'))
    init_files(log, LOGFILESNAME)

    log.info(T('reading input paths'))
    pdbs = FileReaderIterator(pdb_files, ext='.pdb')
    log.info(S('done'))

    bond_geo_db = defaultdict(dict)

    for pdb_path, pdb_bites in pdbs:
        log.info(S(f'Reading... {pdb_path.stem}'))
        try:
            read_trimer_torsion_planar_angles(pdb_bites, bond_geo_db)
        except PDBFormatError as err:
            log.info(str(err))
            continue

    if convert:
        log.info(S('Converting'))
        converted = convert_bond_geo_lib(bond_geo_db)
        save_dict_to_json(converted, output='bgeo.json')

    else:
        save_dict_to_json(bond_geo_db, output='bgeo.json')

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
