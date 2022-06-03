"""
Appends exact bond lengths (angstroms) and angles (radians) for each residue
in the IDPConformerGenerator database.

PROTOCOL:

1. Reads bacbkone coordinates (N, CA, C) from PDB/mmCIF files.
2. Calculates bond lengths and bend angles from the backbone.
3. Saves results to a JSON dictionary where keys are the input file names
    and the value is a dictionary containing 'N_CA', 'CA_C', 'C_Np1' for bond lengths
    and 'Cm1_N_CA', 'N_CA_C', 'CA_C_Np1' for bond angles with the order of usage in
    OMEGA, PHI, PSI, respectively per residue.
4. If `source` JSON file is given, updates that file with the new information.
    Pre-existing keys are deleted.
    
USAGE:
    $ idpconfgen bgeodb [PDBS]
    $ idpconfgen bgeodb [PDBS] -sc file.json
    $ idpconfgen bgeodb [PDBS] -sc file.json -n
    $ idpconfgen bgeodb [PDBS] -sc file.json -o my_new_bgeodb.json -n
"""
import argparse
import numpy as np
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    FileReaderIterator,
    read_dictionary_from_disk,
    save_dict_to_json,
)
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libhigherlevel import (
    read_trimer_torsion_planar_angles, 
    convert_bond_geo_lib,
)
from idpconfgen.libs.libparse import pop_difference_with_log
from idpconfgen.logger import S, T, init_files, report_on_crash

LOGFILESNAME = '.idpconfgen_bgeodb'

_name = 'bgeodb'
_help = 'Calculate bond lengths and angles for PDB files.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

