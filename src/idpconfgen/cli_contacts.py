"""
Finds CA contacts from PDBs.

PROTOCOL:

1. Reads backbone coordinates (N, CA, C) from PDB files.
2. Calculates CA-CA distances within a chain and between chains.
3. If CA-CA distances are found to be between 5A (default) and the atoms are not
    within 5 residues apart (default). Residues from either side of the CA are saved.
    * If consecutive residues observe 5A interactions, take the residue with
    the closest interaction as the central residue and expand on both sides.
4. Saves results to a JSON dictionary where keys are the input file names
    and the value is a dictionary containing three lists: 'RES1', 'RES2'
    'CA_DIST', 'TOR1', 'TOR2'. Where torsions contain the torsion angles
    from each set of residues.
5. If 'source' JSON file is given, updates that file with the new information.
    Preexisting keys (input file names) are deleted.

CONTROLLED CHECKS:

The same checks are as implemented in `torsions`.
Failed PDBs are registered in `.rpr_on_crash` files and are ignored.

USAGE:
    $ idpconfgen contacts [PDBS]
    $ idpconfgen contacts [PDBS] -sc file.json -o mycontacts.json -n
"""
import argparse
from functools import partial

import numpy as np

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    FileReaderIterator,
    extract_from_tar,
    read_path_bundle,
    read_dictionary_from_disk,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libparse import pop_difference_with_log, values_to_dict
from idpconfgen.libs.libstructure import (
    Structure,
    cols_coords,
    col_name,
    col_resSeq,
    )
from idpconfgen.logger import S, T, init_files, report_on_crash



LOGFILESNAME = '.idpconfgen_contacts'
TMPDIR = '__tmpcontacts__'

_name = 'contacts'
_help = 'Calculates CA contact distances within 5A for PDB files.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-dist',
    '--distance',
    help=(
        "Maximum distance in angstroms allowed for a contact between "
        "residues. Defaults to 5."
        ),
    default=5,
    type=int,
    )

ap.add_argument(
    '-mlen',
    '--max-length',
    help=(
        "Maximum length of consecutive residues stored with the "
        "distance restraint for one point of contact."
        "Defaults to 5."
        ),
    default=5,
    type=int,
    )

libcli.add_argument_source(ap)
libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)

ap.add_argument(
    '-tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )


def ca_distance_matrix(pdb, bound=6):
    """
    Find the residues that are in contact with each other.

    Parameters
    ----------
    pdb : Path
        Path to a PDB file of interest.
    
    bound : int or float
        Maximum distance allowed between CA to be considered a contact.

    Returns
    -------
    compiled_contacts : dict
        Dictionary of residues, torsion angles, and CA distances
    """
    pdb_struc = Structure(pdb)
    pdb_struc.build()
    ca_arr = np.array(
        [pdb_struc.data_array[i] \
        for i, data in enumerate(pdb_struc.data_array[:, col_name]) \
        if data == 'CA']
        )
    ca_coordinates = ca_arr[:, cols_coords].astype(float)
    
    num_residues = len(ca_coordinates)
    distance_matrix = []
    for i in range(num_residues):
        for j in range(i+1, num_residues):
            ca_dist = np.linalg.norm(np.array(ca_coordinates[i]) - np.array(ca_coordinates[j]))
            # Euclidian distance must be within range (default is 6 A)
            # residues must be at least 5 apart
            # TODO: for more than one match, pick CA with shortest contact and
            # take 2 (X) resides worth of information from either side of it
            if ca_dist <= bound and j > i + 4:
                distance_matrix.append((i, j))
    
    return distance_matrix


def main(
    pdb_files,
    distance=5,
    max_length=5,
    source=None,
    output=None,
    ncores=1,
    tmpdir=TMPDIR,
    func=None,
    ):
    """
    Execute main client logic.

    Parameters
    ----------
    pdb_files : str or Path, required
        Location for PDB files to operate on, can be within a folder or inside .TAR file

    distance : int, optional
        Maximum distance in angstroms to be considered a contact.
        Defaults to 5.
    
    max_length : int, optional
        Maximum length of residues stored per identified contact.
        Defaults to 5.
    
    source : string or Path, optional
        If given, updates a preexisting torsions.JSON file.
        Defaults to `None`.

    output : string or Path, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.

    ncores : int
        The numbers of cores to use.
    """
    # validates before performing time consuming calculations
    if source and not source.suffix == '.json':
        raise ValueError('Source file should have `.json` extension.')

    output = output or 'contacts.json'
    if not output.endswith('.json'):
        raise ValueError('Output file should have `.json` extension.')

    log.info(T('Extracting torsion angles'))
    init_files(log, LOGFILESNAME)
    
    if source:
        database_dict = read_dictionary_from_disk(source)

    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))

if __name__ == '__main__':
    libcli.maincli(ap, main)
