"""
Calculates exact bond lengths (angstroms) and angles (radians) for the idpcg database.

PROTOCOL:

1. Reads bacbkone coordinates (N, CA, C) from PDB/mmCIF files.
2. Calculates bond lengths and bend angles from the backbone.
3. Saves results to a JSON dictionary where keys are the input file names
    and the value is a dictionary containing 'N_CA', 'CA_C', 'C_Np1', 'C_O' 
    for bond lengths and 'Cm1_N_CA', 'N_CA_C', 'CA_C_Np1', 'Ca_C_O' 
    for bond angles.
4. If `source` JSON file is given, updates that file with the 
    new information. Pre-existing keys are deleted.
    
CONTROLLED CHECKS:

For each PDB/mmCIF analyzed, fails if:

1. The PDB does not start with an N atom.
2. Carbonyl backbone atoms do not follow
    'CA', 'C', 'O', 'CA', 'C', 'O', 'CA' pattern.

Failed PDBs are registered in `.rpr_on_crash` files and ignored.
    
USAGE:
    $ idpconfgen bgeodb [PDBS]
    $ idpconfgen bgeodb [PDBS] -sc file.json
    $ idpconfgen bgeodb [PDBS] -sc file.json -n
    $ idpconfgen bgeodb [PDBS] -sc file.json -o my_new_bgeodb.json -n
"""
import argparse
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    FileReaderIterator,
    read_dictionary_from_disk,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libhigherlevel import cli_helper_calc_bgeo
from idpconfgen.libs.libparse import pop_difference_with_log
from idpconfgen.logger import S, T, init_files, report_on_crash

LOGFILESNAME = '.idpconfgen_bgeodb'

_name = 'bgeodb'
_help = 'Calculate bond lengths and angles per residue for PDB files.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)
libcli.add_argument_source(ap)
libcli.add_argument_output(ap)
libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        source=None,
        output=None,
        ncores=1,
        func=None,
        ):
    """
    Run main script logic.

    Parameters
    ----------
    pdb_files : str or Path, required
        Location for PDB files to operate on, can be within a folder
        or inside a .TAR file
        
    source : string or Path, optional
        If given, updates a preexisting .JSON file.
        Defaults to `None`.
        
    output : string or Path, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.
        
    ncores : int
        The numbers of cores to use.
    """
    if source and not source.suffix == '.json':
        raise ValueError('Source file should have `.json` extension.')
    
    output = output or 'bgeodb.json'
    if not output.endswith('.json'):
        raise ValueError('Output file should have `.json` extension.')

    init_files(log, LOGFILESNAME)
    log.info(T('Extracting bond geometries'))
    
    if source:
        database_dict = read_dictionary_from_disk(source)
    
    log.info(T('reading input paths'))
    pdbs = FileReaderIterator(pdb_files, ext='.pdb')
    log.info(S('done'))
    
    consume = partial(starunpack, cli_helper_calc_bgeo)

    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_prefix=_name,
        )

    execute_pool = pool_function(execute, pdbs, ncores=ncores)
    
    bgeo_result = {
        Path(pdbid).stem: lengths_and_angles
        for pdbid, lengths_and_angles in execute_pool
        }
    
    if source:
        pop_difference_with_log(database_dict, bgeo_result)
        popped_prior = []
        for key, value in bgeo_result.items():
            # where value is a dictionary {'Ca_C_Np1':, 'Ca_C_O':, ...}
            try:
                database_dict[key].update(value)
            except KeyError as e:
                popped_prior.append(str(e)[1:-1])
                
        log.info(S(f'PDB IDs popped during previous steps to initialize the database: {popped_prior}'))
        
        save_dict_to_json(database_dict, output=output)
        
    else:
        save_dict_to_json(bgeo_result, output=output)
    
    log.info(S('done'))
    
    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
