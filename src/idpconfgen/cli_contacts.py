"""# noqa: D205 D400 D411
Find intramolecular and intermolecular CA contacts from PDBs.

PROTOCOL:
1. Reads backbone coordinates (N, CA, C) from PDB files.
2. Calculates CA-CA distances within a chain and between chains.
3. If CA-CA distances are found to be between 6A (default) and the atoms are
    not within 5 residues apart (default). Residues from either side of the
    CA are saved.
    * If consecutive residues observe 5A interactions, take the residue with
    the closest interaction as the central residue and expand on both sides.
4. Saves results to a JSON dictionary where keys are the input file names
    and the values are a list of tuples: [(SEQ1, SEQ2, ARR1, ARR2), ()]
    where the values are an array of information for the specific sequence
    in PDB format.
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
import shutil
from functools import partial

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libcomplex import contact_type
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_dictionary_from_disk,
    read_path_bundle,
    save_dict_to_json,
    )
from idpconfgen.libs.libmultichain import (
    calc_interchain_ca_contacts,
    calc_intrachain_ca_contacts,
    group_chains,
    )
from idpconfgen.libs.libmulticore import pool_function
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
        "residues. Defaults to 6 for intrachain and 12 for interchain."
        "Example: --distance 5 10 or --distance 0 10 for only interchain."
        ),
    type=int,
    nargs="+",
    default=[6, 12],
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


def main(
        pdb_files,
        distance=(6, 12),
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
        Location for PDB files to operate on, can be within a folder
        or inside .TAR file.

    distance : tuple, optional
        Maximum distance in angstroms to be considered a contact.
        Defaults to 6 for intrachain and 12 for interchain.
        
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
    
    if len(distance) >= 2:
        log.info(S(
            f"Taking the first two values of {distance} Å for intramolecular"
            " and intermolecular distance limit respectively."
            ))
        inter_dist = distance[1]
    elif len(distance) < 2:
        inter_dist = 12
    intra_dist = distance[0]
    
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
    
    log.info(T(f'Finding intramolecular contacts of Cα within {intra_dist} Å'))
    consume = partial(
        calc_intrachain_ca_contacts,
        max_dist=intra_dist,
        )
    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_prefix=_name,
        )
    contacts_results = {}
    count = 0
    execute_pool = pool_function(execute, pdbs2operate, ncores=ncores)
    for result in execute_pool:
        if result is False:
            continue
        pdbid = result[0]
        contacts = result[1]
        count += result[2]
        contacts_results[pdbid] = contacts
        contacts_results[pdbid] = {}
        contacts_results[pdbid][contact_type[0]] = contacts
    log.info(f'Total intramolecular contacts: {count} found.')
    
    log.info(T(f'Finding intermolecular contacts of Cα within {inter_dist} Å'))
    pdb_grouped = group_chains(pdbs2operate)
    consume = partial(
        calc_interchain_ca_contacts,
        max_dist=inter_dist,
        )
    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_prefix=_name,
        )
    execute_pool = pool_function(execute, pdb_grouped, ncores=ncores)
    count = 0

    for result in execute_pool:
        if result is False:
            continue
        pdbids = result[0]
        count += result[2]
        for _, ids in enumerate(pdbids):
            if ids not in contacts_results:
                contacts_results[ids] = {}
            contacts_results[ids][contact_type[1]] = result[1][ids]
    log.info(f'Total intermolecular contacts: {count} found.')
    log.info(S('done'))

    log.info('Writing new database to disk...')
    if source:
        for key, value in contacts_results.items():
            try:
                database_dict[key].update(value)
            except KeyError:
                continue
        save_dict_to_json(database_dict, output=output)
    else:
        save_dict_to_json(contacts_results, output=output)
    log.info(S('done'))

    if _istarfile:
        shutil.rmtree(tmpdir)


if __name__ == '__main__':
    libcli.maincli(ap, main)
