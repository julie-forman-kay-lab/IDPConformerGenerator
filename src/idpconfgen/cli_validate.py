"""
Conformer validator.

USAGE:
    $ idpconfgen validator FOLDER
"""
import argparse
import sys
from functools import partial

import numpy as np

from idpconfgen import log
from idpconfgen.core.definitions import heavy_atoms, vdW_radii_dict
from idpconfgen.libs import libcli
from idpconfgen.libs.libcli import CSV2Tuple
from idpconfgen.libs.libio import FileReaderIterator
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libvalidate import (
    validate_bb_bonds_len_from_disk,
    validate_conformer_from_disk,
    bb_bond_length,
    )
from idpconfgen.libs.libstructure import Structure, col_resSeq, col_name, col_resName
from idpconfgen.libs.libcalc import calculate_sequential_bond_distances

from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.validate'
VALIDATION_PROTOCOLS = ('vdw', 'bbl', 'bbd')

_name = 'validate'
_help = 'Validate IDP conformers according to criteria.'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-atm',
    '--atoms-to-consider',
    help=(
        'Which atoms to consider in the validation process. '
        f'Defaults to heavy atoms: {heavy_atoms!r}.'
        ),
    nargs='+',
    default=None,
    )

ap.add_argument(
    '-ele',
    '--elements-to-consider',
    help=(
        'Which elements to consider in the validation process. '
        f'Defaults to heavy atoms: {heavy_atoms!r}.'
        ),
    nargs='+',
    default=tuple(heavy_atoms),
    )

ap.add_argument(
    '-vdw',
    '--vdW-radii',
    help='The van der Waals radii set.',
    choices=list(vdW_radii_dict.keys()),
    default='tsai1999',
    type=str,
    )

ap.add_argument(
    '-vo',
    '--vdW-overlap',
    help='VDW overlap in all atoms.',
    default=0.0,
    type=float,
    )

ap.add_argument(
    '-ra',
    '--residues-apart',
    help=(
        'How many residues apart to evaluate for clashes. '
        'Defaults to 3.'
        ),
    default=3,
    type=int,
    )

ap.add_argument(
    '-bt',
    '--bond-tolerance',
    help='Bond tolerance in angstroms.',
    type=float,
    default=0.01,
    )

ap.add_argument(
    '-v',
    '--validations',
    help=(
        'Which validations to perform. Input as a list. '
        f'Available options are: {VALIDATION_PROTOCOLS}.'
        ),
    nargs='+',
    default=VALIDATION_PROTOCOLS,
    )


libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        atoms_to_consider=None,
        elements_to_consider=None,
        func=None,
        ncores=1,
        residues_apart=3,
        vdW_radii='tsai1999',
        vdW_overlap=0.0,
        bond_tolerance=0.01,
        validations=VALIDATION_PROTOCOLS,
        ):
    """Perform main logic."""
    log.info(T('Validating conformers'))
    init_files(log, LOGFILESNAME)

    if 'vdw' in validations:
        validate_vdW_clashes(
            pdb_files,
            atoms_to_consider,
            elements_to_consider,
            residues_apart,
            vdW_radii,
            vdW_overlap,
            ncores,
            )

    if 'bbl' in validations:
        validate_bond_lengths(
            pdb_files,
            bond_tolerance,
            ncores,
            )

    if 'bbd' in validations:
        bb_bond_length_distribution(
            pdb_files,
            ncores,
            )


def bb_bond_length_distribution(pdb_files, ncores):
    """."""

    pdbs2operate = FileReaderIterator(pdb_files, ext='*.pdb')

    name, pdb_data = next(pdbs2operate)
    s = Structure(pdb_data)
    s.build()
    s.add_filter_backbone(minimal=True)
    fa = s.filtered_atoms
    number_of_bb_atoms = fa.shape[0]

    # ddata from distribution data
    # the number of bonds is the number of atoms -1
    ddata = np.zeros(
        (len(pdbs2operate), number_of_bb_atoms - 1),
        dtype=np.float32,
        )

    ddata[0, :] = calculate_sequential_bond_distances(s.coords)

    labels = np.array([
        f'{s}{r}{a}'
        for r, s, a in fa[:-1, [col_resSeq, col_resName, col_name]]
        ])

    print(labels.shape)

    del s

    execute = partial(report_on_crash, bb_bond_length)
    execute_pool = pool_function(
        partial(starunpack, execute),
        pdbs2operate,
        ncores=ncores,
        )

    for i, lengths in enumerate(execute_pool, start=1):
        ddata[i, :] = lengths

    mean = np.mean(ddata, axis=0)
    std = np.std(ddata, axis=0)
    median = np.median(ddata, axis=0)
    var = np.var(ddata, axis=0)

    results = np.stack([mean, std, median, var], axis=-1)

    table = np.append(labels[:, np.newaxis], results.astype('|S7'), axis=1)

    print(table.shape)

    np.savetxt(
        'validate_bb_bond_distribution.txt',
        table,
        fmt='%s',# + ['%.3f'] * 4,
        delimiter='\t',
        comments='#',
        header=(
            f' Backbone bond lengths distributions for conformers in:\n'
            f' {pdb_files}\n'
            f' analyzed a total of {ddata.shape[0]} conformers\n'
            f' with {ddata.shape[1]} bonds.\n'
            f' mean, std, median, variance'
            ),
        )

    return




def validate_bond_lengths(
        pdb_files,
        bond_tolerance,
        ncores,
        ):
    """."""
    log.info(T('validation for bond lengths'))

    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')

    execute = partial(
        report_on_crash,
        validate_bb_bonds_len_from_disk,
        ROC_prefix=_name,
        tolerance=bond_tolerance,
        )

    execute_pool = pool_function(
        partial(starunpack, execute),
        pdbs2operate,
        ncores=ncores,
        )

    results = list(execute_pool)

    counts, reports = generate_general_validate_report(results)

    with open('validation_report_bb_bond_length.txt', 'w') as fout:
        cstring = '\n'.join(counts)
        rstring = '\n'.join(reports)
        fout.write(
            '# Used parameters:\n'
            f'tolerance (angstroms): {bond_tolerance}\n'
            '\n# Number of conflicts per conformer:\n'
            f"{cstring}"
            '\n\nSpecific conflicts per conformer:\n'
            f"{rstring}"
            )

    return


def validate_vdW_clashes(
        pdb_files,
        atoms_to_consider,
        elements_to_consider,
        residues_apart,
        vdW_radii,
        vdW_overlap,
        ncores,
        ):
    """."""
    log.info(T('validation for van der waals clashes'))

    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')

    # consider placing this at argparse level?
    if len(pdbs2operate) == 0:
        log.error(T('invalid input.'))
        log.error(S(
            'There are no PDB files in the specified folders. '
            'Nothing to do...\n'
            ))
        sys.exit(2)

    execute = partial(
        report_on_crash,
        validate_conformer_from_disk,
        ROC_prefix=_name,
        atoms_to_consider=atoms_to_consider,
        elements_to_consider=elements_to_consider,
        residues_apart=residues_apart,
        vdW_radii=vdW_radii,
        vdW_overlap=vdW_overlap,
        )

    execute_pool = pool_function(
        partial(starunpack, execute),
        pdbs2operate,
        ncores=ncores,
        )

    results = list(execute_pool)

    counts, reports = generate_general_validate_report(results)

    with open('validation_report_vdW_clashes.txt', 'w') as fout:
        cstring = '\n'.join(counts)
        rstring = '\n'.join(reports)
        fout.write(
            '# Used parameters:\n'
            f'Atoms considered: {atoms_to_consider}\n'
            f'Elements considered: {elements_to_consider}\n'
            f'Residue spacing: {residues_apart}\n'
            f'vdW_radii: {vdW_radii}\n'
            '\n# Number of vdW clashes per conformer:\n'
            f"{cstring}"
            '\n\nSpecific clashes per conformer:\n'
            f"{rstring}"
            )


def generate_general_validate_report(results):
    """."""
    results.sort(key=lambda x: x[0])
    max_length = max(len(x[0].stem) for x in results)

    counts = []
    reports = []
    for name, error_number, report_details in results:
        counts.append(
            f'{name}:{" " * (max_length - len(name.stem))} '
            f'{error_number}'
            )
        reports.append(f'# {name}\n{report_details}\n')

    return counts, reports


if __name__ == '__main__':
    libcli.maincli(ap, main)
