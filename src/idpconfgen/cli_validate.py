"""
Conformer validator.

USAGE:
    $ idpconfgen validator FOLDER
"""
import argparse
import sys
from functools import partial
from pathlib import Path

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
    bb_bond_length_dist,
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
    # the number of bonds is the number of atoms -1
    number_of_bb_atoms = fa.shape[0] - 1

    # ddata from distribution data
    ddata = np.zeros(
        (len(pdbs2operate), number_of_bb_atoms),
        dtype=np.float32,
        )

    # masks
    N_mask = fa[:, col_name] == 'N'
    CA_mask = fa[:, col_name] == 'CA'
    C_mask = fa[:, col_name] == 'C'

    # preparing labels
    # 5 is 3 from 3-letter aa code and 2 from 'CA' label length
    max_len = len(str(np.sum(N_mask))) + 5
    fmt_res = '{:<' + str(max_len) + '}'

    N_CA_labels = generate_residue_labels(
        fa[:, [col_resSeq, col_resName, col_name]][N_mask],
        fa[:, [col_resSeq, col_resName, col_name]][CA_mask],
        fmt=fmt_res,
        )

    CA_C_labels = generate_residue_labels(
        fa[:, [col_resSeq, col_resName, col_name]][CA_mask],
        fa[:, [col_resSeq, col_resName, col_name]][C_mask],
        fmt=fmt_res,
        )


    C_N_labels = generate_residue_labels(
        fa[:, [col_resSeq, col_resName, col_name]][C_mask][:-1],
        fa[:, [col_resSeq, col_resName, col_name]][N_mask][1:],
        fmt=fmt_res,
        )

    max_label_len = max(len(i) for i in N_CA_labels)

    labels = np.zeros(number_of_bb_atoms, dtype=f'<U{max_label_len}')
    labels[0::3] = N_CA_labels
    labels[1::3] = CA_C_labels
    labels[2::3] = C_N_labels

    del pdbs2operate
    del s

    pdbs2operate = FileReaderIterator(pdb_files, ext='*.pdb')

    execute = partial(
        report_on_crash,
        bb_bond_length_dist,
        ROC_prefix=f'{_name}_bb_dist',
        )
    execute_pool = pool_function(
        partial(starunpack, execute),
        pdbs2operate,
        ncores=ncores,
        )

    for i, lengths in enumerate(execute_pool, start=0):
        ddata[i, :] = lengths

    mean = np.mean(ddata, axis=0)
    std = np.std(ddata, axis=0)
    median = np.median(ddata, axis=0)
    var = np.var(ddata, axis=0)

    results = np.stack([mean, std, median, var], axis=-1)
    table = np.append(labels[:, np.newaxis], results.astype('|S7'), axis=1)

    # results for bond types
    #N_CA_conformer_means = mean[::3]
    #N_CA_mean = np.mean(N_CA_conformer_means)
    #N_CA_std = np.std(N_CA_conformer_means)
    #N_CA_median = np.median(N_CA_conformer_means)
    #N_CA_variance = np.var(N_CA_conformer_means)

    #CA_C_conformer_means = mean[1::3]
    #CA_C_mean = np.mean(CA_C_conformer_means)
    #CA_C_std = np.std(CA_C_conformer_means)
    #CA_C_median = np.median(CA_C_conformer_means)
    #CA_C_variance = np.var(CA_C_conformer_means)

    #C_N_conformer_means = mean[2::3]
    #C_N_mean = np.mean(C_N_conformer_means)
    #C_N_std = np.std(C_N_conformer_means)
    #C_N_median = np.median(C_N_conformer_means)
    #C_N_variance = np.var(C_N_conformer_means)

    #bond_type_report = (
    #    f' N  - CA: {N_CA_mean:.3f}\t{N_CA_std:.3f}\t{N_CA_median:.3f}\t{N_CA_variance:.3f}\n'
    #    f' CA - C : {CA_C_mean:.3f}\t{CA_C_std:.3f}\t{CA_C_median:.3f}\t{CA_C_variance:.3f}\n'
    #    f' C  - N : {C_N_mean:.3f}\t{C_N_std:.3f}\t{C_N_median:.3f}\t{C_N_variance:.3f}\n'
    #    )

    # results for bond types
    N_CA_conformer_means = ddata[:, ::3]
    N_CA_mean = np.mean(N_CA_conformer_means)
    N_CA_std = np.std(N_CA_conformer_means)
    N_CA_median = np.median(N_CA_conformer_means)
    N_CA_variance = np.var(N_CA_conformer_means)

    CA_C_conformer_means = ddata[:, 1::3]
    CA_C_mean = np.mean(CA_C_conformer_means)
    CA_C_std = np.std(CA_C_conformer_means)
    CA_C_median = np.median(CA_C_conformer_means)
    CA_C_variance = np.var(CA_C_conformer_means)

    C_N_conformer_means = ddata[:, 2::3]
    C_N_mean = np.mean(C_N_conformer_means)
    C_N_std = np.std(C_N_conformer_means)
    C_N_median = np.median(C_N_conformer_means)
    C_N_variance = np.var(C_N_conformer_means)

    bond_type_report = (
        f' N  - CA: {N_CA_mean:.5f}\t{N_CA_std:.5f}\t{N_CA_median:.5f}\t{N_CA_variance:.5f}\n'
        f' CA - C : {CA_C_mean:.5f}\t{CA_C_std:.5f}\t{CA_C_median:.5f}\t{CA_C_variance:.5f}\n'
        f' C  - N : {C_N_mean:.5f}\t{C_N_std:.5f}\t{C_N_median:.5f}\t{C_N_variance:.5f}\n'
        )

    np.savetxt(
        'validate_bb_bond_distribution.txt',
        table,
        fmt='%s',
        delimiter='\t',
        comments='#',
        header=(
            f' Backbone bond lengths distributions for conformers in:\n'
            f' {[Path(p).resolve() for p in pdb_files]}\n'
            f' analyzed a total of {ddata.shape[0]} conformers\n'
            f' with {ddata.shape[1]} bonds.\n'
            f' Results for bond type:\n'
            f'{bond_type_report}\n'
            f' mean, std, median, variance'
            ),
        )

    return


def generate_residue_labels(labels1, labels2, fmt='{:<15}'):
    """."""
    empty_join = ''.join
    labels = []
    LA = labels.append
    for R1, R2 in zip(labels1, labels2):
        r1l = fmt.format(empty_join(R1))
        r2l = fmt.format(empty_join(R2))
        LA(f'{r1l} - {r2l}')

    return labels


    #return [
    #    f'{s1}{r1}{a1} - {s2}{r2}{a2}'
    #    for (s1, r1, a1), (s2, r2, a2) in zip(
    #        labels1, labels2)
    #        #data_array[:, [col_resSeq, col_resName, col_name]][mask1],
    #        #data_array[:, [col_resSeq, col_resName, col_name]][mask2],
    #        #)
    #    ]


def validate_bond_lengths(
        pdb_files,
        bond_tolerance,
        ncores,
        ):
    """."""
    log.info(T('validation for bond lengths'))

    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')

    execute = partial(
        report_on_crash, validate_bb_bonds_len_from_disk,
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
