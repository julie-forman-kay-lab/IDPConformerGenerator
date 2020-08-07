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
from idpconfgen.core.definitions import (
    heavy_atoms, vdW_radii_dict, aa1to3,
    distance_N_CA, distance_CA_C, distance_C_Np1,
    distance_N_CA_std, distance_CA_C_std, distance_C_Np1_std
    )
from idpconfgen.libs import libcli
from idpconfgen.libs.libcli import CSV2Tuple
from idpconfgen.libs.libio import FileReaderIterator
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libvalidate import (
    validate_bb_bonds_len_from_disk,
    evaluate_vdw_clash_by_threshold_from_disk,
    eval_bb_bond_length_distribution,
    )
from idpconfgen.libs.libstructure import (
    Structure,
    col_resSeq,
    col_name,
    col_resName,
    generate_residue_labels,
    )
from idpconfgen.libs.libcalc import calculate_sequential_bond_distances
from idpconfgen.libs.libplot import plot_distribution, plot_distribution_list

from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.validate'
VALIDATION_PROTOCOLS = ('vdw', 'bbl', 'bbd')

_name = 'validate'
_help = 'Validate IDP conformers according to criteria.'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

# CLI parameters receive parameters for all the validation tests.
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
        'Which atoms to consider in the VDW clash validation process. '
        f'Defaults to heavy atoms: {heavy_atoms!r}. '
        'Works for --validations vdw.'
        ),
    nargs='+',
    default=tuple(heavy_atoms),
    )

ap.add_argument(
    '-ele',
    '--elements-to-consider',
    help=(
        'Which elements to consider in the VDW validation process. '
        f'Defaults to heavy atoms: {heavy_atoms!r}.'
        'Works for --validations vdw.'
        ),
    nargs='+',
    default=tuple(heavy_atoms),
    )

ap.add_argument(
    '-vdw',
    '--vdW-radii',
    help=(
        'The van der Waals radii set. '
        'Works for --validations vdw.'
        ),
    choices=list(vdW_radii_dict.keys()),
    default='tsai1999',
    type=str,
    )

ap.add_argument(
    '-vo',
    '--vdW-overlap',
    help=(
        'VDW overlap in all atoms. '
        'Works for --validations vdw.'
        ),
    default=0.0,
    type=float,
    )

ap.add_argument(
    '-ra',
    '--residues-apart',
    help=(
        'How many residues apart to start evaluating for VDW clashes. '
        'Defaults to 2. '
        'Works for --validations vdw.'
        ),
    default=2,
    type=int,
    )

ap.add_argument(
    '-bt',
    '--bond-tolerance',
    help=(
        'Bond tolerance in angstroms. Defaults to 0.01. '
        'Works for --validations bbl'
        ),
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
        bond_tolerance=0.01,
        elements_to_consider=None,
        func=None,
        ncores=1,
        residues_apart=3,
        validations=VALIDATION_PROTOCOLS,
        vdW_overlap=0.0,
        vdW_radii='tsai1999',
        ):
    """Validate conformers according to the criteria implemented."""
    log.info(T('Validating conformers'))
    init_files(log, LOGFILESNAME)

    # confirms there are PDBs in the input
    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')
    if len(pdbs2operate) == 0:
        log.error(T('invalid input.'))
        log.error(S(
            'There are no PDB files in the specified folders. '
            'Nothing to do...\n'
            ))
        sys.exit(2)
    del pdbs2operate

    # preparares available validation tests
    available_validations = {
        'vdw': partial(
            validate_vdw_clashes,
            pdb_files,
            atoms_to_consider,
            elements_to_consider,
            residues_apart,
            vdW_radii,
            vdW_overlap,
            ncores,
            ),

        'bbl': partial(
            validate_backbone_bond_lengths,
            pdb_files,
            bond_tolerance,
            ncores,
            ),

        'bbd': partial(eval_bb_bond_length_distribution, pdb_files, ncores),
        }

    # executes requested validations
    for validation in validations:
        available_validations[validation]()


def eval_bb_bond_length_distribution(pdb_files, ncores):
    """Evaluate backbone bond length distribution."""
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

    # residue types
    residue_masks =  {}
    for residue in aa1to3.values():
        mask_ = fa[:, col_resName] == residue
        # consider only residues that exist in the protein
        if np.any(mask_):
            residue_masks[residue] = mask_

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
        eval_bb_bond_length_distribution,
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

    plot_distribution_list(
        [ddata[:, i::3] for i in range(3)],
        title=['N - CA bond', 'CA - C bond', 'C - N bond'],
        subplots=3,
        labels=[labels[i::3] for i in range(3)],
        vert=False,
        figname='bb_bond_dist.pdf',
        ylabel='Bond',
        xlabel='bond length ($\AA$)',
        usermean=[distance_N_CA, distance_CA_C, distance_C_Np1],
        userstd=[distance_N_CA_std, distance_CA_C_std, distance_C_Np1_std],
        )

    res_N = []
    res_CA = []
    res_C = []
    for residue, mask in residue_masks.items():
        res_N.append(ddata[:, mask[:-1]][:, ::3].ravel())
        res_CA.append(ddata[:, mask[:-1]][:, 1::3].ravel())
        res_C.append(ddata[:, mask[:-1]][:, 2::3].ravel())

    res_labels = list(residue_masks.keys())
    res_type_dist = [res_N, res_CA, res_C]

    plot_distribution_list(
        res_type_dist,
        title=['N - CA bond', 'CA - C bond', 'C - N bond'],
        subplots=3,
        labels=[res_labels] * 3,
        vert=False,
        figname='bb_bond_restype_dist.pdf',
        ylabel='Bond',
        xlabel='bond length ($\AA$)',
        usermean=[distance_N_CA, distance_CA_C, distance_C_Np1],
        userstd=[distance_N_CA_std, distance_CA_C_std, distance_C_Np1_std],
        )

    return




def validate_backbone_bond_lengths(
        pdb_files,
        bond_tolerance,
        ncores,
        report_file_name='validation_report_bb_bond_length.txt',
        ):
    """
    Validate backbone bond lengths.

    Acceptable bond lengths are according to definitions in IDPConfGen.
    """
    log.info(T('validating backbone bond lengths'))

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

    counts, reports = generate_general_validation_report(results)
    cstring = '\n'.join(counts)
    rstring = '\n'.join(reports)

    with open(report_file_name, 'w') as fout:
        fout.write(
            '# Used parameters:\n'
            f'tolerance (angstroms): {bond_tolerance}\n'
            '\n# Number of conflicts per conformer:\n'
            f"{cstring}"
            '\n\nSpecific conflicts per conformer:\n'
            f"{rstring}"
            )

    return


def validate_vdw_clashes(
        pdb_files,
        atoms_to_consider,
        elements_to_consider,
        residues_apart,
        vdW_radii,
        vdW_overlap,
        ncores,
        report_file_name='validation_report_vdW_clashes.txt',
        ):
    """
    Validate conformers based on VDW clashes considering a threshold.

    At the current implementation is this function is a helper of the
    the `cli_validate`, therefore all parameters are mandatory.
    """
    log.info(T('validation for van der waals clashes'))

    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')

    execute = partial(
        report_on_crash,
        evaluate_vdW_clashes_from_disk,
        ROC_prefix=f'{_name}_vdw_clash_threshold',
        atoms_to_consider=atoms_to_consider,
        elements_to_consider=elements_to_consider,
        residues_apart=residues_apart,
        vdW_radii=vdW_radii,
        vdW_overlap=vdW_overlap,
        )

    execute_pool = pool_function(
        # starunpack is needed because FileReaderIterator
        partial(starunpack, execute),
        pdbs2operate,
        ncores=ncores,
        )

    results = list(execute_pool)

    counts, reports = generate_general_validation_report(results)
    cstring = '\n'.join(counts)
    rstring = '\n'.join(reports)

    with open(report_file_name, 'w') as fout:
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

    return


def generate_general_validation_report(results):
    """
    Generate a per conformer validation report.

    Firstly, conformers are listed followed by the number of errors
        identified for each.

    Secondly, specific errors are listed for each conformer.

    Returns
    -------
    tuple of lists
        Each list contains the strings for each report type.
    """
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
