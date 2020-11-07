"""
Conformer validator.

USAGE:
    $ idpconfgen validator FOLDER
"""
import argparse
import sys
from functools import partial
from pathlib import Path
from time import time

import numpy as np

from idpconfgen import log
from idpconfgen.cli_build import generate_vdW_data
from idpconfgen.core.definitions import (
    aa1to3,
    heavy_atoms,
    vdW_radii_dict,
    )

from idpconfgen.core.build_definitions import (
    average_distance_C_Np1,
    std_distance_C_Np1,
    average_distance_CA_C,
    std_distance_CA_C,
    average_distance_N_CA,
    std_distance_N_CA,
    )

from idpconfgen.libs import libcli
from idpconfgen.libs.libcalc import (
    calc_MSMV,
    calc_all_vs_all_dists_square,
    )
from idpconfgen.libs.libio import FileReaderIterator
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libplot import plot_distribution_list
from idpconfgen.libs.libstructure import (
    Structure,
    col_resName,
    generate_backbone_pairs_labels,
    )
from idpconfgen.libs.libvalidate import (
    eval_bb_bond_length_distribution,
    evaluate_vdw_clash_by_threshold_from_disk,
    validate_bb_bonds_len_from_disk,
    )
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.validate'
VALIDATION_PROTOCOLS = ('vdw2', 'vdw', 'bbl', 'bbd')

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
    default=None,
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

    # preparares available validation subroutines
    available_validations = {
        'vdw2': partial(
            vdw2,
            pdb_files,
            vdW_radii_dict[vdW_radii],
            ),
        'vdw': partial(
            eval_vdw_clashes_cli,
            pdb_files,
            atoms_to_consider,
            elements_to_consider,
            residues_apart,
            vdW_radii,
            vdW_overlap,
            ncores,
            ),

        'bbl': partial(
            eval_backbone_bond_lengths_cli,
            pdb_files,
            bond_tolerance,
            ncores,
            ),

        'bbd': partial(eval_bb_bond_length_distribution_cli, pdb_files, ncores),
        }

    # executes requested validations
    for validation in validations:
        available_validations[validation]()


def vdw2(pdb_files, vdwR):
    """Inspect clashes within the structure."""
    count_nonzero = np.count_nonzero
    logical_and = np.logical_and
    cva = calc_all_vs_all_dists_square
    validation_folder = Path('vdW_validation_results')
    validation_folder.mkdir(parents=False, exist_ok=True)

    s = Structure(Path(pdb_files[0]))
    s.build()

    atoms = s.data_array[:, 2]
    res = s.data_array[:, 6]
    res_nums = np.array([int(i) for i in res])
    res_labels = s.data_array[:, 4]
    coords = s.coords

    vdW_sums, vdw_valid = generate_vdW_data(atoms, res_nums, res_labels, vdwR)

    del s

    for pdb in pdb_files:

        pdb_path = Path(pdb)
        s = Structure(pdb_path)
        s.build()
        coords = s.coords

        results = cva(coords)
        clash = logical_and(results < vdW_sums, vdw_valid)
        has_clash = count_nonzero(clash)
        print(f'** Found {has_clash} clashes for {pdb}')

        # lazy operation
        if has_clash:
            clash_summary = []
            clash_summary_append = clash_summary.append
            i = 0
            j = 0
            for r1, t1, a1 in zip(res_nums, res_labels, atoms):
                i += 1
                for r2, t2, a2 in zip(res_nums[i:], res_labels[i:], atoms[i:]):
                    if clash[j]:
                        line = f'{r1} {t1} {a1} ** {r2} {t2} {a2}'
                        clash_summary_append(line)
                    j += 1

            fout = validation_folder.joinpath(
                f'validation_clash_report_{pdb_path.stem}.txt'
                )
            fout.write_text('\n'.join(clash_summary))


# The following are the subroutine functions that prepare for each
# validation routine

def eval_backbone_bond_lengths_cli(
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


def eval_bb_bond_length_distribution_cli(pdb_files, ncores):
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

    # residue types
    residue_masks = {}
    for residue in aa1to3.values():
        mask_ = fa[:, col_resName] == residue
        # consider only residues that exist in the protein
        if np.any(mask_):
            residue_masks[residue] = mask_

    labels = generate_backbone_pairs_labels(fa)

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
    N_CA_vals = calc_MSMV(ddata[:, ::3])
    CA_C_vals = calc_MSMV(ddata[:, 1::3])
    C_N_vals = calc_MSMV(ddata[:, 2::3])

    # float format
    ffmt = '{:.5f}'
    FF = ffmt.format
    NCAs = '\t'.join(FF(_) for _ in N_CA_vals)
    CACs = '\t'.join(FF(_) for _ in CA_C_vals)
    CNs = '\t'.join(FF(_) for _ in C_N_vals)

    bond_type_report = (
        f' N  - CA: {NCAs}\n'
        f' CA - C : {CACs}\n'
        f' C  - N : {CNs}\n'
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

    # prepares data for plotting

    plt_kwargs = {
        'title': ['N - CA bond', 'CA - C bond', 'C - N bond'],
        'subplots': 3,
        'xlabel': r'bond length ($\AA$)',
        'ylabel': 'Bond',
        'vert': False,
        'usermean': [average_distance_N_CA, average_distance_CA_C, average_distance_C_Np1],
        'userstd': [std_distance_N_CA, std_distance_CA_C, std_distance_C_Np1],
        }

    plot_distribution_list(
        [ddata[:, i::3] for i in range(3)],
        labels=[labels[i::3] for i in range(3)],
        figname='bb_bond_dist.pdf',
        **plt_kwargs,
        )

    res_N = []
    res_CA = []
    res_C = []
    for _residue, mask in residue_masks.items():
        res_N.append(ddata[:, mask[:-1]][:, ::3].ravel())
        res_CA.append(ddata[:, mask[:-1]][:, 1::3].ravel())
        res_C.append(ddata[:, mask[:-1]][:, 2::3].ravel())

    res_labels = list(residue_masks.keys())
    res_type_dist = [res_N, res_CA, res_C]

    plot_distribution_list(
        res_type_dist,
        labels=[res_labels] * 3,
        figname='bb_bond_restype_dist.pdf',
        **plt_kwargs,
        )

    return


def eval_vdw_clashes_cli(
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
        evaluate_vdw_clash_by_threshold_from_disk,
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
            f'vdW tolerance overlap: {vdW_overlap}\n'
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
