"""
Client for using the expended database to build complexes.

Build from a database of CA contacts, torsion angles,
and secondary structure information. Database is as
created by `idpconfgen contacts` CLI.

USAGE:
    $ idpconfgen complex -db contacts.json -seq sequence.fasta --plot
"""
import argparse
from functools import partial
from itertools import combinations, product

import numpy as np

from idpconfgen import Path, log
from idpconfgen.components.bgeo_strategies import (
    add_bgeo_strategy_arg,
    bgeo_strategies_default,
    )
from idpconfgen.components.energy_threshold_type import add_et_type_arg
from idpconfgen.components.plots.plotfuncs import plot_contacts_matrix
from idpconfgen.components.residue_tolerance import add_res_tolerance_groups
from idpconfgen.components.sidechain_packing import (
    DEFAULT_SDM,
    add_mcsce_subparser,
    add_sidechain_method,
    )
from idpconfgen.components.xmer_probs import add_xmer_arg
from idpconfgen.core.build_definitions import forcefields
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import make_folder_or_cwd, read_dictionary_from_disk
from idpconfgen.libs.libmultichain import (
    contact_matrix,
    electropotential_matrix,
    extract_interpairs_from_db,
    extract_intrapairs_from_db,
    find_sa_residues,
    )
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libstructure import Structure
from idpconfgen.logger import S, T, init_files, report_on_crash


_file = Path(__file__).myparents()
LOGFILESNAME = '.idpconfgen_complex'

_name = 'complex'
_help = 'Using the extended database to build intra- and intermolecular complexes.'  # noqa: E501

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_idb(ap)
libcli.add_argument_seq(ap)

ap.add_argument(
    '-fld',
    '--folded-structure',
    help=(
        "Input .PDB or .CIF file for folded structure of interest. "
        "If given folded-structure, only surface residues will be "
        "considered for intermolecular contacts with the sequence. "
        "Please remove all waters before input."
        ),
    default=None,
    )

ap.add_argument(
    '-nc',
    '--nconfs',
    help='Number of conformers to build.',
    default=1,
    type=int,
    )

ap.add_argument(
    '-ph',
    help=(
        'The pH condition for calculating electrostatic potentials.'
        'Defaults to 7.'
        ),
    default=7.0,
    type=float,
    )

ap.add_argument(
    '--blend-weight',
    help=(
        'Integer weight from 0-100. Where 100 only uses the electrostatic '
        'potential contacts frequency and a value of 0 only uses the '
        'sequence-based contacts frequency. '
        'Defaults to 50.'
        ),
    default=50,
    type=int,
    )

#########################################
libcli.add_argument_dloopoff(ap)
libcli.add_argument_dhelix(ap)
libcli.add_argument_dstrand(ap)
libcli.add_argument_dany(ap)
libcli.add_argument_duser(ap)
#########################################

_ffchoice = list(forcefields.keys())
FFDEFAULT = _ffchoice[0]
ap.add_argument(
    '-ff',
    '--forcefield',
    help=(
        'Forcefield parameters and atom labels. '
        f'Defaults to {_ffchoice[0]}.'
        ),
    choices=_ffchoice,
    default=FFDEFAULT,
    )

add_bgeo_strategy_arg(ap)

ap.add_argument(
    '-etbb',
    '--energy-threshold-backbone',
    help=(
        'The energy (kJ) threshold above which fragments will be rejected '
        'when building the BACKBONE atoms. Defaults to 100.'
        ),
    default=100.0,
    type=float,
    )

ap.add_argument(
    '-etss',
    '--energy-threshold-sidechains',
    help=(
        'The energy (kJ) threshold above which conformers will be rejected '
        'after packing the sidechains (ignored if `-dsd`). '
        'Defaults to 250.'
        ),
    default=250.0,
    type=float,
    )

add_et_type_arg(ap)
add_xmer_arg(ap)
add_res_tolerance_groups(ap)

ap.add_argument(
    '-el',
    '--energy-log',
    help='File where to save the energy value of each conformer.',
    type=Path,
    default='energies.log',
    )


add_sidechain_method(ap)
add_mcsce_subparser(ap)
libcli.add_argument_output_folder(ap)
ap.add_argument(
    '--plot-name',
    nargs='?',
    help='Change the name of the heatmap plot saved on disk.',
    default="complex_heatmap.png",
    )
libcli.add_argument_random_seed(ap)
libcli.add_argument_ncores(ap)


def main(
        input_seq,
        database,
        folded_structure=None,
        dloop_off=False,
        dstrand=False,
        dhelix=False,
        duser=False,
        dany=False,
        func=None,
        forcefield=FFDEFAULT,
        bgeo_strategy=bgeo_strategies_default,
        bgeo_path=None,
        residue_tolerance=None,
        nconfs=1,
        ncores=1,
        ph=7,
        blend_weight=50,
        random_seed=0,
        xmer_probs=None,
        output_folder=None,
        energy_log='energies.log',
        sidechain_method=DEFAULT_SDM,
        plot_name="complex_heatmap.png",
        **kwargs,  # other kwargs target energy function, for example.
        ):
    """
    Execute main client logic.

    Distributes over processors.
    """
    db = read_dictionary_from_disk(database)
    db_lst = [{k: v} for k, v in db.items()]
    db_inter_lst = []
    for seg, values in db.items():
        if "inter" in values:
            inter_dict = {}
            inter_dict[seg] = values
            for contact in values["inter"]:
                s1 = next(iter(contact))
                s2 = contact[s1][0]
                inter_dict[s2] = db[s2]
            db_inter_lst.append(inter_dict)
        continue
    
    output_folder = make_folder_or_cwd(output_folder)
    init_files(log, Path(output_folder, LOGFILESNAME))
    
    if type(input_seq) is dict:
        if len(input_seq) > 1:
            log.info(T('multiple sequences detected. assuming they are in a multi-chain complex'))  # noqa: E501
            for chain in input_seq:
                log.info(S(f'{chain}: {input_seq[chain]}'))
        else:
            seq = list(input_seq.values())[0]
            log.info(S(f'input sequence: {seq}'))
    else:
        log.info(S(f'input sequence: {input_seq}'))
        input_seq = {"IDP": input_seq}
    
    log.info(T("analyzing intramolecular contacts from the database"))
    intra_consume = partial(extract_intrapairs_from_db)
    intra_execute = partial(
        report_on_crash,
        intra_consume,
        ROC_exception=Exception,
        ROC_folder=output_folder,
        ROC_prefix=_name,
        )
    intra_execute_pool = pool_function(intra_execute, db_lst, ncores=ncores)
    intra_contacts = []
    for result in intra_execute_pool:
        if result is False:
            continue
        intra_contacts.append(result)
    log.info(S("done"))
    
    log.info(T("analyzing intermolecular contacts from the database"))
    inter_consume = partial(extract_interpairs_from_db)
    inter_execute = partial(
        report_on_crash,
        inter_consume,
        ROC_exception=Exception,
        ROC_folder=output_folder,
        ROC_prefix=_name,
        )
    inter_execute_pool = pool_function(
        inter_execute,
        db_inter_lst,
        ncores=ncores,
        )
    inter_contacts = []
    for result in inter_execute_pool:
        if result is False:
            continue
        inter_contacts.append(result)
    log.info(S("done"))
    
    all_contacts = intra_contacts + inter_contacts
    all_contacts_filtered = [c for c in all_contacts if c]
    
    in_seqs = list(input_seq.values())
    combo_seqs = [c for c in combinations(in_seqs, 2)]
    if folded_structure:
        log.info(T("Folded structure found, identifying surface accessible resiudes"))  # noqa: E501
        assert folded_structure.endswith('.pdb') or \
            folded_structure.endswith('.cif')

        fld_struc = Structure(Path(folded_structure))
        fld_struc.build()
        fld_seqs = fld_struc.fasta
        sa_idx = find_sa_residues(folded_structure)
        
        sa_seqs = {}
        for chain, combos in sa_idx.items():
            temp_seq = []
            for c in combos:
                res = ""
                for i in c:
                    res += fld_seqs[chain][i]
                temp_seq.append(res)
            sa_seqs[chain] = temp_seq

        fld_sa_seqs = list(sa_seqs.values())
        log.info("Found the folowing sequences to be surface accessible:")
        for chain in sa_seqs:
            log.info(f"Chain {chain}: {sa_seqs[chain]}")
        log.info(S("done"))
        combo_seqs = list(product(fld_sa_seqs, in_seqs)) + combo_seqs

    if len(all_contacts_filtered) == 0:
        log.info("WARNING: No contacts found. Your database is invalid.")
        log.info("Only using electropotential contact frequencies.")
    else:
        log.info(T("Calculating all-contacts probability matrix"))
        inter_mtxs = []
        if len(combo_seqs) >= 1:
            inter_seqs = []
            # TODO need a way to delineate the indices and chain IDs
            # from folded template
            for combo in combo_seqs:
                c1 = combo[0]
                c2 = combo[1]
                if type(c1) is list:
                    temp_mtxs = np.zeros(shape=(1, len(c2)))
                    temp_seq = ""
                    for seq in c1:
                        log.info(S(f"calculating intermolecular probabilities between {c2} and {seq}"))  # noqa: E501
                        matrix_consume = partial(contact_matrix, sequence=[seq, c2])  # noqa: E501
                        matrix_execute = partial(
                            report_on_crash,
                            matrix_consume,
                            ROC_exception=Exception,
                            ROC_folder=output_folder,
                            ROC_prefix=_name,
                            )
                        matrix_execute_pool = pool_function(
                            matrix_execute,
                            all_contacts_filtered,
                            ncores=ncores,
                            )
                        mtx = np.zeros((len(seq), len(c2)))
                        location_mtx = {}
                        for id, hit, loc in matrix_execute_pool:
                            mtx = np.add(mtx, hit)
                            location_mtx[id] = loc
                        # TODO figure out what to do with the relative
                        # locations of data
                        # temp_mtxs.append((mtx, location_mtx))
                        temp_seq += seq
                        temp_mtxs = np.vstack((mtx, temp_mtxs))
                        log.info(S("done"))
                    # First row will be all zeroes
                    mtx = temp_mtxs[:-1]
                    norm_mtx = (mtx - np.min(mtx)) / (np.max(mtx) - np.min(mtx))
                    inter_mtxs.append(norm_mtx)
                    inter_seqs.append([c2, temp_seq])
                else:
                    log.info(S(f"calculating intermolecular probabilities between {c1} and {c2}"))  # noqa: E501
                    matrix_consume = partial(contact_matrix, sequence=[c1, c2])
                    matrix_execute = partial(
                        report_on_crash,
                        matrix_consume,
                        ROC_exception=Exception,
                        ROC_folder=output_folder,
                        ROC_prefix=_name,
                        )
                    matrix_execute_pool = pool_function(
                        matrix_execute,
                        all_contacts_filtered,
                        ncores=ncores,
                        )
                    mtx = np.zeros((len(c1), len(c2)))
                    location_mtx = {}
                    for id, hit, loc in matrix_execute_pool:
                        mtx = np.add(mtx, hit)
                        location_mtx[id] = loc
                    norm_mtx = (mtx - np.min(mtx)) / (np.max(mtx) - np.min(mtx))  # noqa: E501
                    inter_mtxs.append(norm_mtx)
                    inter_seqs.append([c2, c1])
                    log.info(S("done"))
        
        log.info(T("Calculating intramolecular contacts probability matrix"))
        intra_mtxs = []
        # Don't need `intra_seqs` variable here because `in_seqs` is the same
        for seq in in_seqs:
            matrix_consume = partial(contact_matrix, sequence=seq)
            matrix_execute = partial(
                report_on_crash,
                matrix_consume,
                ROC_exception=Exception,
                ROC_folder=output_folder,
                ROC_prefix=_name,
                )
            matrix_execute_pool = pool_function(
                matrix_execute,
                all_contacts_filtered,
                ncores=ncores,
                )
            mtx = np.zeros((len(seq), len(seq)))
            location_mtx = {}
            for id, hit, loc in matrix_execute_pool:
                mtx = np.add(mtx, hit)
                location_mtx[id] = loc
            norm_mtx = (mtx - np.min(mtx)) / (np.max(mtx) - np.min(mtx))
            intra_mtxs.append(norm_mtx)
    
    log.info(T("Calculating electropotential contact matrix"))
    inter_electro_mtx = []
    if len(combo_seqs) >= 1:
        for combo in combo_seqs:
            c1 = combo[0]
            c2 = combo[1]
            if type(c1) is list:
                tmp_electro_mtx = np.zeros(shape=(1, len(c2)))
                for seq in c1:
                    mtx = electropotential_matrix([seq, c2], ph)
                    tmp_electro_mtx = np.vstack((mtx, tmp_electro_mtx))
                mtx = tmp_electro_mtx[:-1]
                norm_electro_mtx = (mtx - np.min(mtx)) / (np.max(mtx) - np.min(mtx))  # noqa: E501
                inter_electro_mtx.append(norm_electro_mtx)
            else:
                mtx = electropotential_matrix([c1, c2], ph)
                norm_electro_mtx = (mtx - np.min(mtx)) / (np.max(mtx) - np.min(mtx))  # noqa: E501
                inter_electro_mtx.append(norm_electro_mtx)
    intra_electro_mtx = []
    for seq in in_seqs:
        mtx = electropotential_matrix(seq, ph)
        norm_electro_mtx = (mtx - np.min(mtx)) / (np.max(mtx) - np.min(mtx))  # noqa: E501
        intra_electro_mtx.append(norm_electro_mtx)
    log.info(S('done'))

    log.info(T('saving heatmap distribution plots to output directory'))
    if not (plot_name.endswith(".png") or plot_name.endswith(".jpg")):
        plot_path = Path(plot_name)
        plot_name = plot_path.with_suffix('.png')

    if len(all_contacts_filtered) > 0:
        if len(inter_mtxs) >= 1:
            for i, mtx in enumerate(inter_mtxs):
                plot_out = str(output_folder) + f'/inter_{i}_' + plot_name
                plot_contacts_matrix(mtx, inter_seqs[i], plot_out)
        for i, mtx in enumerate(intra_mtxs):
            plot_out = str(output_folder) + f'/intra_{i}_' + plot_name
            plot_contacts_matrix(mtx, in_seqs[i], plot_out)

    if len(inter_electro_mtx) >= 1:
        for e, mtx in enumerate(inter_electro_mtx):
            plot_electro_out = str(output_folder) + f'/electro_inter_{e}_' + plot_name  # noqa: E501
            plot_contacts_matrix(
                mtx,
                inter_seqs[e],
                plot_electro_out,
                title="Contacts Frequency Heatmap (Electrostatic)"
                )
    for e, mtx in enumerate(intra_electro_mtx):
        plot_electro_out = str(output_folder) + f'/electro_intra_{e}_' + plot_name  # noqa: E501
        plot_contacts_matrix(
            mtx,
            in_seqs[e],
            plot_electro_out,
            title="Contacts Frequency Heatmap (Electrostatic)"
            )
    log.info(S('done'))

    log.info(T('Blending contact map frequencies and plotting'))
    blend_weight = np.clip(blend_weight, 0, 100)  # ensure it's between 0-100
    minimized_blend_weight = blend_weight / 100.0
    if len(all_contacts_filtered) > 0:
        if len(inter_mtxs) >= 1:
            for i, mtx in enumerate(inter_mtxs):
                e_mtx = inter_electro_mtx[i]
                plot_out = str(output_folder) + f'/blend_inter_{i}_' + plot_name
                blended_mtx = (1 - minimized_blend_weight) * mtx + minimized_blend_weight * e_mtx  # noqa: E501
                norm_blended_mtx = (blended_mtx - np.min(blended_mtx)) / (np.max(blended_mtx) - np.min(blended_mtx))  # noqa: E501
                plot_contacts_matrix(
                    norm_blended_mtx,
                    inter_seqs[i],
                    plot_out,
                    title=f"Inter Contacts Frequency Heatmap (Blended {100 - blend_weight}:{blend_weight})"  # noqa: #501
                    )
        for i, mtx in enumerate(intra_mtxs):
            e_mtx = intra_electro_mtx[i]
            plot_out = str(output_folder) + f'/blend_intra_{i}_' + plot_name
            blended_mtx = (1 - minimized_blend_weight) * mtx + minimized_blend_weight * e_mtx  # noqa: E501
            norm_blended_mtx = (blended_mtx - np.min(blended_mtx)) / (np.max(blended_mtx) - np.min(blended_mtx))  # noqa: E501
            plot_contacts_matrix(
                norm_blended_mtx,
                in_seqs[i],
                plot_out,
                title=f"Intra Contacts Frequency Heatmap (Blended {100 - blend_weight}:{blend_weight})"  # noqa: #501
                )
    log.info(S('done'))
    '''
    # Calculates the maximum number of contacts for each input sequence
    # TODO modify this to accept multiple sequences
    max_contacts = calculate_max_contacts(input_seq)
    preselected_contacts = {}
    # Pre-select contacts for each sequence based on maximum number and matrix
    for chain in max_contacts:
        chosen_contacts = []
        max_num = max_contacts[chain]
        for _ in range(nconfs):
            chosen_contacts.append(pick_point_from_heatmap(norm_blended_mtx, max_num))  # noqa: E501
        preselected_contacts[chain] = chosen_contacts
    '''


if __name__ == "__main__":
    libcli.maincli(ap, main)
