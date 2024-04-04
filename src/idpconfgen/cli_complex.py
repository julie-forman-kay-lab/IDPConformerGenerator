"""
Client for using the expended database to build complexes.

Build from a database of CA contacts, torsion angles,
and secondary structure information. Database is as
created by `idpconfgen contacts` CLI.

USAGE:
    $ idpconfgen complex -db contacts.json -seq sequence.fasta --plot
"""
import argparse
import re
from collections import defaultdict
from functools import partial
from itertools import combinations, product
from random import randint, random

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
    contact_type,
    electropotential_matrix,
    extract_interpairs_from_db,
    extract_intrapairs_from_db,
    find_sa_residues,
    process_custom_contacts,
    select_contacts,
    select_custom_contacts,
    )
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libparse import update_chars_lower
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
    '--phos',
    help=(
        "Indicates which residues on which chain in the FASTA file is "
        "phosphorylated, thus adding more negative charges depending on "
        "the pH. Chains are denoted by colons and phosphorylated residues "
        "for that chain are separated by commas. Additional chains are "
        "delimited by slash and pattern must end at a slash. The name "
        "of the chain will correspond to the >Name in the FASTA file. "
        "For e.g. --phos Name:12,14,15/B:13,10/"
        ),
    nargs='?',
    default=None,
    )

ap.add_argument(
    '--custom-contacts',
    help=(
        "Input text (.txt) file containing known contacts between certain "
        "residues in the provided protein sequences/templates. "
        "Each new line in the format file will be a new known contact. "
        "Intra- and inter-contacts are separated by a single slash (/). "
        "The name of the IDP chain will correspond to the >Name in the FASTA "
        "file while single letter chains suggest a folded template. "
        "For example: (B:13,14,15/sic1:1,3,4,5) will indicate an "
        "intermolecular contact between the IDP sequence >sic1 and chain B "
        "of the folded protein template."
        "Folded chains MUST be before IDP sequence as indicated above."
        "If custom-contacts are given, these contacts will be selected with a "
        "probability of 0.9 and knowledge-based + electrostatic potential "
        "contacts will be selected with a probability of 0.1."
        ),
    default=None,
    )

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
    '--blended-contacts-weight',
    help=(
        'Integer weight from 0-100. Where 100 only uses the electrostatic '
        'potential contacts frequency and a value of 0 only uses the '
        'sequence-based contacts frequency. '
        'Defaults to 60.'
        ),
    default=60,
    type=int,
    )

ap.add_argument(
    '--custom-contacts-weight',
    help=(
        'Integer weight from 0-100. Where 100 only uses the custom-contacts '
        'given by the user and a value of 0 only uses the sequence-based '
        'contacts frequency. '
        'Defaults to 90.'
        ),
    default=90,
    type=int,
    )

ap.add_argument(
    '-max',
    '--max-contacts',
    help=(
        'Integer maximum of contacts built between chains. '
        'Defaults to 10.'
        ),
    default=10,
    type=int,
    )

ap.add_argument(
    '--ignore-sasa',
    help="Ignores SASA limit of 31.65 Å² for custom contacts.",
    action='store_true',
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
        phos=None,
        folded_structure=None,
        custom_contacts=None,
        ignore_sasa=False,
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
        blended_contacts_weight=60,
        custom_contacts_weight=90,
        max_contacts=10,
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
    output_folder = make_folder_or_cwd(output_folder)
    init_files(log, Path(output_folder, LOGFILESNAME))
    
    log.info(T('Reading extended database from disk'))
    db = read_dictionary_from_disk(database)
    db_lst = [{k: v} for k, v in db.items()]
    db_inter_lst = []
    for seg, values in db.items():
        if contact_type[1] in values:
            inter_dict = {}
            inter_dict[seg] = values
            for contact in values[contact_type[1]]:
                s1 = next(iter(contact))
                s2 = contact[s1][0]
                inter_dict[s2] = db[s2]
            db_inter_lst.append(inter_dict)
        continue
    log.info(S('done'))
    
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
        input_seq = {1: input_seq}
    
    if phos:
        REGEX_PATTERN = re.compile(r"(.*:([0-9,])*[/])+")
        phos_ptm = {}
        if REGEX_PATTERN.match(phos):
            phos_residues = phos.split('/')
            phos_residues.pop()  # last element should be empty
            for c in phos_residues:
                chain = c.split(':')
                residues = chain[1].split(',')
                residues.pop()
                if len([input_seq.values()]) == 1:
                    phos_ptm[1] = [int(r) for r in residues]
                else:
                    phos_ptm[f">{chain[0]}"] = [int(r) for r in residues]
        else:
            log.info(S('Incorrect pattern input for --phos.'))
            log.info(S('Pattern is as follows: A:1,3,4,/B:3,9,5,/'))
         
        mod_in_seqs = []
        for key, sequence in input_seq.items():
            phos_idx = [x - 1 for x in phos_ptm[key]]
            mod_in_seqs.append(update_chars_lower(sequence, phos_idx))
        combo_mod_seqs = [c for c in combinations(mod_in_seqs, 2)]
    
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
    in_res = []
    for seq in in_seqs:
        in_res.append([i for i in range(1, len(seq))])
    in_chains = list(input_seq.keys())
    combo_res = [c for c in combinations(in_res, 2)]
    combo_seqs = [c for c in combinations(in_seqs, 2)]
    combo_chains = [c for c in combinations(in_chains, 2)]
    if folded_structure:
        log.info(T("Folded structure found, identifying surface accessible resiudes"))  # noqa: E501
        assert folded_structure.endswith('.pdb') or \
            folded_structure.endswith('.cif')

        fld_struc = Structure(Path(folded_structure))
        fld_struc.build()
        fld_res = fld_struc.residues_splitted
        fld_seqs = fld_struc.fasta
        sa_idx = find_sa_residues(folded_structure)
        
        sa_seqs = {}
        sa_res = {}
        for chain, combos in sa_idx.items():
            temp_seq = []
            temp_res = []
            for c in combos:
                res = ""
                real_res = []
                for i in c:
                    res += fld_seqs[chain][i]
                    real_res.append(fld_res[chain][i])
                temp_seq.append(res)
                temp_res.append(real_res)
            sa_seqs[chain] = temp_seq
            sa_res[chain] = temp_res
        fld_sa_res = list(sa_res.values())
        fld_sa_seqs = list(sa_seqs.values())
        fld_sa_chains = list(sa_seqs.keys())
        log.info("Found the folowing sequences to be surface accessible:")
        for chain in sa_seqs:
            log.info(f"Chain {chain}: {sa_seqs[chain]}")
        log.info(S("done"))
        # Combo seqs and combo chains are aligned to give locations
        combo_res = list(product(fld_sa_res, in_res)) + combo_res
        combo_seqs = list(product(fld_sa_seqs, in_seqs)) + combo_seqs
        combo_chains = list(product(fld_sa_chains, in_chains)) + combo_chains
        if phos:
            combo_mod_seqs = list(product(fld_sa_seqs, mod_in_seqs)) + combo_mod_seqs  # noqa: E501
    # NOTE: based on chain information, see if any of the residues in
    # custom contacts are in the ones with folded protein chains.
    # If custom contacts not within surface residues, log it for the user.
    
    # We do not change "combo_chains" since that's our anchor for alignment.
    # If we do not have custom contacts for specific combo chains, store None.
    
    # We select later if we have custom contacts, then 90% of the time we will
    # select a custom_res/custom_seq. If not, then 10% of the time we select
    # from our knowledge driven database.
    if custom_contacts is not None:
        cus_inter_res, cus_intra_res = process_custom_contacts(
            custom_contacts,
            combo_chains,
            combo_res,
            input_seq,
            ignore_sasa=ignore_sasa,
            )
    input_seq_keys = list(input_seq.keys())

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
            log.info(S(f"calculating intramolecular contact probabilities of {seq}"))  # noqa: E501
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
    if phos:
        combo_seqs = combo_mod_seqs
        in_seqs = mod_in_seqs
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
            plot_contacts_matrix(mtx, in_seqs[i].upper(), plot_out)

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
            in_seqs[e].upper(),
            plot_electro_out,
            title="Contacts Frequency Heatmap (Electrostatic)"
            )
    log.info(S('done'))

    log.info(T('Blending contact map frequencies and plotting'))
    # Ensure it's between 0-100
    blended_inter_mtxs = []
    blended_intra_mtxs = []
    blend_weight = np.clip(blended_contacts_weight, 0, 100)
    minimized_blend_weight = blend_weight / 100.0
    if len(all_contacts_filtered) > 0:
        if len(inter_mtxs) >= 1:
            for i, mtx in enumerate(inter_mtxs):
                e_mtx = inter_electro_mtx[i]
                plot_out = str(output_folder) + f'/blend_inter_{i}_' + plot_name
                blended_mtx = (1 - minimized_blend_weight) * mtx + minimized_blend_weight * e_mtx  # noqa: E501
                norm_blended_mtx = (blended_mtx - np.min(blended_mtx)) / (np.max(blended_mtx) - np.min(blended_mtx))  # noqa: E501
                blended_inter_mtxs.append(norm_blended_mtx)
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
            blended_intra_mtxs.append(norm_blended_mtx)
            plot_contacts_matrix(
                norm_blended_mtx,
                in_seqs[i].upper(),
                plot_out,
                title=f"Intra Contacts Frequency Heatmap (Blended {100 - blend_weight}:{blend_weight})"  # noqa: #501
                )
    log.info(S('done'))
    
    # Custom and knowledge based contacts are delineated by datatype.
    # Custom contacts are lists while knolwedge based ones are tuple.
    # TODO function to return a distribution of distances given
    # residue pairs and the information we have in our database.
    selected_contacts = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))  # noqa: E501
    if custom_contacts:
        custom_contacts_weight = np.clip(custom_contacts_weight, 0, 100)
        min_contacts_weight = custom_contacts_weight / 100.0
        log.info(T(f'Choosing at most {max_contacts} contacts for every conformer. Custom-contacts will be chosen with a probability of {custom_contacts_weight} %'))  # noqa: E501
        for _ in range(nconfs):
            for i, norm_inter_mtx in enumerate(blended_inter_mtxs):
                case = combo_chains[i]
                contacts_counter = randint(1, max_contacts)
                while contacts_counter - 1 > 0:
                    custom = random() < min_contacts_weight
                    if not custom or cus_inter_res[i] is None:
                        x_coords, y_coords = select_contacts(
                            coords=norm_inter_mtx,
                            max_num_points=contacts_counter
                            )
                        selected_contacts['X'][contact_type[1]][case].append(x_coords)  # noqa: E501
                        selected_contacts['Y'][contact_type[1]][case].append(y_coords)  # noqa: E501
                        contacts_counter -= len(x_coords)
                    elif custom and cus_inter_res[i] is not None:
                        pair1, pair2 = select_custom_contacts(
                            contacts=cus_inter_res,
                            idx=i,
                            c_type=contact_type[1],
                            max_contacts=contacts_counter,
                            )
                        selected_contacts['X'][contact_type[1]][case].append(pair1)  # noqa: E501
                        selected_contacts['Y'][contact_type[1]][case].append(pair2)  # noqa: E501
                        contacts_counter -= len(pair1)
            
            for i, norm_intra_mtx in enumerate(blended_intra_mtxs):
                contacts_counter = randint(1, max_contacts)
                case = input_seq_keys[i]
                while contacts_counter - 1 > 0:
                    custom = random() < min_contacts_weight
                    if not custom or list(cus_intra_res.values())[i] is None:
                        x_coords, y_coords = select_contacts(
                            coords=norm_intra_mtx,
                            max_num_points=contacts_counter
                            )
                        selected_contacts['X'][contact_type[0]][case].append(x_coords)  # noqa: E501
                        selected_contacts['Y'][contact_type[0]][case].append(y_coords)  # noqa: E501
                        contacts_counter -= len(x_coords)
                    elif custom and list(cus_intra_res.values())[i] is not None:
                        pair1, pair2 = select_custom_contacts(
                            contacts=cus_intra_res,
                            idx=i,
                            c_type=contact_type[0],
                            max_contacts=contacts_counter,
                            )
                        selected_contacts['X'][contact_type[0]][case].append(pair1)  # noqa: E501
                        selected_contacts['Y'][contact_type[0]][case].append(pair2)  # noqa: E501
                        contacts_counter -= len(pair1)
        log.info(S('done'))
    else:
        log.info(T(f'Choosing at most {max_contacts} contacts for every conformer from the blended contact heatmap.'))  # noqa: E501
        for _ in range(nconfs):
            for i, norm_inter_mtx in enumerate(blended_inter_mtxs):
                contacts_counter = randint(1, max_contacts)
                case = combo_chains[i]
                while contacts_counter > 0:
                    x_coords, y_coords = select_contacts(
                        coords=norm_inter_mtx,
                        max_num_points=contacts_counter
                        )
                    selected_contacts['X'][contact_type[1]][case].append(x_coords)  # noqa: E501
                    selected_contacts['Y'][contact_type[1]][case].append(y_coords)  # noqa: E501
                    contacts_counter -= len(x_coords)
            for i, norm_intra_mtx in enumerate(blended_intra_mtxs):
                contacts_counter = randint(1, max_contacts)
                case = input_seq_keys[i]
                while contacts_counter > 0:
                    x_coords, y_coords = select_contacts(
                        coords=norm_intra_mtx,
                        max_num_points=contacts_counter
                        )
                    selected_contacts["X"][contact_type[0]][case].append(x_coords)  # noqa: E501
                    selected_contacts["Y"][contact_type[0]][case].append(y_coords)  # noqa: E501
                    contacts_counter -= len(x_coords)
        log.info(S('done'))


if __name__ == "__main__":
    libcli.maincli(ap, main)
