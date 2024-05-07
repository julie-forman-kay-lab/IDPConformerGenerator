"""
Client for using the expended database to build complexes.

Build from a database of CA contacts, torsion angles,
and secondary structure information. Database is as
created by `idpconfgen contacts` CLI.

USAGE:
    $ idpconfgen complex -db contacts.json -seq sequence.fasta --plot
"""
import argparse
import pickle
import re
from collections import defaultdict
from functools import partial
from itertools import combinations, cycle, product
from multiprocessing import Pool, Queue
from random import randint, random

import numpy as np

from idpconfgen import Path, log
from idpconfgen.cli_build import (
    EnergyLogSaver,
    gen_PDB_from_conformer,
    get_adjacent_angles,
    )
from idpconfgen.components.bgeo_strategies import (
    add_bgeo_strategy_arg,
    bgeo_error_msg,
    bgeo_exact_name,
    bgeo_fixed_name,
    bgeo_int2cart_name,
    bgeo_sampling_name,
    bgeo_strategies,
    bgeo_strategies_default,
    )
from idpconfgen.components.bgeo_strategies.fixed import get_cycle_bend_angles
from idpconfgen.components.energy_threshold_type import add_et_type_arg
from idpconfgen.components.plots.plotfuncs import plot_contacts_matrix
from idpconfgen.components.residue_tolerance import add_res_tolerance_groups
from idpconfgen.components.sidechain_packing import (
    DEFAULT_SDM,
    add_mcsce_subparser,
    add_sidechain_method,
    get_sidechain_packing_parameters,
    sidechain_packing_methods,
    )
from idpconfgen.components.xmer_probs import (
    add_xmer_arg,
    compress_xmer_to_key,
    prepare_xmer_probs,
    )
from idpconfgen.core.build_definitions import (
    backbone_atoms,
    build_bend_CA_C_O,
    build_bend_H_N_C,
    distance_C_O,
    distance_H_N,
    forcefields,
    n_proline_h_coord_at_origin,
    n_terminal_h_coords_at_origin,
    sidechain_templates,
    )
from idpconfgen.core.definitions import dssp_ss_keys
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs import libcli
from idpconfgen.libs.libbuild import (
    create_sidechains_masks_per_residue,
    get_cycle_bond_type,
    get_cycle_distances_backbone,
    init_conflabels,
    init_confmasks,
    prepare_energy_function,
    prepare_slice_dict,
    )
from idpconfgen.libs.libcalc import (
    calc_residue_num_from_index,
    calc_torsion_angles,
    make_coord_Q,
    make_coord_Q_COO,
    make_coord_Q_planar,
    place_sidechain_template,
    rotate_coordinates_Q_njit,
    rrd10_njit,
    )
from idpconfgen.libs.libcomplex import (
    contact_matrix,
    contact_type,
    electropotential_matrix,
    extract_interpairs_from_db,
    extract_intrapairs_from_db,
    get_contact_distances,
    process_custom_contacts,
    select_contacts,
    select_custom_contacts,
    update_distance_distribution_matrix,
    )
from idpconfgen.libs.libfilter import aligndb
from idpconfgen.libs.libhigherlevel import bgeo_reduce
from idpconfgen.libs.libio import (
    change_extension,
    file_exists,
    make_folder_or_cwd,
    read_dictionary_from_disk,
    )
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libparse import (
    get_trimer_seq_njit,
    remap_sequence,
    remove_empty_keys,
    translate_seq_to_3l,
    update_chars_lower,
    )
from idpconfgen.libs.libstructure import Structure, find_sa_residues
from idpconfgen.logger import S, T, init_files, pre_msg, report_on_crash


_file = Path(__file__).myparents()
LOGFILESNAME = '.idpconfgen_complex'

BGEO_full = {}
BGEO_trimer = {}
BGEO_res = {}

INT2CART = None

ANGLES = None
BEND_ANGS = None
BOND_LENS = None
SLICEDICT_XMERS = None
XMERPROBS = None
GET_ADJ = None

MULTICHAIN = False

CONF_NUMBER = Queue()
RANDOMSEEDS = Queue()

ALL_ATOM_LABELS = None
ALL_ATOM_MASKS = None
ALL_ATOM_EFUNC = None
TEMPLATE_LABELS = None
TEMPLATE_MASKS = None
TEMPLATE_EFUNC = None


def are_globals(bgeo_strategy):
    """Assess if global variables needed for building are populated."""
    if bgeo_strategy == bgeo_sampling_name:
        return all((
            ALL_ATOM_LABELS,
            ALL_ATOM_MASKS,
            ALL_ATOM_EFUNC,
            TEMPLATE_LABELS,
            TEMPLATE_MASKS,
            TEMPLATE_EFUNC,
            BGEO_full,
            BGEO_trimer,
            BGEO_res,
            ))

    elif bgeo_strategy in (bgeo_exact_name, bgeo_fixed_name):
        return all((
            ALL_ATOM_LABELS,
            ALL_ATOM_MASKS,
            ALL_ATOM_EFUNC,
            TEMPLATE_LABELS,
            TEMPLATE_MASKS,
            TEMPLATE_EFUNC,
            ))

    elif bgeo_strategy == bgeo_int2cart_name:
        return all((
            ALL_ATOM_LABELS,
            ALL_ATOM_MASKS,
            ALL_ATOM_EFUNC,
            TEMPLATE_LABELS,
            TEMPLATE_MASKS,
            TEMPLATE_EFUNC,
            BGEO_full,
            BGEO_trimer,
            BGEO_res,
            INT2CART,
            ))
    else:
        raise AssertionError(bgeo_error_msg.format(bgeo_strategy))


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

ap.add_argument(
    '--contacts-checkpoint',
    help=(
        'Path to the intermediate binary .PKL database of chain & residue '
        'combinations, intra- & inter- contact frequencies, and respective '
        'intra- & inter- contact distances. For use with the same input '
        'sequence/folded protein. Defaults to contacts_checkpoint.pkl'
        ),
    default="contacts_checkpoint.pkl",
    )

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

ap.add_argument(
    '--ignore-intra',
    help=(
        "Will not build ensembles with intramolecular contacts. "
        "Intramolecular contacts from the database will still be used. "
        "May dramatically increase processing speed."
        ),
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


ENERGYLOGSAVER = EnergyLogSaver()


def main(
        input_seq,
        database,
        contacts_checkpoint="contacts_checkpoint.pkl",
        phos=None,
        folded_structure=None,
        custom_contacts=None,
        ignore_sasa=False,
        ignore_intra=False,
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
    global ANGLES, BEND_ANGS, BOND_LENS, SLICEDICT_XMERS, XMERPROBS, GET_ADJ
    
    if residue_tolerance is not None:
        _restol = str(residue_tolerance)[1:-1]
        log.info(S(f"Building with residue tolerances: {_restol}"))
    
    all_valid_ss_codes = ''.join(dssp_ss_keys.valid)
    # ensuring some parameters do not overlap
    dloop = not dloop_off
    any_def_loops = any((dloop, dhelix, dstrand))
    non_overlapping_parameters = (any_def_loops, dany, duser)
    _sum = sum(map(bool, non_overlapping_parameters))

    if _sum > 1:
        emsg = (
            'Note (dloop, dstrand, dhelix), dany, and duser '
            'are mutually exclusive.'
            )
        raise ValueError(emsg)
    elif _sum < 1:
        raise ValueError("Give at least one sampling option.")

    del _sum
    del non_overlapping_parameters
    # done
    if dany:
        # will sample the database disregarding the SS annotation
        dssp_regexes = [all_valid_ss_codes]

    elif any((dloop, dhelix, dstrand)):
        dssp_regexes = []
        if dloop:
            dssp_regexes.append("L")
        if dhelix:
            dssp_regexes.append("H")
        if dstrand:
            dssp_regexes.append("E")

    elif duser:
        # this is a very advanced option,
        # users should know what they are doing :-)
        dssp_regexes = duser

    else:
        raise AssertionError("One option is missing. Code shouldn't be here.")

    assert isinstance(dssp_regexes, list), \
        f"`dssp_regexes` should be a list at this point: {type(dssp_regexes)}"
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
    
    if bgeo_strategy == bgeo_exact_name:
        try:
            _, ANGLES, BEND_ANGS, BOND_LENS, secondary, primary = aligndb(db, True)  # noqa: E501
        except KeyError:
            log.info(S('!!!!!!!!!!!!!!!'))
            log.info(S(
                'DATABASE ERROR: '
                'the `database` requested is invalid. Please give the '
                'database generated with `bgeodb`. See the usage '
                'documentation for details while using '
                '`--bgeo-strategy exact`.'
                ))
            return
    else:
        _, ANGLES, secondary, primary = aligndb(db)
        
    del db
    
    xmer_probs_tmp = prepare_xmer_probs(xmer_probs)
    
    log.info(S('done'))
    
    log.info(T('Checking validity of --contacts-checkpoint'))
    process_extended_db = True
    if not contacts_checkpoint.endswith(".pkl"):
        log.info(S(f"Incorrect formatting for contacts checkpoint: {contacts_checkpoint}"))  # noqa: E501
    else:
        contacts_checkpoint = change_extension(contacts_checkpoint, "pkl")
        if file_exists(contacts_checkpoint):
            log.info(S(f"Found a checkpoint for contacts information: {contacts_checkpoint}"))  # noqa: E501
            process_extended_db = False
        else:
            log.info(S("No checkpoint for contacts information found, will save a checkpoint after processing"))  # noqa: E501
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
    
    if process_extended_db:
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
        intra_dists = []
        for contact, dist in intra_execute_pool:
            if contact is False:
                continue
            intra_contacts.append(contact)
            intra_dists.append(dist)
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
        inter_dists = []
        for contact, dist in inter_execute_pool:
            if contact is False:
                continue
            inter_contacts.append(contact)
            inter_dists.append(dist)
        log.info(S("done"))

        all_contacts = intra_contacts + inter_contacts
        all_dists = intra_dists + inter_dists
        all_contacts_filtered = [c for c in all_contacts if c]
        all_dists_filtered = [d for d in all_dists if d]  # noqa: F841

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
            combo_chains = list(product(fld_sa_chains, in_chains)) + combo_chains  # noqa: E501
            if phos:
                combo_mod_seqs = list(product(fld_sa_seqs, mod_in_seqs)) + combo_mod_seqs  # noqa: E501

        if len(all_contacts_filtered) == 0:
            log.info("WARNING: No contacts found. Your database is invalid.")
            log.info("Please use `idpconfgen contacts` to make the extended database.")  # noqa: E501
            return
        else:
            all_db_merged = [(c, d) for c, d in zip(all_contacts_filtered, all_dists_filtered)]  # noqa: E501

            log.info(T("Calculating all-contacts probability matrix"))
            inter_c_mtxs = []
            inter_d_mtxs = []
            if len(combo_seqs) >= 1:
                inter_seqs = []
                for combo in combo_seqs:
                    c1 = combo[0]
                    c2 = combo[1]
                    if type(c1) is list:
                        temp_c_mtxs = np.zeros(shape=(1, len(c2)))
                        temp_d_mtxs = np.empty(shape=(1, len(c2)))
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
                                all_db_merged,
                                ncores=ncores,
                                )
                            c_mtx = np.zeros((len(seq), len(c2)))
                            d_mtx = np.empty((len(seq), len(c2)), dtype=object)
                            d_mtx.fill({})
                            for hit, dist in matrix_execute_pool:
                                c_mtx = np.add(c_mtx, hit)
                                d_mtx = update_distance_distribution_matrix(dist, d_mtx)  # noqa: E501
                            temp_seq += seq
                            temp_c_mtxs = np.vstack((c_mtx, temp_c_mtxs))
                            temp_d_mtxs = np.vstack((d_mtx, temp_d_mtxs))
                            log.info(S("done"))
                        # First row will be all zeroes, so remove them
                        c_mtx = temp_c_mtxs[:-1]
                        d_mtx = temp_d_mtxs[:-1]
                        norm_mtx = (c_mtx - np.min(c_mtx)) / (np.max(c_mtx) - np.min(c_mtx))  # noqa: E501
                        inter_c_mtxs.append(norm_mtx)
                        inter_d_mtxs.append(d_mtx)
                        inter_seqs.append([c2, temp_seq])
                    else:
                        log.info(S(f"calculating intermolecular probabilities between {c1} and {c2}"))  # noqa: E501
                        matrix_consume = partial(contact_matrix, sequence=[c1, c2])  # noqa: E501
                        matrix_execute = partial(
                            report_on_crash,
                            matrix_consume,
                            ROC_exception=Exception,
                            ROC_folder=output_folder,
                            ROC_prefix=_name,
                            )
                        matrix_execute_pool = pool_function(
                            matrix_execute,
                            all_db_merged,
                            ncores=ncores,
                            )
                        c_mtx = np.zeros((len(c1), len(c2)))
                        d_mtx = np.empty((len(c1), len(c2)), dtype=object)
                        d_mtx.fill({})
                        for hit, dist in matrix_execute_pool:
                            c_mtx = np.add(c_mtx, hit)
                            d_mtx = update_distance_distribution_matrix(dist, d_mtx)  # noqa: E501
                        norm_mtx = (c_mtx - np.min(c_mtx)) / (np.max(c_mtx) - np.min(c_mtx))  # noqa: E501
                        inter_c_mtxs.append(norm_mtx)
                        inter_d_mtxs.append(d_mtx)
                        inter_seqs.append([c2, c1])
                        log.info(S("done"))

            if ignore_intra:
                log.info("Since --ignore-intra has been turned on. We will not build intramolecular complexes.")  # noqa: E501
            else:
                log.info(T("Calculating intramolecular contacts probability matrix"))  # noqa: E501
                intra_c_mtxs = []
                intra_d_mtxs = []
                # Don't need `intra_seqs` here because `in_seqs` is the same
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
                        all_db_merged,
                        ncores=ncores,
                        )
                    c_mtx = np.zeros((len(seq), len(seq)))
                    d_mtx = np.empty((len(seq), len(seq)), dtype=object)
                    d_mtx.fill({})
                    for hit, dist in matrix_execute_pool:
                        c_mtx = np.add(c_mtx, hit)
                        d_mtx = update_distance_distribution_matrix(dist, d_mtx)
                    norm_mtx = (c_mtx - np.min(c_mtx)) / (np.max(c_mtx) - np.min(c_mtx))  # noqa: E501
                    intra_c_mtxs.append(norm_mtx)
                    intra_d_mtxs.append(d_mtx)

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
        if ignore_intra is False:
            intra_electro_mtx = []
            for seq in in_seqs:
                mtx = electropotential_matrix(seq, ph)
                norm_electro_mtx = (mtx - np.min(mtx)) / (np.max(mtx) - np.min(mtx))  # noqa: E501
                intra_electro_mtx.append(norm_electro_mtx)
        log.info(S('done'))

        log.info(T('saving heatmap distribution plots to output directory'))
        output_plots = make_folder_or_cwd(str(output_folder) + "/plots")
        if not (plot_name.endswith(".png") or plot_name.endswith(".jpg")):
            plot_path = Path(plot_name)
            plot_name = plot_path.with_suffix('.png')

        if len(all_contacts_filtered) > 0:
            if len(inter_c_mtxs) >= 1:
                for i, mtx in enumerate(inter_c_mtxs):
                    plot_out = str(output_plots) + f'/inter_{i}_' + plot_name
                    plot_contacts_matrix(mtx, inter_seqs[i], plot_out)
            if ignore_intra is False:
                for i, mtx in enumerate(intra_c_mtxs):
                    plot_out = str(output_plots) + f'/intra_{i}_' + plot_name
                    plot_contacts_matrix(mtx, in_seqs[i].upper(), plot_out)

        if len(inter_electro_mtx) >= 1:
            for e, mtx in enumerate(inter_electro_mtx):
                plot_electro_out = str(output_plots) + f'/electro_inter_{e}_' + plot_name  # noqa: E501
                plot_contacts_matrix(
                    mtx,
                    inter_seqs[e],
                    plot_electro_out,
                    title="Contacts Frequency Heatmap (Electrostatic)"
                    )
        if ignore_intra is False:
            for e, mtx in enumerate(intra_electro_mtx):
                plot_electro_out = str(output_plots) + f'/electro_intra_{e}_' + plot_name  # noqa: E501
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
            if len(inter_c_mtxs) >= 1:
                for i, mtx in enumerate(inter_c_mtxs):
                    e_mtx = inter_electro_mtx[i]
                    plot_out = str(output_plots) + f'/blend_inter_{i}_' + plot_name  # noqa: E501
                    blended_mtx = (1 - minimized_blend_weight) * mtx + minimized_blend_weight * e_mtx  # noqa: E501
                    norm_blended_mtx = (blended_mtx - np.min(blended_mtx)) / (np.max(blended_mtx) - np.min(blended_mtx))  # noqa: E501
                    blended_inter_mtxs.append(norm_blended_mtx)
                    plot_contacts_matrix(
                        norm_blended_mtx,
                        inter_seqs[i],
                        plot_out,
                        title=f"Inter Contacts Frequency Heatmap (Blended {100 - blend_weight}:{blend_weight})"  # noqa: #501
                        )
            if ignore_intra is False:
                for i, mtx in enumerate(intra_c_mtxs):
                    e_mtx = intra_electro_mtx[i]
                    plot_out = str(output_plots) + f'/blend_intra_{i}_' + plot_name  # noqa: E501
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
    
        # Make a savepoint here (--contacts-checkpoint)
        # Contains all information about chain combinations,
        # residue numbers, contact frequencies, and CA distances
        with open(contacts_checkpoint, 'wb') as f:
            pickle.dump(combo_chains, f)
            pickle.dump(combo_res, f)
            pickle.dump(blended_inter_mtxs, f)
            pickle.dump(blended_intra_mtxs, f)
            pickle.dump(inter_d_mtxs, f)
            pickle.dump(intra_d_mtxs, f)
    else:
        # Extract all necessary information from checkpoint
        with open(contacts_checkpoint, 'rb') as f:
            combo_chains = pickle.load(f)
            combo_res = pickle.load(f)
            blended_inter_mtxs = pickle.load(f)
            blended_intra_mtxs = pickle.load(f)
            inter_d_mtxs = pickle.load(f)
            intra_d_mtxs = pickle.load(f)
    
    # We do not change "combo_chains" since that's our anchor for alignment.
    # If we do not have custom contacts for specific combo chains, store None.
    if custom_contacts is not None:
        cus_inter_res, cus_intra_res = process_custom_contacts(
            custom_contacts,
            combo_chains,
            combo_res,
            input_seq,
            ignore_sasa=ignore_sasa,
            )
    input_seq_keys = list(input_seq.keys())
    
    # Custom and knowledge based contacts are delineated by datatype.
    # Custom contacts are lists while knolwedge based ones are tuple.
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
            if ignore_intra is False:
                for i, norm_intra_mtx in enumerate(blended_intra_mtxs):
                    contacts_counter = randint(1, max_contacts)
                    case = input_seq_keys[i]
                    while contacts_counter - 1 > 0:
                        custom = random() < min_contacts_weight
                        if not custom or list(cus_intra_res.values())[i] is None:  # noqa: E501
                            x_coords, y_coords = select_contacts(
                                coords=norm_intra_mtx,
                                max_num_points=contacts_counter
                                )
                            selected_contacts['X'][contact_type[0]][case].append(x_coords)  # noqa: E501
                            selected_contacts['Y'][contact_type[0]][case].append(y_coords)  # noqa: E501
                            contacts_counter -= len(x_coords)
                        elif custom and list(cus_intra_res.values())[i] is not None:  # noqa: E501
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
            if ignore_intra is False:
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
    
    # TODO extracting distance distributions from database
    # For custom-contacts we would need to rescan the database for
    # residue pairs
    # - Make a d_mtx for every custom contact and align it with
    #   `cus_inter_res` and `cus_intra_res`
    
    # NOTE work with generalizable inter- for IDP-Folded and IDP-IDP before
    # algorithm for intramolecular contacts
    for conf in range(nconfs):
        log.info(T(f"Generating fragments for conformer {conf + 1}"))
        conf_out = make_folder_or_cwd(str(output_folder) + f"/conformer_{conf + 1}")  # noqa: E501
        
        for idx, chains in enumerate(combo_chains):
            res = combo_res[idx]
            inter_mtx = inter_d_mtxs[idx]
            inter_x_coords = selected_contacts["X"][contact_type[1]][chains]
            inter_y_coords = selected_contacts["Y"][contact_type[1]][chains]

            for i, x_coords in enumerate(inter_x_coords):
                xy = []
                y_coords = inter_y_coords[i]
                for j, x in enumerate(x_coords):
                    xy.append((x, y_coords[j]))
                res_combos = []
                distances = []
                for coords in xy:
                    d, r = get_contact_distances(
                        coords,
                        res,
                        inter_mtx,
                        folded=True,
                        )
                    res_combos.append(r)
                    distances.append(d)

            # We want to only build fragments of IDP
            seq1_id = chains[0]
            seq2_id = chains[1]
            # We have two IDPs
            if seq1_id[0] == ">" and seq2_id[0] == ">":
                pass
            # IDP-fld case
            elif seq1_id[0] != ">" and seq2_id[0] == ">":
                idp_sequences = []
                for res_pair in res_combos:
                    idp_seq = ""
                    idp_res = res_pair[1]
                    for r in idp_res:
                        idp_seq += input_seq[seq2_id][r]
                    idp_sequences.append(idp_seq)
            
            SLICEDICT_XMERS = []
            XMERPROBS = []
            GET_ADJ = []
            
            for i, seq in enumerate(idp_sequences):
                log.info(S(f"Preparing database for fragment #{i + 1}: {seq}"))
                
                SLICEDICT_XMERS.append(prepare_slice_dict(
                    primary,
                    seq,
                    dssp_regexes=dssp_regexes,
                    secondary=secondary,
                    mers_size=xmer_probs_tmp.size,
                    res_tolerance=residue_tolerance,
                    ncores=ncores,
                    ))
                
                remove_empty_keys(SLICEDICT_XMERS[i])
                # updates user defined fragment sizes and probabilities to the
                # ones actually observed
                _ = compress_xmer_to_key(xmer_probs_tmp, sorted(SLICEDICT_XMERS[i].keys()))  # noqa: E501
                XMERPROBS.append(_.probs)

                GET_ADJ.append(get_adjacent_angles(
                    sorted(SLICEDICT_XMERS[i].keys()),
                    XMERPROBS[i],
                    seq,
                    ANGLES,
                    bgeo_strategy,
                    SLICEDICT_XMERS[i],
                    csss=None,
                    residue_tolerance=residue_tolerance,
                    ))

                log.info(S("done"))
                
            ENERGYLOGSAVER.start(output_folder.joinpath(energy_log))
            
            sidechain_parameters = \
                get_sidechain_packing_parameters(kwargs, sidechain_method)
            for i, seq in enumerate(idp_sequences):
                consume = partial(
                    _build_conformers,
                    index=i,
                    input_seq=seq,
                    output_folder=conf_out,
                    sidechain_parameters=sidechain_parameters,
                    sidechain_method=sidechain_method,
                    bgeo_strategy=bgeo_strategy,
                    **kwargs,
                    )
                
                execute = partial(
                    report_on_crash,
                    consume,
                    ROC_exception=Exception,
                    ROC_folder=output_folder,
                    ROC_prefix=_name,
                    )
                
                with Pool(ncores) as pool:
                    imap = pool.imap(execute, range(ncores))
                    for _ in imap:
                        pass


def populate_globals(
        *,
        input_seq=None,
        bgeo_strategy=bgeo_strategies_default,
        bgeo_path=None,
        forcefield=None,
        **efunc_kwargs):
    """
    Populate global variables needed for building.
    
    Refer to `cli_build.py` for documentation.
    """
    if not isinstance(input_seq, str):
        raise ValueError(
            '`input_seq` not valid. '
            f'Expected string found {type(input_seq)}'
            )

    if bgeo_strategy not in bgeo_strategies:
        raise AssertionError(bgeo_error_msg.format(bgeo_strategy))

    if bgeo_strategy in (bgeo_sampling_name, bgeo_int2cart_name, bgeo_exact_name):  # noqa: E501
        from idpconfgen.components.bgeo_strategies.sampling import bgeo_sampling_path  # noqa: E501  # isort:skip

        if bgeo_path is None:
            bgeo_path = bgeo_sampling_path

        global BGEO_full, BGEO_trimer, BGEO_res
        BGEO_full.update(read_dictionary_from_disk(bgeo_sampling_path))
        _1, _2 = bgeo_reduce(BGEO_full)
        BGEO_trimer.update(_1)
        BGEO_res.update(_2)
        del _1, _2
        assert BGEO_full
        assert BGEO_trimer
        assert BGEO_res
        # this asserts only the first layer of keys
        assert list(BGEO_full.keys()) == list(BGEO_trimer.keys()) == list(BGEO_res.keys())  # noqa: E501

    # Also prepare BGEO_int2cart when needed
    if bgeo_strategy == bgeo_int2cart_name:
        global INT2CART
        from idpconfgen.components.bgeo_strategies.int2cart.bgeo_int2cart import BGEO_Int2Cart  # noqa: E501  # isort:skip
        try:
            INT2CART = BGEO_Int2Cart()
        except RuntimeError as e:
            log.info(S(
                "WARNING: please use CUDA compatible GPUs while running"
                "--bgeo_strategy int2cart."
                ))
            log.info(S(f"Error: {e}"))

    # populates the labels
    global ALL_ATOM_LABELS, ALL_ATOM_MASKS, ALL_ATOM_EFUNC
    global TEMPLATE_LABELS, TEMPLATE_MASKS, TEMPLATE_EFUNC

    topobj = forcefield(add_OXT=True, add_Nterminal_H=True)

    ALL_ATOM_LABELS = init_conflabels(input_seq, topobj.atom_names)
    TEMPLATE_LABELS = init_conflabels(remap_sequence(input_seq), topobj.atom_names)  # noqa: E501

    ALL_ATOM_MASKS = init_confmasks(ALL_ATOM_LABELS.atom_labels)
    TEMPLATE_MASKS = init_confmasks(TEMPLATE_LABELS.atom_labels)

    ALL_ATOM_EFUNC = prepare_energy_function(
        ALL_ATOM_LABELS.atom_labels,
        ALL_ATOM_LABELS.res_nums,
        ALL_ATOM_LABELS.res_labels,
        topobj,
        **efunc_kwargs)

    TEMPLATE_EFUNC = prepare_energy_function(
        TEMPLATE_LABELS.atom_labels,
        TEMPLATE_LABELS.res_nums,
        TEMPLATE_LABELS.res_labels,
        topobj,
        **efunc_kwargs)

    del topobj
    return


# private function because it depends on the global `CONF_NUMBER`
# which is assembled in `main()`
def _build_conformers(
        *args,
        index=None,
        input_seq=None,
        conformer_name='fragment',
        output_folder=None,
        sidechain_parameters=None,
        bgeo_strategy=bgeo_strategies_default,
        **kwargs,
        ):
    """Arrange building of conformers and saves them to PDB files."""
    ROUND = np.round

    # TODO: this has to be parametrized for the different HIS types
    input_seq_3_letters = translate_seq_to_3l(input_seq)
  
    builder = conformer_generator(
        index=index,
        input_seq=input_seq,
        random_seed=RANDOMSEEDS.get(),
        sidechain_parameters=sidechain_parameters,
        bgeo_strategy=bgeo_strategy,
        **kwargs)

    atom_labels, residue_numbers, _residue_labels = next(builder)

    energy, coords = next(builder)
    
    pdb_string = gen_PDB_from_conformer(
        input_seq_3_letters,
        atom_labels,
        residue_numbers,
        ROUND(coords, decimals=3),
        )
    
    fname = f'{conformer_name}_{index + 1}.pdb'
    
    with open(Path(output_folder, fname), 'w') as fout:
        fout.write(pdb_string)
    ENERGYLOGSAVER.save(fname, energy)

    del builder
    return


def conformer_generator(
        *,
        index=None,
        input_seq=None,
        generative_function=None,
        disable_sidechains=True,
        sidechain_method='faspr',
        energy_threshold_backbone=10,
        energy_threshold_sidechains=1000,
        bgeo_strategy=bgeo_strategies_default,
        bgeo_path=None,
        forcefield=None,
        random_seed=0,
        sidechain_parameters=None,
        **energy_funcs_kwargs,
        ):
    """
    Build conformers.

    Refer to documentation in `cli_build.py`.
    """
    if not isinstance(input_seq, str):
        raise ValueError(f'`input_seq` must be given! {input_seq}')
    if sidechain_method not in sidechain_packing_methods:
        raise ValueError(
            f'{sidechain_method} not in `sidechain_packing_methods`. '
            f'Expected {list(sidechain_packing_methods.keys())}.'
            )

    log.info(f'random seed: {random_seed}')
    np.random.seed(random_seed)
    seed_report = pre_msg(f'seed {random_seed}', sep=' - ')

    # prepares protein sequences
    all_atom_input_seq = input_seq
    template_input_seq = remap_sequence(all_atom_input_seq)
    template_seq_3l = translate_seq_to_3l(template_input_seq)

    ANY = np.any
    BUILD_BEND_H_N_C = build_bend_H_N_C
    CALC_TORSION_ANGLES = calc_torsion_angles
    DISTANCE_NH = distance_H_N
    DISTANCE_C_O = distance_C_O
    ISNAN = np.isnan
    GET_TRIMER_SEQ = get_trimer_seq_njit
    MAKE_COORD_Q_COO_LOCAL = make_coord_Q_COO
    MAKE_COORD_Q_PLANAR = make_coord_Q_planar
    MAKE_COORD_Q_LOCAL = make_coord_Q
    NAN = np.nan
    NORM = np.linalg.norm
    # the N terminal Hs are three for all atoms but only two for Proline
    # depending whether the first residue is a Proline, we use one template
    # or another.
    N_TERMINAL_H = n_proline_h_coord_at_origin if input_seq[0] == "P" else n_terminal_h_coords_at_origin  # noqa: E501
    PI2 = np.pi * 2
    PLACE_SIDECHAIN_TEMPLATE = place_sidechain_template
    RAD_60 = np.radians(60)
    RC = np.random.choice
    RINT = randint
    ROT_COORDINATES = rotate_coordinates_Q_njit
    RRD10 = rrd10_njit
    SIDECHAIN_TEMPLATES = sidechain_templates
    SUM = np.nansum
    global BGEO_full
    global BGEO_trimer
    global BGEO_res
    global ALL_ATOM_LABELS
    global ALL_ATOM_MASKS
    global ALL_ATOM_EFUNC
    global TEMPLATE_LABELS
    global TEMPLATE_MASKS
    global TEMPLATE_EFUNC
    global XMERPROBS
    global SLICEDICT_MONOMERS
    global SLICEDICT_XMERS
    global GET_ADJ

    del input_seq

    # these flags exist to populate the global variables in case they were not
    # populated yet. Global variables are populated through the main() function
    # if the script runs as CLI. Otherwise, if conformer_generator() is imported
    # and used directly, the global variables need to be configured here.
    if not are_globals(bgeo_strategy):
        if forcefield not in forcefields:
            raise ValueError(
                f'{forcefield} not in `forcefields`. '
                f'Expected {list(forcefields.keys())}.'
                )
        populate_globals(
            input_seq=all_atom_input_seq,
            bgeo_strategy=bgeo_strategy,
            bgeo_path=bgeo_path,
            forcefield=forcefields[forcefield],
            **energy_funcs_kwargs,
            )

    # semantic exchange for speed al readibility
    with_sidechains = not disable_sidechains

    if with_sidechains:
        log.info(S(f"configuring sidechain method: {sidechain_method}"))
        # we use named arguments here to allow ignored non needed parameters
        # with **kwargs
        build_sidechains = sidechain_packing_methods[sidechain_method](
            input_seq=all_atom_input_seq,
            template_masks=TEMPLATE_MASKS,
            all_atom_masks=ALL_ATOM_MASKS,
            user_parameters=sidechain_parameters,
            )

    # tests generative function complies with implementation requirements
    if generative_function:
        try:
            generative_function(nres=1, cres=0)
        except Exception as err:  # this is generic Exception on purpose
            errmsg = (
                'The `generative_function` provided is not compatible with '
                'the building process. Please read `build_conformers` docstring'
                ' for more details.'
                )
            raise IDPConfGenException(errmsg) from err

    # yields atom labels
    # all conformers generated will share these labels
    yield (
        ALL_ATOM_LABELS.atom_labels,
        ALL_ATOM_LABELS.res_nums,
        ALL_ATOM_LABELS.res_labels,
        )
    all_atom_num_atoms = len(ALL_ATOM_LABELS.atom_labels)
    template_num_atoms = len(TEMPLATE_LABELS.atom_labels)

    all_atom_coords = np.full((all_atom_num_atoms, 3), NAN, dtype=np.float64)
    template_coords = np.full((template_num_atoms, 3), NAN, dtype=np.float64)

    # +2 because of the dummy coordinates required to start building.
    # see later adding dummy coordinates to the structure seed
    bb = np.full((TEMPLATE_MASKS.bb3.size + 2, 3), NAN, dtype=np.float64)
    bb_real = bb[2:, :]  # backbone coordinates without the dummies

    # coordinates for the carbonyl oxygen atoms
    bb_CO = np.full((TEMPLATE_MASKS.COs.size, 3), NAN, dtype=np.float64)

    # notice that NHydrogen_mask does not see Prolines
    bb_NH = np.full((TEMPLATE_MASKS.NHs.size, 3), NAN, dtype=np.float64)
    bb_NH_idx = np.arange(len(bb_NH))
    # Creates masks and indexes for the `for` loop used to place NHs.
    # The first residue has no NH, prolines have no NH.
    non_pro = np.array(list(template_input_seq)[1:]) != 'P'
    # NHs index numbers in bb_real
    bb_NH_nums = np.arange(3, (len(template_input_seq) - 1) * 3 + 1, 3)[non_pro]  # noqa: E501
    bb_NH_nums_p1 = bb_NH_nums + 1
    assert bb_NH.shape[0] == bb_NH_nums.size == bb_NH_idx.size

    # sidechain masks
    # this is sidechain agnostic, works for every sidechain, yet here we
    # use only ALA, PRO, GLY - Mon Feb 15 17:29:20 2021
    ss_masks = create_sidechains_masks_per_residue(
        TEMPLATE_LABELS.res_nums,
        TEMPLATE_LABELS.atom_labels,
        backbone_atoms,
        )
    # ?

    # /
    # creates seed coordinates:
    # because the first torsion angle of a residue is the omega, we need
    # to prepare 2 dummy atoms to simulate the residue -1, so that the
    # first omega can be placed. There is no need to setup specific
    # positions, just to create a place upon which the build atom
    # routine can create a new atom from a torsion.
    dummy_CA_m1_coord = np.array((0.0, 1.0, 1.0))
    dummy_C_m1_coord = np.array((0.0, 1.0, 0.0))
    n_terminal_N_coord = np.array((0.0, 0.0, 0.0))

    # seed coordinates array
    seed_coords = np.array((
        dummy_CA_m1_coord,
        dummy_C_m1_coord,
        n_terminal_N_coord,
        ))
    # ?

    # /
    # prepares method binding
    bbi0_register = []
    bbi0_R_APPEND = bbi0_register.append
    bbi0_R_POP = bbi0_register.pop
    bbi0_R_CLEAR = bbi0_register.clear

    COi0_register = []
    COi0_R_APPEND = COi0_register.append
    COi0_R_POP = COi0_register.pop
    COi0_R_CLEAR = COi0_register.clear

    res_R = []  # residue number register
    res_R_APPEND = res_R.append
    res_R_POP = res_R.pop
    res_R_CLEAR = res_R.clear
    # ?

    # /
    # required inits
    broke_on_start_attempt = False
    start_attempts = 0
    max_start_attempts = 500  # maximum attempts to start a conformer
    # because we are building from a experimental database there can be
    # some angle combinations that fail on our validation process from start
    # if this happens more than `max_start_attemps` the production is canceled.
    # ?

    # /
    # STARTS BUILDING
    conf_n = 1
    while 1:
        # prepares cycles for building process
        bond_lens = get_cycle_distances_backbone()
        bond_type = get_cycle_bond_type()

        if bgeo_strategy == bgeo_fixed_name:
            bend_angles = get_cycle_bend_angles()

        # in the first run of the loop this is unnecessary, but is better to
        # just do it once than flag it the whole time
        template_coords[:, :] = NAN
        bb[:, :] = NAN
        bb_CO[:, :] = NAN
        bb_NH[:, :] = NAN
        for _mask, _coords in ss_masks:
            _coords[:, :] = NAN

        bb[:3, :] = seed_coords  # this contains a dummy coord at position 0

        # add N-terminal hydrogens to the origin

        bbi = 1  # starts at 1 because there are two dummy atoms
        bbi0_R_CLEAR()
        bbi0_R_APPEND(bbi)

        COi = 0  # carbonyl atoms
        COi0_R_CLEAR()
        COi0_R_APPEND(COi)

        # residue integer number
        current_res_number = 0
        res_R_CLEAR()
        res_R_APPEND(current_res_number)

        backbone_done = False
        number_of_trials = 0
        # TODO: use or not to use number_of_trials2? To evaluate in future.
        number_of_trials2 = 0
        number_of_trials3 = 0

        # used only if bgeo_strategy == int2cart
        torsion_records = []

        # run this loop until a specific BREAK is triggered
        while 1:  # 1 is faster than True :-)

            # I decided to use an if-statement here instead of polymorph
            # the else clause to a `generative_function` variable because
            # the resulting overhead from the extra function call and
            # **kwargs handling was greater then the if-statement processing
            # https://pythonicthoughtssnippets.github.io/2020/10/21/PTS14-quick-in-if-vs-polymorphism.html  # noqa: E501
            if generative_function:
                primer_template, agls = generative_function(
                    nres=RINT(1, 6),
                    cres=calc_residue_num_from_index(bbi)
                    )

            else:
                # algorithm for adjacent building
                # TODO
                # primer_template here is used temporarily, and needs to be
                # removed when get_adj becomes an option
                if index is not None:
                    get_adj = GET_ADJ[index]
                else:
                    get_adj = GET_ADJ
                if bgeo_strategy == bgeo_exact_name:
                    primer_template, agls, bangs, blens = get_adj(bbi - 1)  # noqa: E501
                else:
                    primer_template, agls = get_adj(bbi - 1)
                
            # index at the start of the current cycle
            PRIMER = cycle(primer_template)

            try:
                for (omg, phi, psi) in zip(agls[0::3], agls[1::3], agls[2::3]):

                    current_res_number = calc_residue_num_from_index(bbi - 1)

                    # assert the residue being built is of the same nature as
                    # the one in the angles
                    # TODO: remove this assert
                    n_ = next(PRIMER)
                    assert all_atom_input_seq[current_res_number] == n_, \
                        (all_atom_input_seq[current_res_number], n_)

                    curr_res, tpair = GET_TRIMER_SEQ(
                        all_atom_input_seq,
                        current_res_number,
                        )
                    torpair = f'{RRD10(phi)},{RRD10(psi)}'

                    if bgeo_strategy == bgeo_int2cart_name:
                        torsion_records.append((omg, phi, psi))
                        seq = all_atom_input_seq[:current_res_number + 1]

                        tors = np.array(torsion_records)  # omega, phi, psi

                        # phi, psi, omega
                        tors = np.hstack([tors[:, 1:], tors[:, :1]])

                        _ = INT2CART.get_internal_coords(seq, tors)
                        d1, d2, d3, theta1, theta2, theta3 = _

                        bend_angles = [theta3, theta1, theta2]
                        bond_lens = [d1, d2, d3]

                    for torsion_idx, torsion_angle in enumerate((omg, phi, psi)):  # noqa: E501

                        if bgeo_strategy == bgeo_int2cart_name:
                            # needed for correctly calculating Q
                            _bend_angle = (np.pi - bend_angles[torsion_idx]) / 2
                            _bond_lens = bond_lens[torsion_idx]

                        elif bgeo_strategy == bgeo_exact_name:
                            _bend_angle = bangs[torsion_idx]
                            _bond_lens = blens[torsion_idx]

                        elif bgeo_strategy == bgeo_sampling_name:
                            _bt = next(bond_type)

                            try:
                                _bend_angle = RC(BGEO_full[_bt][curr_res][tpair][torpair])  # noqa: E501
                            except KeyError:
                                try:
                                    _bend_angle = RC(BGEO_trimer[_bt][curr_res][tpair])  # noqa: E501
                                except KeyError:
                                    _bend_angle = RC(BGEO_res[_bt][curr_res])

                            _bond_lens = next(bond_lens)[curr_res]

                        elif bgeo_strategy == bgeo_fixed_name:
                            _bend_angle = next(bend_angles)[curr_res]
                            _bond_lens = next(bond_lens)[curr_res]

                        bb_real[bbi, :] = MAKE_COORD_Q_LOCAL(
                            bb[bbi - 1, :],
                            bb[bbi, :],
                            bb[bbi + 1, :],
                            _bond_lens,
                            _bend_angle,
                            torsion_angle,
                            )
                        bbi += 1

                    if bgeo_strategy in (bgeo_int2cart_name, bgeo_sampling_name):  # noqa: E501

                        try:
                            co_bend = RC(BGEO_full['Ca_C_O'][curr_res][tpair][torpair])  # noqa: E501
                        except KeyError:
                            try:
                                co_bend = RC(BGEO_trimer['Ca_C_O'][curr_res][tpair])  # noqa: E501
                            except KeyError:
                                co_bend = RC(BGEO_res['Ca_C_O'][curr_res])

                    elif bgeo_strategy == bgeo_fixed_name:
                        co_bend = build_bend_CA_C_O

                    else:
                        co_bend = bangs[3]
                        DISTANCE_C_O = blens[3]

                    bb_CO[COi, :] = MAKE_COORD_Q_PLANAR(
                        bb_real[bbi - 3, :],
                        bb_real[bbi - 2, :],
                        bb_real[bbi - 1, :],
                        distance=DISTANCE_C_O,
                        bend=co_bend
                        )
                    COi += 1

            except IndexError:
                # IndexError happens when the backbone is complete
                # in this protocol the last atom build was a carbonyl C
                # bbi is the last index of bb + 1, and the last index of
                # bb_real + 2

                # activate flag to finish loop at the end
                backbone_done = True

                # add the carboxyls
                template_coords[TEMPLATE_MASKS.cterm] = \
                    MAKE_COORD_Q_COO_LOCAL(bb[-2, :], bb[-1, :])

            # Adds N-H Hydrogens
            # Not a perfect loop. It repeats for Hs already placed.
            # However, was a simpler solution than matching the indexes
            # and the time cost is not a bottle neck.
            _ = ~ISNAN(bb_real[bb_NH_nums_p1, 0])
            for k, j in zip(bb_NH_nums[_], bb_NH_idx[_]):

                bb_NH[j, :] = MAKE_COORD_Q_PLANAR(
                    bb_real[k + 1, :],
                    bb_real[k, :],
                    bb_real[k - 1, :],
                    distance=DISTANCE_NH,
                    bend=BUILD_BEND_H_N_C,
                    )

            # Adds sidechain template structures
            for res_i in range(res_R[-1], current_res_number + 1):  # noqa: E501

                _sstemplate, _sidechain_idxs = \
                    SIDECHAIN_TEMPLATES[template_seq_3l[res_i]]

                sscoords = PLACE_SIDECHAIN_TEMPLATE(
                    bb_real[res_i * 3:res_i * 3 + 3, :],  # from N to C
                    _sstemplate,
                    )

                ss_masks[res_i][1][:, :] = sscoords[_sidechain_idxs]

            # Transfers coords to the main coord array
            for _smask, _sidecoords in ss_masks[:current_res_number + 1]:
                template_coords[_smask] = _sidecoords

            # / Place coordinates for energy calculation
            #
            # use `bb_real` to do not consider the initial dummy atom
            template_coords[TEMPLATE_MASKS.bb3] = bb_real
            template_coords[TEMPLATE_MASKS.COs] = bb_CO
            template_coords[TEMPLATE_MASKS.NHs] = bb_NH

            if len(bbi0_register) == 1:
                # places the N-terminal Hs only if it is the first
                # fragment being built
                _ = PLACE_SIDECHAIN_TEMPLATE(bb_real[0:3, :], N_TERMINAL_H)
                template_coords[TEMPLATE_MASKS.Hterm, :] = _[3:, :]
                current_Hterm_coords = _[3:, :]
                del _

                # rotating the N-term H's is not needed for G and P
                if template_input_seq[0] not in ('G', 'P'):
                    # rotates only if the first residue is not an
                    # alanie

                    # measure torsion angle reference H1 - HA
                    _h1_ha_angle = CALC_TORSION_ANGLES(
                        template_coords[TEMPLATE_MASKS.H2_N_CA_CB, :]
                        )[0]

                    # given any angle calculated along an axis, calculate how
                    # much to rotate along that axis to place the
                    # angle at 60 degrees
                    _rot_angle = _h1_ha_angle % PI2 - RAD_60

                    current_Hterm_coords = ROT_COORDINATES(
                        template_coords[TEMPLATE_MASKS.Hterm, :],
                        template_coords[1] / NORM(template_coords[1]),
                        _rot_angle,
                        )

                    template_coords[TEMPLATE_MASKS.Hterm, :] = current_Hterm_coords  # noqa: E501
            # ?

            total_energy = TEMPLATE_EFUNC(template_coords)

            if ANY(total_energy > energy_threshold_backbone):
                # reset coordinates to the original value
                # before the last fragment added

                # reset the same fragment maximum 5 times,
                # after that reset also the fragment before
                try:
                    if number_of_trials > 30:
                        bbi0_R_POP()
                        COi0_R_POP()
                        res_R_POP()
                        number_of_trials = 0
                        number_of_trials2 += 1

                    if number_of_trials2 > 5:
                        bbi0_R_POP()
                        COi0_R_POP()
                        res_R_POP()
                        number_of_trials2 = 0
                        number_of_trials3 += 1

                    if number_of_trials3 > 5:
                        bbi0_R_POP()
                        COi0_R_POP()
                        res_R_POP()
                        number_of_trials3 = 0

                    _bbi0 = bbi0_register[-1]
                    _COi0 = COi0_register[-1]
                    _resi0 = res_R[-1]
                except IndexError:
                    # if this point is reached,
                    # we erased until the beginning of the conformer
                    # discard conformer, something went really wrong
                    broke_on_start_attempt = True
                    break  # conformer while loop, starts conformer from scratch

                # clean previously built protein fragment
                bb_real[_bbi0:bbi, :] = NAN
                bb_CO[_COi0:COi, :] = NAN

                # reset also indexes
                bbi = _bbi0
                COi = _COi0
                current_res_number = _resi0

                # remove torsion angle records for this chunk
                if bgeo_strategy == bgeo_int2cart_name:
                    torsion_records = torsion_records[:current_res_number + 1]

                # coords needs to be reset because size of protein next
                # fragments may not be equal
                template_coords[:, :] = NAN
                template_coords[TEMPLATE_MASKS.Hterm, :] = current_Hterm_coords

                # prepares cycles for building process
                # this is required because the last fragment created may have
                # been the final part of the conformer
                if backbone_done:
                    bond_lens = get_cycle_distances_backbone()
                    bond_type = get_cycle_bond_type()

                # we do not know if the next fragment will finish the protein
                # or not
                backbone_done = False
                number_of_trials += 1
                continue  # send back to the fragment while loop

            # if the conformer is valid
            number_of_trials = 0
            bbi0_R_APPEND(bbi)
            COi0_R_APPEND(COi)
            # the residue where the build process stopped
            res_R_APPEND(current_res_number)

            if backbone_done:
                # this point guarantees all protein atoms are built
                break  # fragment while loop
        # END of fragment while loop, go up and build the next fragment

        if broke_on_start_attempt:
            start_attempts += 1
            if start_attempts > max_start_attempts:
                log.error(
                    'Reached maximum amount of re-starts. Canceling... '
                    f'Built a total of {conf_n} conformers.'
                    )
                return
            broke_on_start_attempt = False
            continue  # send back to the fragment while loop

        # we do not want sidechains at this point
        all_atom_coords[ALL_ATOM_MASKS.bb4] = template_coords[TEMPLATE_MASKS.bb4]  # noqa: E501
        all_atom_coords[ALL_ATOM_MASKS.NHs] = template_coords[TEMPLATE_MASKS.NHs]  # noqa: E501
        all_atom_coords[ALL_ATOM_MASKS.Hterm] = template_coords[TEMPLATE_MASKS.Hterm]  # noqa: E501
        all_atom_coords[ALL_ATOM_MASKS.cterm, :] = template_coords[TEMPLATE_MASKS.cterm, :]  # noqa: E501

        if with_sidechains:

            # this is uniformed API for all build_sidechains
            _mask, _new_sd_coords = build_sidechains(template_coords)

            if _new_sd_coords is None:
                _emsg = (
                    "Could not find a solution for sidechains, "
                    "discarding the conformer...")
                log.info(seed_report(_emsg))
                continue

            all_atom_coords[_mask] = _new_sd_coords

            if ALL_ATOM_EFUNC is None:
                total_energy = 0
            else:
                total_energy = ALL_ATOM_EFUNC(all_atom_coords)

            if ANY(total_energy > energy_threshold_sidechains):
                _msg = (
                    'Conformer with energy higher than allowed threshold '
                    '- discarded.'
                    )
                log.info(seed_report(_msg))
                continue

        _total_energy = np.nansum(total_energy)
        _msg = f'finished conf: {conf_n} with energy {_total_energy}'
        log.info(seed_report(_msg))

        yield SUM(total_energy), all_atom_coords
        conf_n += 1


if __name__ == "__main__":
    libcli.maincli(ap, main)
