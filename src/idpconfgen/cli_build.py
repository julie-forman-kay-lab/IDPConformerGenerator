"""
Builds IDP conformers.

Build from a database of torsion angles and secondary structure
information. Database is as created by `idpconfgen torsions` CLI.

USAGE:
    $ idpconfgen build -db torsions.json -seq MMMMMMM...

"""
import argparse
import os
import re
from functools import partial
from itertools import cycle
from multiprocessing import Pool, Queue
from random import randint
from time import time

import numpy as np

from idpconfgen import Path, log
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
from idpconfgen.ldrs_helper import (
    align_coords,
    count_clashes,
    disorder_cases,
    psurgeon,
    )
from idpconfgen.libs import libcli
from idpconfgen.libs.libbuild import (
    build_regex_substitutions,
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
    make_seq_probabilities,
    place_sidechain_template,
    rotate_coordinates_Q_njit,
    rrd10_njit,
    )
from idpconfgen.libs.libfilter import aligndb
from idpconfgen.libs.libhigherlevel import bgeo_reduce
from idpconfgen.libs.libio import (
    make_folder_or_cwd,
    read_dict_from_json,
    read_dictionary_from_disk,
    )
from idpconfgen.libs.libparse import (
    fill_list,
    get_seq_chunk_njit,
    get_trimer_seq_njit,
    remap_sequence,
    remove_empty_keys,
    split_by_ranges,
    split_into_chunks,
    translate_seq_to_3l,
    )
from idpconfgen.libs.libpdb import atom_line_formatter
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resSeq,
    cols_coords,
    parse_pdb_to_array,
    structure_to_pdb,
    )
from idpconfgen.logger import S, T, init_files, pre_msg, report_on_crash


_file = Path(__file__).myparents()
LOGFILESNAME = '.idpconfgen_build'

# Global variables needed to build conformers.
# Why are global variables needed?
# I use global variables to facilitate distributing conformer creation
# processes across multiple cores. In this way cores can read global variables
# fast and with non-significant overhead.

# Bond Geometry library variables
# if __name__ == '__main__', these will be populated in main()
# else will be populated in conformer_generator
# populate_globals() populates these variables once called.
# sampling globals
BGEO_full = {}
BGEO_trimer = {}
BGEO_res = {}
# int2cart globals
INT2CART = None

# ANGLES will be populated in main() with the torsion angles.
# it is not expected SLICES or ANGLES to be populated anywhere else.
# The slice objects from where the builder will feed to extract torsion
# fragments from ANGLES.
ANGLES = None
BEND_ANGS = None
BOND_LENS = None
SLICEDICT_XMERS = None
XMERPROBS = None
GET_ADJ = None

# keeps a record of the conformer numbers written to disk across the different
# cores
CONF_NUMBER = Queue()
RANDOMSEEDS = Queue()

# The conformer building process needs data structures for two different
# identities: the all-atom representation of the input sequence, and the
# corresponding Ala/Gly/Pro template uppon which the coordinates will be built.
# These variables are defined at the module level so they serve as global
# variables to be read by the different process during multiprocessing. Reading
# from global variables is performant in Python multiprocessing. This is the
# same strategy as applied for SLICES and ANGLES.
ALL_ATOM_LABELS = None
ALL_ATOM_MASKS = None
ALL_ATOM_EFUNC = None
TEMPLATE_LABELS = None
TEMPLATE_MASKS = None
TEMPLATE_EFUNC = None

# Global variables for enabling the "long" feature for building
# extended IDPs (e.g. > 300 AA)
GET_ADJ_LONG = None
LONG_FRAGMENTS = None


class _BuildPreparation:
    pass


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


# CLI argument parser parameters
_name = 'build'
_help = 'Builds conformers from database.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)


ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

libcli.add_argument_idb(ap)
libcli.add_argument_seq(ap)

ap.add_argument(
    '-nc',
    '--nconfs',
    help='Number of conformers to build.',
    default=1,
    type=int,
    )

ap.add_argument(
    '--long',
    help=(
        'Switch to enable building long IDPs. '
        'Note this will NOT automatically enable if you have IDPs '
        'longer than 300 AA but it is recommended to turn this on. '
        'Defaults to True.'
        ),
    action="store_true",
    )

ap.add_argument(
    '--long-ranges',
    help=(
        "Custom ranges of residues to build fragmentally for a long IDP "
        "is denoted by dashes for residue numbers and commas for different "
        "ranges. Note that ALL patterns MUST end "
        "at the last residue with a comma. "
        "For ex. a 301 AA long IDP: --long-ranges 1-113,114-210,211-301,"
        "Optional flag. If left empty, generate IDP with fragments of length "
        "150 AA at a time. E.g. the same as --long-ranges 1-150,151-300,"
        ),
    nargs='?',
    )


#########################################
libcli.add_argument_dloopoff(ap)
libcli.add_argument_dhelix(ap)
libcli.add_argument_dstrand(ap)
libcli.add_argument_dany(ap)
libcli.add_argument_duser(ap)
#########################################

ap.add_argument(
    '-csss',
    '--custom-sampling',
    help=(
        'Input .JSON file for probabilistic CSSS. '
        'Will use DSSP codes in this .JSON instead of --dhelix, --dstrand, '
        '--dany. Requires --dloop-off. CSSS.JSON file is as created by the '
        '`idpconfgen csssconv` or `idpconfgen makecsss` command.'
        ),
    default=None,
    )

ap.add_argument(
    '-dsd',
    '--disable-sidechains',
    help='Whether or not to compute sidechains. Defaults to True.',
    action='store_true',
    )

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
        'The energy threshold above which fragments will be rejected '
        'when building the BACKBONE atoms. Defaults to 100 kJ.'
        ),
    default=100.0,
    type=float,
    )

ap.add_argument(
    '-etss',
    '--energy-threshold-sidechains',
    help=(
        'The energy threshold above which conformers will be rejected '
        'after packing the sidechains (ignored if `-dsd`). '
        'Defaults to 250 kJ.'
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
libcli.add_argument_random_seed(ap)
libcli.add_argument_ncores(ap)


class EnergyLogSaver:
    """
    Save conformer energies to a log file.

    This object is intended to be used by the client only.
    It can accommodate calls from different processors, but it is not
    sent to the different processors, it is managed by the main()
    function.
    """

    def start(self, path):
        """Open file for writing."""
        self.dest = open(path, 'w')

    def save(self, confname, energy):
        """Save conformer name and energy to file."""
        self.dest.write(f'{confname},{energy}\n')
        self.dest.flush()

    def close(self):
        """Close file."""
        self.dest.close()


ENERGYLOGSAVER = EnergyLogSaver()


def parse_CSSS(path2csss):
    """
    Prepare CSSS.JSON dictionary for the conformer building process.

    The secondary structure keys are identified.
    The probabilities for each residue are normalized to 1, that is:
    (1 2 2) results in (0.2 0.4 0.4).

    Parameters
    ----------
    path2csss : string
        Path to where the csss_[ID].json file is containing ss_regexes and
        their respective probabilities.

    Returns
    -------
    dict
        First key layer indicats residue number position, second key layer
        indicates the DSSP regex to search for and the values are the
        probabilities.

    set
        A set with all the different secondary structure keys identified in the
        CSSS.JSON file.
    """
    # this function was originally done by @menoliu
    # @joaomcteixeira gave it a touch
    csss_dict = read_dict_from_json(path2csss)
    all_dssps = set()

    # we can use this implementation because dictionaries are sorted by default
    for _resid, dssps in csss_dict.items():
        probabilities = list(dssps.values())
        all_dssps.update(dssps.keys())
        prob_normalized = make_seq_probabilities(probabilities)
        for dssp_code, prob_n in zip(dssps.keys(), prob_normalized):
            dssps[dssp_code] = prob_n

    return csss_dict, all_dssps


def main(
        input_seq,
        database,
        custom_sampling,
        long=False,
        long_ranges=None,
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
        random_seed=0,
        xmer_probs=None,
        output_folder=None,
        energy_log='energies.log',
        sidechain_method=DEFAULT_SDM,
        **kwargs,  # other kwargs target energy function, for example.
        ):
    """
    Execute main client logic.

    Distributes over processors.
    """
    # ensuring some parameters do not overlap
    dloop = not dloop_off
    any_def_loops = any((dloop, dhelix, dstrand))
    non_overlapping_parameters = (any_def_loops, dany, duser, bool(custom_sampling))  # noqa: E501
    _sum = sum(map(bool, non_overlapping_parameters))

    if _sum > 1:
        emsg = (
            'Note (dloop, dstrand, dhelix), dany, duser, and '
            'custom_sampling are mutually exclusive.'
            )
        raise ValueError(emsg)
    elif _sum < 1:
        raise ValueError("Give at least one sampling option.")

    del _sum
    del non_overlapping_parameters
    # done

    output_folder = make_folder_or_cwd(output_folder)
    init_files(log, Path(output_folder, LOGFILESNAME))
    log.info(T('starting the building process'))
    # We only are interested in the first (one) sequence
    if type(input_seq) is dict:
        input_seq = list(input_seq.values())[0]
    log.info(S(f'input sequence: {input_seq}'))
    
    if len(input_seq) > 300:
        if long is False:
            log.info(
                "TIP: if your IDP is longer than ~300 residues, consider "
                "enabling the `--long` flag for faster generation."
                )
        else:
            global LONG_FRAGMENTS
            if long_ranges:
                range_regex = re.compile(r"\d+-\d+,")
                if range_regex.match(long_ranges):
                    ranges = long_ranges.split(',')
                    ranges.pop()  # last element should be empty
                    idx_ranges = []
                    for r in ranges:
                        parts = r.split("-")
                        if len(parts) >= 2:
                            idx_ranges.append(int(parts[1]))
                    long_fragments = split_by_ranges(input_seq, idx_ranges)
                else:
                    log.info(S('Incorrect pattern input. Resorting to default.'))  # noqa: E501
                    log.info(S('Sample pattern is as follows: 1-89,90-191,'))
            else:
                long_fragments = split_into_chunks(input_seq)
            
            for i in range(len(long_fragments) - 1):
                j = i + 1
                # Add overlapping residues to subsequent fragments
                long_fragments[j] = long_fragments[i][-2:] + long_fragments[j]
                
            LONG_FRAGMENTS = long_fragments
    else:
        long = False
        long_ranges = None
    
    # Calculates how many conformers are built per core
    if nconfs < ncores:
        ncores = 1
        conformers_per_core = nconfs
        remaining_confs = 0
    else:
        conformers_per_core = nconfs // ncores
        # in case nconfs is not multiple of ncores, builds the remaining confs
        # at the end
        remaining_confs = nconfs % ncores

    log.info(
        f'running in {ncores} cores with '
        f'{remaining_confs} remaining confs'
        )

    # we use a dictionary because fragments will be evaluated to exact match
    global ANGLES, BEND_ANGS, BOND_LENS, SLICEDICT_XMERS, XMERPROBS, GET_ADJ
    
    xmer_probs_tmp = prepare_xmer_probs(xmer_probs)

    # set up the information from CSSS.JSON files
    csss_dict = False
    csss_dssp_regexes = None

    all_valid_ss_codes = ''.join(dssp_ss_keys.valid)

    # There are four possibilities of sampling:
    # 1) Sampling loops and/or helix and/or strands, where the found fragments
    #    are all of the same secondary structure
    # 2) sample "any". Disregards any secondary structure annotated
    # 3) custom sample given by the user
    # 4) advanced sampling
    #
    # The following if/else block creates the needed variables according to each
    # scenario.

    if dany:
        # will sample the database disregarding the SS annotation
        dssp_regexes = [all_valid_ss_codes]

    elif custom_sampling:
        csss_dict, csss_dssp_regexes = parse_CSSS(custom_sampling)

        # If the user wants to sample "any" for some residues
        # users can have "X" in the CSSS.JSON but that will be converted
        # internally below
        if "X" in csss_dssp_regexes:
            csss_dssp_regexes.remove("X")
            csss_dssp_regexes.add(all_valid_ss_codes)
            for _k, _v in csss_dict.items():
                # X means any SS.
                if "X" in _v:
                    _v[all_valid_ss_codes] = _v.pop("X")

        dssp_regexes = list(csss_dssp_regexes)

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

    db = read_dictionary_from_disk(database)

    if bgeo_strategy == bgeo_exact_name:
        try:
            _, ANGLES, BEND_ANGS, BOND_LENS, secondary, primary = aligndb(db, True)  # noqa: E501
        except KeyError:
            log.info(S('!!!!!!!!!!!!!!!'))
            log.info(S(
                'DATABASE ERROR: '
                'the `database` requested is invalid. Please give the database '
                'generated with `bgeodb`. See the usage documentation for '
                'details while using `--bgeo-strategy exact`.'
                ))
            return

    else:
        _, ANGLES, secondary, primary = aligndb(db)

    del db

    if residue_tolerance is not None:
        _restol = str(residue_tolerance)[1:-1]
        log.info(S(f"Building with residue tolerances: {_restol}"))
    
    # create different random seeds for the different cores
    # seeds created to the cores based on main seed are predictable
    for i in range(ncores + bool(remaining_confs)):
        RANDOMSEEDS.put(random_seed + i)
        
    # creates a queue of numbers that will serve all subprocesses.
    # Used to name the output files, conformer_1, conformer_2, ...
    for i in range(1, nconfs + 1):
        CONF_NUMBER.put(i)
    
    # get sidechain dedicated parameters
    sidechain_parameters = \
        get_sidechain_packing_parameters(kwargs, sidechain_method)
    
    if long:
        csss_multi_dict = []  # list of csss_dict for each frag
        if custom_sampling:
            counter = 1
            for seq in LONG_FRAGMENTS:
                temp_csss = {}
                for res, _ in enumerate(seq):
                    temp_csss[str(res + 1)] = csss_dict[str(counter)]
                    counter += 1
                counter -= 2
                csss_multi_dict.append(temp_csss)
                assert len(temp_csss) == len(seq)
        else:
            for _ in LONG_FRAGMENTS:
                csss_multi_dict.append(False)
                    
        for idx, seq in enumerate(LONG_FRAGMENTS):
            log.info(S(f"Preparing database for sequence: {seq}"))
            SLICEDICT_XMERS = prepare_slice_dict(
                primary,
                seq,
                csss=bool(csss_multi_dict[idx]),
                dssp_regexes=dssp_regexes,
                secondary=secondary,
                mers_size=xmer_probs_tmp.sizes,
                res_tolerance=residue_tolerance,
                ncores=ncores,
                )
            remove_empty_keys(SLICEDICT_XMERS)
            # updates user defined fragment sizes and probabilities to the
            # ones actually observed
            _ = compress_xmer_to_key(xmer_probs_tmp, sorted(SLICEDICT_XMERS.keys()))  # noqa: E501
            XMERPROBS = _.probs
            
            GET_ADJ = get_adjacent_angles(
                sorted(SLICEDICT_XMERS.keys()),
                XMERPROBS,
                seq,
                ANGLES,
                bgeo_strategy,
                SLICEDICT_XMERS,
                csss=csss_multi_dict[idx],
                residue_tolerance=residue_tolerance,
                )
        
            log.info(S("done"))
            
            populate_globals(
                input_seq=seq,
                bgeo_strategy=bgeo_strategy,
                bgeo_path=bgeo_path,
                forcefield=forcefields[forcefield],
                **kwargs)
            
            ENERGYLOGSAVER.start(output_folder.joinpath(energy_log))
            
            # first run, need to generate chains first
            if seq == LONG_FRAGMENTS[0]:
                # prepars execution function
                consume = partial(
                    _build_conformers,
                    input_seq=seq,  # string
                    output_folder=output_folder,
                    nconfs=conformers_per_core,  # int
                    sidechain_parameters=sidechain_parameters,
                    sidechain_method=sidechain_method,  # goes back to kwards
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
                start = time()
            else:
                consume = partial(
                    _build_conformers,
                    input_seq=seq,  # string
                    output_folder=output_folder,
                    long=long,
                    nconfs=conformers_per_core,  # int
                    sidechain_parameters=sidechain_parameters,
                    sidechain_method=sidechain_method,  # goes back to kwards
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
                
            if remaining_confs:
                execute(conformers_per_core * ncores, nconfs=remaining_confs)
            
            # reinitialize queues so reiteration doesn't crash
            for i in range(ncores + bool(remaining_confs)):
                RANDOMSEEDS.put(random_seed + i)
            for i in range(1, nconfs + 1):
                CONF_NUMBER.put(i)
    else:
        # these are the slices with which to sample the ANGLES array
        SLICEDICT_XMERS = prepare_slice_dict(
            primary,
            input_seq,
            csss=bool(csss_dict),
            dssp_regexes=dssp_regexes,
            secondary=secondary,
            mers_size=xmer_probs_tmp.sizes,
            res_tolerance=residue_tolerance,
            ncores=ncores,
            )

        remove_empty_keys(SLICEDICT_XMERS)
        # updates user defined fragment sizes and probabilities
        # to the ones actually observed
        _ = compress_xmer_to_key(xmer_probs_tmp, sorted(SLICEDICT_XMERS.keys()))  # noqa: E501
        XMERPROBS = _.probs

        GET_ADJ = get_adjacent_angles(
            sorted(SLICEDICT_XMERS.keys()),
            XMERPROBS,
            input_seq,
            ANGLES,
            bgeo_strategy,
            SLICEDICT_XMERS,
            csss_dict,
            residue_tolerance=residue_tolerance,
            )

        populate_globals(
            input_seq=input_seq,
            bgeo_strategy=bgeo_strategy,
            bgeo_path=bgeo_path,
            forcefield=forcefields[forcefield],
            **kwargs)
        
        ENERGYLOGSAVER.start(output_folder.joinpath(energy_log))
        
        # prepars execution function
        consume = partial(
            _build_conformers,
            input_seq=input_seq,  # string
            output_folder=output_folder,
            long=long,
            nconfs=conformers_per_core,  # int
            sidechain_parameters=sidechain_parameters,
            sidechain_method=sidechain_method,  # goes back to kwards
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

        start = time()
        with Pool(ncores) as pool:
            imap = pool.imap(execute, range(ncores))
            for _ in imap:
                pass
            
        if remaining_confs:
            execute(conformers_per_core * ncores, nconfs=remaining_confs)

    log.info(f'{nconfs} conformers built in {time() - start:.3f} seconds')
    ENERGYLOGSAVER.close()


def populate_globals(
        *,
        input_seq=None,
        bgeo_strategy=bgeo_strategies_default,
        bgeo_path=None,
        forcefield=None,
        **efunc_kwargs):
    """
    Populate global variables needed for building.

    Currently, global variables include:

    BGEO_full
    BGEO_trimer
    BGEO_res
    ALL_ATOM_LABELS, ALL_ATOM_MASKS, ALL_ATOM_EFUNC
    TEMPLATE_LABELS, TEMPLATE_MASKS, TEMPLATE_EFUNC
    INT2CART

    Parameters
    ----------
    bgeo_strategy : str
        A key from the
        :py:data:`idpconfgen.components.bgeo_strategies.bgeo_strategies`.

    forcefield : str
        A key in the `core.build_definitions.forcefields` dictionary.
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
        input_seq=None,
        conformer_name='conformer',
        output_folder=None,
        long=False,
        nconfs=1,
        sidechain_parameters=None,
        bgeo_strategy=bgeo_strategies_default,
        **kwargs,
        ):
    """Arrange building of conformers and saves them to PDB files."""
    ROUND = np.round
    
    # TODO: this has to be parametrized for the different HIS types
    input_seq_3_letters = translate_seq_to_3l(input_seq)
    
    builder = conformer_generator(
        input_seq=input_seq,
        random_seed=RANDOMSEEDS.get(),
        sidechain_parameters=sidechain_parameters,
        bgeo_strategy=bgeo_strategy,
        **kwargs)
    
    atom_labels, residue_numbers, residue_labels = next(builder)
    
    if long:
        for _ in range(nconfs):
            conf_number = CONF_NUMBER.get()
            prev_struc_name = f'{conformer_name}_{conf_number}.pdb'
            prev_struc = Structure(Path(output_folder, prev_struc_name))
            prev_struc.build()
            atom_names = prev_struc.data_array[:, col_name]
            prev_seq = prev_struc.data_array[:, col_resSeq].astype(int)
            last_seq = prev_seq[-1]
            
            terminal_idx = {}
            for j, _atom in enumerate(atom_names):
                k = len(atom_names) - 1 - j
                curr_seq = prev_seq[k]
                
                if curr_seq == last_seq and atom_names[k] == "N":
                    terminal_idx["N"] = k
                elif curr_seq == last_seq and atom_names[k] == "CA":
                    terminal_idx["CA"] = k
                elif curr_seq == last_seq - 1 and atom_names[k] == "C":
                    terminal_idx["C"] = k
                elif curr_seq == last_seq - 2:
                    break
                
            stitch_Cxyz = prev_struc.data_array[terminal_idx["C"]][cols_coords].astype(float).tolist()  # noqa: E501
            stitch_Nxyz = prev_struc.data_array[terminal_idx["N"]][cols_coords].astype(float).tolist()  # noqa: E501
            stitch_CAxyz = prev_struc.data_array[terminal_idx["CA"]][cols_coords].astype(float).tolist()  # noqa: E501
            # Coordinates of boundary to stitch to later on
            stitch_coords = np.array([stitch_Cxyz, stitch_Nxyz, stitch_CAxyz])

            while 1:
                energy, coords = next(builder)

                pdb_string = gen_PDB_from_conformer(
                    input_seq_3_letters,
                    atom_labels,
                    residue_numbers,
                    ROUND(coords, decimals=3),
                    )
                pdb_arr = parse_pdb_to_array(pdb_string)
                rotated = align_coords(pdb_arr, stitch_coords, disorder_cases[2])  # noqa: E501
                clashes, fragment = count_clashes(
                    rotated,
                    prev_struc,
                    case=disorder_cases[2],
                    max_clash=40,
                    tolerance=0.4,
                    )

                if type(clashes) is int:
                    success_frag = structure_to_pdb(fragment)
                    fname_temp = f'{conformer_name}_{conf_number}_frag.pdb'
                    with open(Path(output_folder, fname_temp), 'w') as fout:
                        for line in success_frag:
                            fout.write(line + "\n")
                    final = psurgeon(
                        {"A": [Path(output_folder, fname_temp)]},
                        Path(output_folder, prev_struc_name),
                        {"A": [disorder_cases[2]]},
                        {"A": [(0, 0)]},  # placeholder
                        )
                    final_struc = structure_to_pdb(final)
                    os.remove(Path(output_folder, fname_temp))
                    os.remove(Path(output_folder, prev_struc_name))
                    with open(Path(output_folder, prev_struc_name), 'w') as fout:  # noqa: E501
                        for line in final_struc:
                            fout.write(line + "\n")
                    ENERGYLOGSAVER.save(prev_struc_name, energy)
                    break
    else:
        for _ in range(nconfs):
            energy, coords = next(builder)

            pdb_string = gen_PDB_from_conformer(
                input_seq_3_letters,
                atom_labels,
                residue_numbers,
                ROUND(coords, decimals=3),
                )

            fname = f'{conformer_name}_{CONF_NUMBER.get()}.pdb'

            with open(Path(output_folder, fname), 'w') as fout:
                fout.write(pdb_string)

            ENERGYLOGSAVER.save(fname, energy)

    del builder
    return


# the name of this function is likely to change in the future
def conformer_generator(
        *,
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

    `conformer_generator` is actually a Python generator. Examples on
    how it works:

    Note that all arguments are **named** arguments.

    >>> builder = conformer_generator(
    >>>    input_seq='MGAETTWSCAAA'  # the primary sequence of the protein
    >>>    )

    `conformer_generator` is a generator, you can instantiate it simply
    providing the residue sequence of your protein of interest.

    The **very first** iteration will return the labels of the protein
    being built. Labels are sorted by all atom models. Likewise,
    `residue_number` and `residue_labels` sample **all atoms**. These
    three are numpy arrays and can be used to index the actual coordinates.

    >>> atom_labels, residue_numbers, residue_labels = next(builder)

    After this point, each iteraction `next(builder)` yields the coordinates
    for a new conformer. There is no limit in the generator.

    >>> new_coords = next(builder)

    `new_coords` is a (N, 3) np.float64 array where N is the number of
    atoms. As expected, atom coordinates are aligned with the labels
    previously generated.

    When no longer needed,

    >>> del builder

    Should delete the builder generator.

    You can gather the coordinates of several conformers in a single
    multi dimensional array with the following:

    >>> builder = conformer_generator(
    >>>     input_seq='MGGGGG...',
    >>>     generative_function=your_function)
    >>>
    >>> atoms, res3letter, resnums = next(builder)
    >>>
    >>> num_of_conformers = 10_000
    >>> shape = (num_of_conformers, len(atoms), 3)
    >>> all_coords = np.empty(shape, dtype=float64)
    >>>
    >>> for i in range(num_of_conformers):
    >>>     all_coords[i, :, :] = next(builder)
    >>>

    Parameters
    ----------
    input_seq : str, mandatory
        The primary sequence of the protein being built in FASTA format.
        `input_seq` will be used to generate the whole conformers' and
        labels arrangement.
        Example: "MAGERDDAPL".

    generative_function : callable, optional
        The generative function used by the builder to retrieve torsion
        angles during the building process.

        The builder expects this function to receive two parameters:
            - `nres`, the residue fragment size to get angles from
            - `cres`, the next residue being built. For example,
                with cres=10, the builder will expect a minimum of three
                torsion angles (phi, psi, omega) for residue 10.

        Depending on the nature of the `generative function` the two
        pameters may be ignored by the function itself (use **kwargs
        for that purpose).

        If `None` provided, the builder will use the internal `SLIDES`
        and `ANGLES` variables and will assume the `cli_build.main` was
        executed priorly, or that ANGLES and SLICES were populated
        properly.

    disable_sidechains : bool
        Disables sidechain creation. Defaults to `False`, computes
        sidechains.

    nconfs : int
        The number of conformers to build.

    sidechain_method : str
        The method used to build/pack sidechains over the backbone
        structure. Defaults to `faspr`.
        Expects a key in
        `components.sidechain_packing.sidechain_packing_methods`.

    bgeo_strategy : str
        The strategy used to generate the bond geometries. Available options
        are: :py:data:`idpconfgen.components.bgeo_strategies.bgeo_strategies`.

    bgeo_path : str of Path
        Path to a bond geometry library as created by `bgeo` CLI.

    Yields
    ------
    First yield: tuple (np.ndarray, np.ndarray, np.ndarray)
        The conformer label arrays.

    Other yields: tuple (float, np.ndarray)
        Energy of the conformer, conformer coordinates.
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
    bb_NH_nums = np.arange(3, (len(template_input_seq) - 1) * 3 + 1, 3)[non_pro]
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
                if bgeo_strategy == bgeo_exact_name:
                    primer_template, agls, bangs, blens = GET_ADJ(bbi - 1)
                else:
                    primer_template, agls = GET_ADJ(bbi - 1)

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


def gen_PDB_from_conformer(
        input_seq_3_letters,
        atom_labels,
        residues,
        coords,
        ALF=atom_line_formatter,
        ):
    """."""
    lines = []
    LINES_APPEND = lines.append
    ALF_FORMAT = ALF.format
    resi = -1

    # this is possible ONLY because there are no DOUBLE CHARS atoms
    # in the atoms that constitute a protein chain
    ATOM_LABEL_FMT = ' {: <3}'.format

    assert len(atom_labels) == coords.shape[0]

    atom_i = 1
    for i in range(len(atom_labels)):

        if np.isnan(coords[i, 0]):
            continue

        if atom_labels[i] == 'N':
            resi += 1
            current_residue = input_seq_3_letters[resi]
            current_resnum = residues[i]

        atm = atom_labels[i].strip()
        ele = atm.lstrip('123')[0]

        if len(atm) < 4:
            atm = ATOM_LABEL_FMT(atm)

        LINES_APPEND(ALF_FORMAT(
            'ATOM',
            atom_i,
            atm,
            '',
            current_residue,
            'A',
            current_resnum,
            '',
            coords[i, 0],
            coords[i, 1],
            coords[i, 2],
            0.0,
            0.0,
            '',
            ele,
            '',
            ))

        atom_i += 1

    return os.linesep.join(lines)


def get_adjacent_angles(
        options,
        probs,
        seq,
        dihedrals_db,
        bgeo_strategy,
        slice_dict,
        csss,
        residue_tolerance=None,
        ):
    """
    Get angles to build the next adjacent protein fragment.

    Parameters
    ----------
    options : list
        The length of the possible fragment sizes.

    probs : list
        A list with the relative probabilites to select from `options`.

    seq : str
        The conformer sequence.

    dihedrals_db : dict-like
        The angle omega/phi/psi database.

    bgeo_strategy : string
        Bond geometry strategy to use.

    slice_dict : dict-like
        A dictionary containing the fragments strings as keys and as values
        lists with slice objects.

    csss : dict-like
        A dictionary containing probabilities of secondary structures per
        amino acid residue position.
    
    residue_tolerance : dict-like
        A dictionary for possible residue tolerances to look for while sampling
        torsion angles of amino acids.
    """
    residue_tolerance = residue_tolerance or {}
    probs = fill_list(probs, 0, len(options))

    # prepares helper lists
    lss = []  # list of possible secondary structures in case `csss` is given
    lssprobs = []  # list of possible ss probabilities in case `csss` is given
    lssE, lssprobsE = lss.extend, lssprobs.extend
    lssC, lssprobsC = lss.clear, lssprobs.clear

    def func(
            aidx,
            CRNFI=calc_residue_num_from_index,
            RC=np.random.choice,
            GSCNJIT=get_seq_chunk_njit,
            BRS=build_regex_substitutions,
            ):

        # calculates the current residue number from the atom index
        cr = CRNFI(aidx)
        
        # chooses the size of the fragment from
        # pre-configured range of sizes
        plen = RC(options, p=probs)
            
        # defines the fragment identity accordingly
        primer_template = GSCNJIT(seq, cr, plen)
        _ori_template = primer_template
        next_residue = GSCNJIT(seq, cr + plen, 1)
        # recalculates the plen to avoid plen/template inconsistencies that
        # occur if the plen is higher then the number of
        # residues until the end of the protein.
        plen = len(primer_template)
        pt_sub = BRS(primer_template, residue_tolerance)

        while plen > 0:
            if next_residue == 'P':
                pt_sub = f'{pt_sub}_P'

            try:
                if csss:
                    # matches current residue
                    # to build with residue number in CSSS
                    cr_plus_1 = str(cr + 1)
                    # clear lists
                    lssC()
                    lssprobsC()
                    # adds possible secondary structure for the residue
                    # the first residue of the fragment
                    lssE(csss[cr_plus_1].keys())
                    # adds SS probabilities for the same residue
                    lssprobsE(csss[cr_plus_1].values())
                    # based on the probabilities,
                    # select a SS for residue in question
                    pcsss = RC(lss, p=lssprobs)
                    _slice = RC(slice_dict[plen][pt_sub][pcsss])
                else:
                    _slice = RC(slice_dict[plen][pt_sub])
                
                dihedrals = dihedrals_db[_slice, :].ravel()
                
                if bgeo_strategy == bgeo_exact_name:
                    bend_angs = BEND_ANGS[_slice, :].ravel()
                    bond_lens = BOND_LENS[_slice, :].ravel()

            except (KeyError, ValueError):
                # walks back one residue
                plen -= 1
                next_residue = primer_template[-1]
                primer_template = primer_template[:-1]
                pt_sub = BRS(primer_template, residue_tolerance)
            else:
                break
        else:
            # raise AssertionError to avoid `python -o` silencing
            _emsg = (
                "The code should not arrive here. "
                "If it does, it may mean no matches were found for fragment "
                f"{_ori_template!r} down to the single residue."
                )
            raise AssertionError(_emsg)

        if next_residue == 'P':
            # because angles have the proline information

            if bgeo_strategy == bgeo_exact_name:
                return primer_template + 'P', dihedrals, bend_angs, bond_lens

            return primer_template + 'P', dihedrals

        else:
            if bgeo_strategy == bgeo_exact_name:
                return primer_template, dihedrals, bend_angs, bond_lens

            return primer_template, dihedrals

    return func


if __name__ == "__main__":
    libcli.maincli(ap, main)
