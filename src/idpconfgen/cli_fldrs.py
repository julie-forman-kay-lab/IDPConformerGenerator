"""
Client for filling in the tails of folded proteins with IDRs.

Build from a database of torsion angles and secondary structure
information. Database is as created by `idpconfgen torsions` CLI.

Future Ideas
------------
- Populate internal disordered regions (code has been scaffolded already)
- Implement CSSS on IDR tails/fragments

USAGE:
    $ idpconfgen fldrs -db torsions.json -seq sequence.fasta -fld folded.pdb
"""
import argparse
import shutil
from functools import partial
from itertools import cycle
from multiprocessing import Pool, Queue
from random import randint
from time import time

import numpy as np

from idpconfgen import Path, log
from idpconfgen.cli_build import gen_PDB_from_conformer, get_adjacent_angles
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
from idpconfgen.fldrs_helper import (
    align_coords,
    break_check,
    consecutive_grouper,
    count_clashes,
    create_combinations,
    disorder_cases,
    psurgeon,
    store_idp_paths,
    tolerance_calculator,
    )
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
from idpconfgen.libs.libfilter import aligndb
from idpconfgen.libs.libhigherlevel import bgeo_reduce
from idpconfgen.libs.libio import make_folder_or_cwd, read_dictionary_from_disk
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.libs.libparse import (
    get_trimer_seq_njit,
    remap_sequence,
    remove_empty_keys,
    translate_seq_to_3l,
    )
from idpconfgen.libs.libpdb import get_fasta_from_PDB
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    cols_coords,
    parse_pdb_to_array,
    structure_to_pdb,
    )
from idpconfgen.logger import S, T, init_files, pre_msg, report_on_crash


_file = Path(__file__).myparents()
LOGFILESNAME = '.idpconfgen_fldrs'

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

# Case of disorder region (1, 2, or 3) determines strategy
# also set defaults to be none assuming full IDP
DISORDER_CASE = None
DISORDER_BOUNDS = None
DISORDER_INDEX = None
DISORDER_SEQS = None

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
_name = 'fldrs'
_help = 'Building IDRs within a given folded protein.'

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

#########################################
libcli.add_argument_dloopoff(ap)
libcli.add_argument_dhelix(ap)
libcli.add_argument_dstrand(ap)
libcli.add_argument_dany(ap)
libcli.add_argument_duser(ap)
#########################################

ap.add_argument(
    '-fld',
    '--folded-structure',
    help="Input .PDB file for folded structure of interest.",
    required=True,
    )

ap.add_argument(
    '-kt',
    '--keep-temporary',
    help=(
        "Switch to keep temporary disordered fragments while building."
        ),
    action='store_true',
    )

ap.add_argument(
    '-tol',
    '--clash-tolerance',
    help=(
        "Float value clash tolerance between 0.0-1.0 where 0.5 is the default "
        "value denoting minimum of 50 spherical clashes, 0.5 Angstroms "
        "of tolerance with a given vdW radii. Where 1.0 allows for 100 "
        "clashes and 1.0 Angstroms of tolerance for a given vdW radii."
        ),
    default=0.5,
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
        'when building the BACKBONE atoms. Defaults to 10.'
        ),
    default=10.0,
    type=float,
    )

ap.add_argument(
    '-etss',
    '--energy-threshold-sidechains',
    help=(
        'The energy threshold above which conformers will be rejected '
        'after packing the sidechains (ignored if `-dsd`). '
        'Defaults to 1000.'
        ),
    default=1000.0,
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


def main(
        input_seq,
        database,
        clash_tolerance=0.5,
        keep_temporary=False,
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
    try:
        clash_tolerance = float(clash_tolerance)
    except ValueError:
        log.info(
            "Please enter a floating point value for clash-tolerances "
            "between 0.0-1.0."
            )
        return
    
    max_clash, dist_tolerance = tolerance_calculator(clash_tolerance)
    
    if keep_temporary:
        TEMP_DIRNAME = "fldrs_temp/"
    else:
        TEMP_DIRNAME = '.fldrs_temp/'
    
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

    output_folder = make_folder_or_cwd(output_folder)
    init_files(log, Path(output_folder, LOGFILESNAME))
    log.info(T('starting the building process'))
    log.info(S(f'input sequence: {input_seq}'))
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
    
    # set the boundaries of disordered regions for folded proteins
    global DISORDER_CASE, DISORDER_BOUNDS, DISORDER_SEQS

    xmer_probs_tmp = prepare_xmer_probs(xmer_probs)

    all_valid_ss_codes = ''.join(dssp_ss_keys.valid)

    # There are four possibilities of sampling:
    # 1) Sampling loops and/or helix and/or strands, where the found fragments
    #    are all of the same secondary structure
    # 2) sample "any". Disregards any secondary structure annotated
    # 3) custom sample given by the user (CURRENTLY NOT IN FLDRS)
    # 4) advanced sampling
    #
    # The following if/else block creates the needed variables according to each
    # scenario.

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
        
    # TODO: in the future, can give a tarball or folder of .PDB
    # randomly select which one to attach disordered regions onto
    assert folded_structure.endswith('.pdb')
    
    log.info(T('Initializing folded domain information'))
    
    with open(folded_structure) as f:
        pdb_raw = f.read()
    _pdbid, fld_fasta = get_fasta_from_PDB([folded_structure, pdb_raw])
    fld_fasta = fld_fasta.replace('X', '')
    
    # Find out what our disordered sequences are
    # Can be C-term, N-term, in the middle, or all the above
    breaks = break_check(pdb_raw)
    mod_input_seq = input_seq
    if breaks:
        for seq in breaks:
            mod_input_seq = mod_input_seq.replace(seq, len(seq) * '*')
            
    mod_input_seq = mod_input_seq.replace(fld_fasta, len(fld_fasta) * '*')
    assert len(mod_input_seq) == len(input_seq)
    # Treats folded residues as "*" to ignore in IDP building process
    disordered_res = []
    for i, res in enumerate(mod_input_seq):
        if res != "*":
            disordered_res.append(i)
    
    DISORDER_BOUNDS = consecutive_grouper(disordered_res)
    DISORDER_SEQS = []
    for i, bounds in enumerate(DISORDER_BOUNDS):
        lower = bounds[0]
        upper = bounds[1]
        dis_seq = input_seq[lower:upper]
        log.info(S(
            f'Disordered region #{i + 1} = residue {lower + 1} to {upper} '
            f'with the sequence {dis_seq}.'
            ))
        
        # Different than what we tell user due to internal processing for
        # disordered bits (need 2 extra residues at the beginning or end)
        if lower == 0:
            dis_seq = input_seq[lower: upper + 2]
        elif upper == len(input_seq):
            dis_seq = input_seq[lower - 2: upper]
        else:
            dis_seq = input_seq[lower: upper]
        
        DISORDER_SEQS.append(dis_seq)
    log.info(S('done'))

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

    # these are the slices with which to sample the ANGLES array
    SLICEDICT_XMERS = []
    XMERPROBS = []
    GET_ADJ = []
    for i, seq in enumerate(DISORDER_SEQS):
        log.info(S(f"Preparing database for sequence: {seq}"))
        
        SLICEDICT_XMERS.append(prepare_slice_dict(
            primary,
            seq,
            dssp_regexes=dssp_regexes,
            secondary=secondary,
            mers_size=xmer_probs_tmp.sizes,
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

    # create different random seeds for the different cores
    # seeds created to the cores based on main seed are predictable
    for i in range(ncores + bool(remaining_confs)):
        RANDOMSEEDS.put(random_seed + i)

    # creates a queue of numbers that will serve all subprocesses.
    # Used to name the output files, conformer_1, conformer_2, ...
    for i in range(1, nconfs + 1):
        CONF_NUMBER.put(i)

    ENERGYLOGSAVER.start(output_folder.joinpath(energy_log))

    # get sidechain dedicated parameters
    sidechain_parameters = \
        get_sidechain_packing_parameters(kwargs, sidechain_method)
    
    fld_struc = Structure(Path(folded_structure))
    fld_struc.build()
    atom_names = fld_struc.data_array[:, col_name]
    # Generate library of conformers for each case
    for i, seq in enumerate(DISORDER_SEQS):
        fld_term_idx = {}
        lower = DISORDER_BOUNDS[i][0]
        upper = DISORDER_BOUNDS[i][1]
        counter = 0
        # Use the second folded residue's backbone C-N-CA
        # for rotation and alignment where the `C` atom
        # is from the previous resiude
        if lower == 0:  # N-IDR case
            DISORDER_CASE = disorder_cases[0]
            for j, atom in enumerate(atom_names):
                if counter == 5:
                    break
                if atom == "N":
                    fld_term_idx["N"] = j
                    counter += 1
                elif atom == "CA":
                    fld_term_idx["CA"] = j
                    counter += 1
                elif atom == "C":
                    fld_term_idx["C"] = j
                    counter += 1
        elif upper == len(input_seq):  # C-IDR case
            DISORDER_CASE = disorder_cases[2]
            for j, _atom in enumerate(atom_names):
                k = len(atom_names) - 1 - j
                if counter == 1 and atom_names[k] == "N":
                    fld_term_idx["N"] = k
                    counter += 1
                elif counter == 0 and atom_names[k] == "CA":
                    fld_term_idx["CA"] = k
                    counter += 1
                elif counter == 2 and atom_names[k] == "C":
                    fld_term_idx["C"] = k
                    break
        else:
            DISORDER_CASE = disorder_cases[1]  # break

        fld_Cxyz = fld_struc.data_array[fld_term_idx["C"]][cols_coords].astype(float).tolist()  # noqa: E501
        fld_Nxyz = fld_struc.data_array[fld_term_idx["N"]][cols_coords].astype(float).tolist()  # noqa: E501
        fld_CAxyz = fld_struc.data_array[fld_term_idx["CA"]][cols_coords].astype(float).tolist()  # noqa: E501
        # Coordinates of boundary to stitch to later on
        fld_coords = np.array([fld_Cxyz, fld_Nxyz, fld_CAxyz])

        populate_globals(
            input_seq=seq,
            bgeo_strategy=bgeo_strategy,
            bgeo_path=bgeo_path,
            forcefield=forcefields[forcefield],
            **kwargs)
        
        # create the temporary folders housing the disordered PDBs
        # this folder will be deleted when everything is grafted together
        temp_of = make_folder_or_cwd(
            output_folder.joinpath(TEMP_DIRNAME + DISORDER_CASE)
            )
        
        log.info(S(f"Generating temporary disordered conformers for: {seq}"))
        log.info(S(
            "Please note that sequence may contain 2 extra residues "
            "to facilitate reconstruction later on."
            ))
        
        fStruct = Structure(Path(folded_structure))
        fStruct.build()
        
        # prepares execution function
        consume = partial(
            _build_conformers,
            fld_xyz=fld_coords,
            fld_struc=fStruct,
            disorder_case=DISORDER_CASE,
            max_clash=max_clash,
            tolerance=dist_tolerance,
            index=i,
            conformer_name="conformer_" + DISORDER_CASE,
            input_seq=seq,  # string
            output_folder=temp_of,
            nconfs=conformers_per_core,  # int
            sidechain_parameters=sidechain_parameters,
            sidechain_method=sidechain_method,  # goes back to kwargs
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
        
        log.info(f'{nconfs} {DISORDER_CASE} conformers built in {time() - start:.3f} seconds')  # noqa: E501
        
        # reinitialize queues so reiteration doesn't crash
        for i in range(ncores + bool(remaining_confs)):
            RANDOMSEEDS.put(random_seed + i)
        for i in range(1, nconfs + 1):
            CONF_NUMBER.put(i)

    DISORDER_CASE = store_idp_paths(output_folder, TEMP_DIRNAME)
    
    # After conformer construction and moving, time to graft proteins together
    # - Keep an eye out for `col_chainID` and `col_segid` -> make all chain A
    # - Note that residue number `col_resSeq` needs to be continuous
    # - For good form, make sure `col_serial` is consistent as well
    # - When grafting remove the tether residue on donor chain
    # - Generate a tuple database of which pairs have already been generated
    if len(DISORDER_CASE) == 1:
        case = next(iter(DISORDER_CASE))
        files = DISORDER_CASE[case]
    elif disorder_cases[0] in DISORDER_CASE and disorder_cases[2] in DISORDER_CASE:  # noqa: E501
        case = disorder_cases[0] + disorder_cases[2]
        files = create_combinations(
            DISORDER_CASE[disorder_cases[0]],
            DISORDER_CASE[disorder_cases[2]],
            nconfs,
            )
        
    log.info("Stitching conformers onto the folded domain...")
    consume = partial(
        psurgeon,
        fld_struc=Path(folded_structure),
        case=case,
        )
    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_folder=output_folder,
        ROC_prefix=_name,
        )
    execute_pool = pool_function(execute, files, ncores=ncores)
    
    for i, conf in enumerate(execute_pool):
        struc = structure_to_pdb(conf)
        output = output_folder.joinpath(f"conformer_{i + 1}.pdb")
        with open(output, 'w') as f:
            for line in struc:
                f.write(line + "\n")
    
    if not keep_temporary:
        shutil.rmtree(output_folder.joinpath(TEMP_DIRNAME))

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
        fld_xyz=None,
        fld_struc=None,
        disorder_case=None,
        max_clash=50,
        tolerance=0.5,
        index=None,
        input_seq=None,
        conformer_name='conformer',
        output_folder=None,
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
        index=index,
        input_seq=input_seq,
        random_seed=RANDOMSEEDS.get(),
        sidechain_parameters=sidechain_parameters,
        bgeo_strategy=bgeo_strategy,
        **kwargs)

    atom_labels, residue_numbers, _residue_labels = next(builder)

    for _ in range(nconfs):
        while 1:
            energy, coords = next(builder)

            pdb_string = gen_PDB_from_conformer(
                input_seq_3_letters,
                atom_labels,
                residue_numbers,
                ROUND(coords, decimals=3),
                )

            pdb_arr = parse_pdb_to_array(pdb_string)
            rotated = align_coords(pdb_arr, fld_xyz, disorder_case)
            clashes, fragment = count_clashes(
                rotated,
                fld_struc,
                disorder_case,
                max_clash,
                tolerance,
                )
            
            if type(clashes) is int:
                final = structure_to_pdb(fragment)
                log.info(f"Succeeded {disorder_case} to folded region clash check!")  # noqa: E501
                break
            else:
                log.info(f"Failed {disorder_case} to folded region clash check... regenerating")  # noqa: E501
        
        fname = f'{conformer_name}_{CONF_NUMBER.get()}.pdb'

        with open(Path(output_folder, fname), 'w') as fout:
            for line in final:
                fout.write(line + "\n")

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
                    primer_template, agls, bangs, blens = GET_ADJ[index](bbi - 1)  # noqa: E501
                else:
                    primer_template, agls = GET_ADJ[index](bbi - 1)

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
