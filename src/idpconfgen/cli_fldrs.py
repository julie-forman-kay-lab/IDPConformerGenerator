"""
Client for filling in the blanks of folded proteins with IDRs.

Build from a database of torsion angles and secondary structure
information. Database is as created by `idpconfgen torsions` CLI.

USAGE:
    $ idpconfgen fldrs -db torsions.json -seq MMMMMMM...

"""
import argparse
import os
from functools import partial
from itertools import cycle
from multiprocessing import Pool, Queue
from random import randint
from time import time

import numpy as np

from idpconfgen import Path, log
from idpconfgen.cli_build import (
    conformer_generator,
    get_adjacent_angles,
    gen_PDB_from_conformer,
    parse_CSSS,
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
    translate_seq_to_3l,
    )
from idpconfgen.libs.libpdb import atom_line_formatter, get_fasta_from_PDB
from idpconfgen.logger import S, T, init_files, pre_msg, report_on_crash


from idpconfgen.fldrs_helper import (
    disorder_cases,
    break_check,
    consecutive_grouper,
    )

_file = Path(__file__).myparents()
LOGFILESNAME = '.idpconfgen_fldrs'
TEMP_DIRNAME = '.fldrs_temp'

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
    '-flds',
    '--folded-structure',
    help="Input .PDB file for folded structure of interest.",
    required=True,
    )

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
        custom_sampling,
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
        first = bounds[0]
        last = bounds[1]
        dis_seq = input_seq[first:last]
        log.info(S(
            f'Disordered region #{i + 1} = residue {first + 1} to {last} '
            f'with the sequence {dis_seq}.'
        ))
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
        SLICEDICT_XMERS.append(prepare_slice_dict(
            primary,
            seq,
            csss=bool(csss_dict),
            dssp_regexes=dssp_regexes,
            secondary=secondary,
            mers_size=xmer_probs_tmp.sizes,
            res_tolerance=residue_tolerance,
            ncores=ncores,
            ))

        remove_empty_keys(SLICEDICT_XMERS[i])
        # updates user defined fragment sizes and probabilities to the ones actually
        # observed
        _ = compress_xmer_to_key(xmer_probs_tmp, sorted(SLICEDICT_XMERS[i].keys()))
        XMERPROBS.append(_.probs)

        GET_ADJ.append(get_adjacent_angles(
            sorted(SLICEDICT_XMERS[i].keys()),
            XMERPROBS[i],
            seq,
            ANGLES,
            bgeo_strategy,
            SLICEDICT_XMERS[i],
            csss_dict,
            fld_slice_dict=None,
            fld_dihedrals_db=None,
            residue_tolerance=residue_tolerance,
            ))

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


    # Generate library of conformers for each case
    for i, seq in enumerate(DISORDER_SEQS):
        populate_globals(
            input_seq=seq,
            bgeo_strategy=bgeo_strategy,
            bgeo_path=bgeo_path,
            forcefield=forcefields[forcefield],
            **kwargs)
        
        if DISORDER_BOUNDS[i][0] == 0:
            DISORDER_CASE = disorder_cases[0] 
        elif DISORDER_BOUNDS[i][1] == len(input_seq):
            DISORDER_CASE = disorder_cases[2]
        else:
            DISORDER_CASE = disorder_cases[1]
        
        # TODO create the temporary folders housing the PDBs like savepoints
        # this folder will be deleted when everything is grafted together
        
        # prepares execution function
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
        nconfs=1,
        sidechain_parameters=None,
        bgeo_strategy=bgeo_strategies_default,
        **kwargs,
        ):
    """Arrange building of conformers and saves them to PDB files."""
    ROUND = np.round

    # TODO: this has to be parametrized for the different HIS types
    input_seq_3_letters = translate_seq_to_3l(input_seq)

    # TODO: strategy for FLDR/S
    # - generate each disordered region independently and try to glue back on
    #   folded PDB using `psurgeon`
    # - we want to minimize modification of `conformer_generator` to retain
    #   native sampling diversity
    # - use coordinates of atoms on folded protein instead of "dummy atoms"
    
    builder = conformer_generator(
        input_seq=input_seq,
        random_seed=RANDOMSEEDS.get(),
        sidechain_parameters=sidechain_parameters,
        bgeo_strategy=bgeo_strategy,
        **kwargs)

    atom_labels, residue_numbers, residue_labels = next(builder)

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


if __name__ == "__main__":
    libcli.maincli(ap, main)
