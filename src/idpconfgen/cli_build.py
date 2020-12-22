"""
Builds IDP conformers.

Build from a database of torsion angles and secondary structure
information. Database is as created by `idpconfgen torsions` CLI.

USAGE:
    $ idpconfgen build -db torsions.json -seq MMMMMMM...

"""
import argparse
import re
import sys
from collections import Counter
from functools import partial
from itertools import cycle
from multiprocessing import Pool
from numbers import Number
from random import choice as randchoice
from random import randint
from time import time

import numpy as np
from numba import njit

from idpconfgen import log
from idpconfgen.core.build_definitions import (
    amber_pdbs,
    atom_labels_amber,
    backbone_atoms,
    bonds_equal_3_inter,
    bonds_le_2_inter,
    build_bend_H_N_C,
    build_bend_angles_CA_C_Np1,
    build_bend_angles_Cm1_N_CA,
    build_bend_angles_N_CA_C,
    distance_H_N,
    distances_CA_C,
    distances_C_Np1,
    distances_N_CA,
    expand_topology_bonds_apart,
    generate_residue_template_topology,
    inter_residue_connectivities,
    read_ff14SB_params,
    sidechain_templates,
    topology_3_bonds_apart,
    )
from idpconfgen.core.definitions import aa1to3, vdW_radii_dict
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.libs import libcli
from idpconfgen.libs.libcalc import (
    calc_all_vs_all_dists_square,
    calc_all_vs_all_dists,
    make_coord_Q,
    make_coord_Q_CO,
    make_coord_Q_COO,
    place_sidechain_template,
    )
from idpconfgen.libs.libfilter import aligndb, regex_search
from idpconfgen.libs.libio import read_dictionary_from_disk
from idpconfgen.libs.libpdb import atom_line_formatter
from idpconfgen.libs.libtimer import timeme


_name = 'build'
_help = 'Builds conformers from database.'


_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    '-db',
    '--database',
    help='The IDPConfGen database.',
    required=True,
    )

ap.add_argument(
    '-seq',
    '--input_seq',
    help='The Conformer residue sequence.',
    required=True
    )

ap.add_argument(
    '-nc',
    '--nconfs',
    help='Number of conformers to build.',
    default=1,
    type=int,
    )

ap.add_argument(
    '-dr',
    '--dssp-regexes',
    help='Regexes used to search in DSSP',
    default='(?=(L{2,6}))',
    nargs='+',
    )

ap.add_argument(
    '-dsd',
    '--disable-sidechains',
    help='Whether or not to compute sidechais. Defaults to True.',
    action='store_true',
    )

# TODO: these three parameters must be discontinued
#libcli.add_argument_vdWb(ap)
#libcli.add_argument_vdWr(ap)
#libcli.add_argument_vdWt(ap)
libcli.add_argument_ncores(ap)


SLICES = []
ANGLES = None


def main(
        input_seq,
        database,
        dssp_regexes=r'(?=(L{2,6}))',
        func=None,
        nconfs=1,
        ncores=1,
        **kwargs,
        ):
    """
    Execute main client logic.

    Distributes over processors.
    """
    global ANGLES
    # Calculates how many conformers are built per core
    if nconfs < ncores:
        ncores = 1
        core_chunks = nconfs
        remaining_chunks = 0
    else:
        core_chunks = nconfs // ncores
        # in case nconfs is not multiple of ncores, builds the remaining confs
        # at the end
        remaining_chunks = nconfs % ncores

    _slices, ANGLES = read_db_to_slices(database, dssp_regexes, ncores=ncores)
    SLICES.extend(_slices)

    # prepars execution function
    execute = partial(
        build_conformers,
        input_seq=input_seq,  # string
        nconfs=core_chunks,  # int
        **kwargs,
        )

    start = time()
    with Pool(ncores) as pool:
        imap = pool.imap(execute, range(ncores))
        for _ in imap:
            pass

    if remaining_chunks:
        execute(core_chunks * ncores, nconfs=remaining_chunks)

    log.info(f'{nconfs} conformers built in {time() - start:.3f} seconds')


def read_db_to_slices(database, dssp_regexes, ncores=1):
    """Create database base of slice and angles."""
    # reads db dictionary from disk
    db = read_dictionary_from_disk(database)
    log.info(f'Read DB with {len(db)} entries')

    # reads and prepares IDPConfGen data base
    timed = partial(timeme, aligndb)
    pdbs, angles, dssp, resseq = timed(db)

    # searchs for slices in secondary structure, according to user requests
    timed = partial(timeme, regex_search, ncores=ncores)
    dssp_regexes = \
        [dssp_regexes] if isinstance(dssp_regexes, str) else dssp_regexes
    slices = []
    for dssp_regex_string in dssp_regexes:
        slices.extend(timed(dssp, dssp_regex_string))
    log.info(f'Found {len(slices)} indexes for {dssp_regexes}')

    return slices, angles


def build_conformers(
        execution_run,
        input_seq,
        conformer_name='conformer',
        generative_function=None,
        # TODO: these parameters must be discontinued
        #vdW_bonds_apart=3,
        #vdW_tolerance=0.4,
        #vdW_radii='tsai1999',
        nconfs=1,
        disable_sidechains=True,
        ):
    """
    Build conformers.

    *Note*: `execution_run` is a positional parameter in order to maintain
    operability with the multiprocessing operations in `main`, sorry for
    that :-)

    Parameters
    ----------
    execution_run : int
        A zero or positive integer. `execution_run` will number the
        conformers accordingly to the CPU job number which generates
        them.

        `execution_run` dictate the numbering suffix of the PDBs saved
        to disk.

        If a single CPU is being used, just use 0 for convenience.

    input_seq : str
        The FASTA sequence of the protein being built. This parameter
        needs to be accurate and is extremely important to create an
        appropriate data structure able to contain all conformer's
        atoms.
        Example: "MAGERDDAPL".

    conformer_name : str, optional
        The prefix name of the PDB file saved to disk.
        Defaults to 'conformer'.

    generative_function : callable, optional
        The generative function used by the builder to retrieve torsion
        angles during the building process.

        The builder expects this function to receive two parameters:
            - `nres`, the residue chunk size to get angles from
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
    """
    BUILD_BEND_H_N_C = build_bend_H_N_C
    #CALC_DISTS = calc_all_vs_all_dists_square
    # TODO
    CALC_DISTS = calc_all_vs_all_dists
    # TODO, maybe this should be sum, if used at all: https://stackoverflow.com/questions/8364674/how-to-count-the-number-of-true-elements-in-a-numpy-bool-array
    DISTANCE_NH = distance_H_N
    LOGICAL_AND = np.logical_and
    MAKE_COORD_Q_COO_LOCAL = make_coord_Q_COO
    MAKE_COORD_Q_CO_LOCAL = make_coord_Q_CO
    MAKE_COORD_Q_LOCAL = make_coord_Q
    NAN = np.nan
    NANSUM = np.nansum
    PLACE_SIDECHAIN_TEMPLATE = place_sidechain_template
    RC = randchoice
    RINT = randint
    ROUND = np.round
    SIDECHAIN_TEMPLATES = sidechain_templates
    angles = ANGLES
    slices = SLICES
    # bellow an optimization for the builder loop

    # TODO: correct for HIS/HIE/HID/HIP
    input_seq_3_letters = [
        'HIP' if _res == 'H' else aa1to3[_res]
        for _res in input_seq
        ]

    # semantic exchange for speed al readibility
    with_sidechains = not(disable_sidechains)

    # tests generative function complies with implementation requirements
    if generative_function:
        try:
            generative_function(nres=1, cres=0)
        except Exception as err:  # this is generic Exception on purpose
            errmsg = (
                'The `generative_function` provided is not compatible with '
                'the building process. Please read `build_conformers` docstring '
                'for more details.'
                )
            raise IDPConfGenException(errmsg) from err

    # Start building process

    # prepares data based on the input sequence
    # considers sidechain all-atoms
    atom_labels = np.array(generate_atom_labels(input_seq, atom_labels_amber))
    print(atom_labels)
    num_atoms = len(atom_labels)
    print('num atoms: ', num_atoms)

    num_ij_pairs = num_atoms * (num_atoms - 1) // 2
    print('num ij pairs: ', num_ij_pairs)

    residue_numbers = np.array(generate_residue_numbers(atom_labels, start=1))
    residue_labels = np.array(generate_residue_labels(input_seq, atom_labels))

    assert len(residue_numbers) == num_atoms
    assert len(residue_labels) == num_atoms, (len(residue_labels), num_atoms)

    # /
    # Preparation for the energy function
    # TODO: parametrize this. the use should be able to chose different forcefields
    print('reading ff parameters')
    ff14SB = read_ff14SB_params()
    res_topology = generate_residue_template_topology(
        amber_pdbs,
        atom_labels_amber,
        add_OXT=True,
        add_Nterminal_H=True,
        )
    bonds_equal_3_intra = topology_3_bonds_apart(res_topology)
    bonds_le_2_intra = expand_topology_bonds_apart(res_topology, 2)

    print('Preparing ii pairs')
    sigmas_ii, epsilons_ii, charges_i = populate_ff_parameters_in_structure(
        atom_labels,
        residue_numbers,
        residue_labels,
        ff14SB,  # the forcefield
        )

    # The mask to identify ij pairs exactly 3 bonds apart is needed for the
    # special scaling factor of Coulomb and LJ equations
    # This mask will be used only aftert the calculation of the CLJ params
    print('Create 1-4 masks')
    bonds_exact_3_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        bonds_equal_3_intra,
        bonds_equal_3_inter,
        )

    print('Create disable masks')
    # this mask will disable calculations in covalently bond atoms and
    # atoms separated 2 bonds apart
    _bonds_ge_3_mask = create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        bonds_le_2_intra,
        bonds_le_2_inter,
        base_bool=True,
        )
    # on lun 21 dic 2020 17:33:31 EST, I tested for 1M sized array
    # the numeric indexing performed better than the boolean indexing
    # 25 ns versus 31 ns.
    bonds_ge_3_mask = np.where(_bonds_ge_3_mask)[0]
    del _bonds_ge_3_mask

    # /
    # Prepares Coulomb and Lennard-Jones pre computed parameters:
    # calculates ij combinations using raw njitted functions because using
    # numpy outer variantes in very large systems overloads memory and
    # reduces performance.
    #
    # sigmas
    sigmas_ij_pre = np.empty(num_ij_pairs, dtype=np.float64)
    njit_calc_sum_upper_diagonal_raw(sigmas_ii, sigmas_ij_pre)
    #
    # epsilons
    epsilons_ij_pre = np.empty(num_ij_pairs, dtype=np.float64)
    njit_calc_multiplication_upper_diagonal_raw(epsilons_ii, epsilons_ij_pre)
    #
    # charges
    charges_ij = np.empty(num_ij_pairs, dtype=np.float64)
    njit_calc_multiplication_upper_diagonal_raw(charges_i, charges_ij)

    epsilons_ij = np.sqrt(epsilons_ij_pre)  # combinatorial rule / sqrt
    del epsilons_ij_pre

    sigmas_ij = sigmas_ij_pre * 0.5  # combinatorial rule / average
    del sigmas_ij_pre

    acoeff = 4 * epsilons_ij * (sigmas_ij ** 12)
    bcoeff = 4 * epsilons_ij * (sigmas_ij ** 6)

    charges_ij *= 0.25  # dielectic constant

    # this is the Lennard-Jones special case, where scaling factors are applied
    # to atoms bonded 3 bonds apart, known as the '14' cases.
    # 0.4 was calibrated manually, until I could find a conformer
    # within 50 trials dom 20 dic 2020 13:16:50 EST
    # I believe, because we are not doing local minimization here, we
    # cannot be that strick with the 14 scaling factor, and a reduction
    # factor of 2 is not enough
    acoeff[bonds_exact_3_mask] *= float(ff14SB['lj14scale']) * 0.4
    bcoeff[bonds_exact_3_mask] *= float(ff14SB['lj14scale']) * 0.4
    charges_ij[bonds_exact_3_mask] *= float(ff14SB['coulomb14scale'])

    # TODO: needed for other energy functions
    #atoms_VDW = 0.5 * sigmas_ii * 2**(1/6)

    # creates index-translated boolean masks
    bb_mask = np.where(np.isin(atom_labels, ('N', 'CA', 'C')))[0]
    carbonyl_mask = np.where(atom_labels == 'O')[0]
    NHydrogen_mask = np.where(atom_labels == 'H')[0]
    OXT_index = np.where(atom_labels == 'OXT')[0][0]

    # the last O value is removed because it is built in the same step
    # OXT is built
    OXT1_index = carbonyl_mask[-1]
    carbonyl_mask = carbonyl_mask[:-1]

    # create coordinates and views
    coords = np.full((num_atoms, 3), NAN, dtype=np.float64)

    # +1 because of the dummy coordinate required to start building.
    # see later
    bb = np.full((bb_mask.size + 1, 3), NAN, dtype=np.float64)
    bb_real = bb[1:, :]  # backbone coordinates without the dummy
    # coordinates for the carbonyl oxigen atoms
    bb_CO = np.full((carbonyl_mask.size, 3), NAN, dtype=np.float64)

    # notice that NHydrogen_mask does not see Prolines
    bb_NH = np.full((NHydrogen_mask.size, 3), NAN, dtype=np.float64)

    # the following array serves to avoid placing HN in Proline residues
    residue_labels_bb_simulating = residue_labels[bb_mask]

    # creates seed coordinates:
    # 1st) a dummy atom at the y-axis to build the first atom
    # 2nd) N-terminal N-atom is at 0, 0, 0
    # 3rd) CA atom of the firs residue is at the x-axis
    # # coordinates are created always from the parameters in the core
    # # definitions of IDPConfGen
    dummy_CA_m1_coord = np.array((0.0, 1.0, 0.0))
    n_terminal_N_coord = np.array((0.0, 0.0, 0.0))
    n_terminal_CA_coord = np.array((distances_N_CA[input_seq[0]], 0.0, 0.0))

    # seed coordinates array
    seed_coords = np.array((
        dummy_CA_m1_coord,
        n_terminal_N_coord,
        n_terminal_CA_coord,
        ))

    ss_masks = create_sidechains_masks_per_residue(
        residue_numbers,
        atom_labels,
        backbone_atoms,
        )

    bbi0_register = []
    bbi0_R_APPEND = bbi0_register.append
    bbi0_R_POP = bbi0_register.pop
    bbi0_R_CLEAR = bbi0_register.clear

    COi0_register = []
    COi0_R_APPEND = COi0_register.append
    COi0_R_POP = COi0_register.pop
    COi0_R_CLEAR = COi0_register.clear

    NHi0_register = []
    NHi0_R_APPEND = NHi0_register.append
    NHi0_R_POP = NHi0_register.pop
    NHi0_R_CLEAR = NHi0_register.clear

    res_R = []  # residue number register
    res_R_APPEND = res_R.append
    res_R_POP = res_R.pop
    res_R_CLEAR = res_R.clear

    broke_on_start_attempt = False
    start_attempts = 0
    max_start_attempts = 50  # maximum attempts to start a conformer
    # because we are building from a experimental database there can be
    # some angle combinations that fail on our validation process from start
    # if this happens more than `max_start_attemps` the production is canceled.

    # STARTS BUILDING
    # aligns the indexes so that conformers can be named properly in
    # multicore operations
    start_conf = nconfs * execution_run
    end_conf = start_conf + nconfs
    conf_n = start_conf
    while conf_n < end_conf:
        print('building conformer: ', conf_n)

        # prepares cycles for building process
        # cycles need to be regenerated every conformer because the first
        # atom build is a C and the last atom built is the CA, which breaks
        # periodicity
        # all these are dictionaries
        bond_lens = cycle((
            distances_CA_C,  # used for PHI
            distances_C_Np1,  # used for PSI
            distances_N_CA,  # used for OMEGA
            ))
        bond_bend = cycle((
            build_bend_angles_N_CA_C,  # used for PHI
            build_bend_angles_CA_C_Np1,  # used for PSI
            build_bend_angles_Cm1_N_CA,  # used for OMEGA
            ))

        # in the first run of the loop this is unnecessary, but is better to
        # just do it once than flag it the whole time
        coords[:, :] = NAN
        bb[:, :] = NAN
        bb_CO[:, :] = NAN
        bb_NH[:, :] = NAN
        for _mask, _coords in ss_masks:
            _coords[:, :] = NAN

        bb[:3, :] = seed_coords  # this contains a dummy coord at position 0

        # IMPLEMENT SIDECHAINS HERE

        bbi = 2  # starts at 2 because the first 3 atoms are already placed
        bbi0_R_CLEAR()
        bbi0_R_APPEND(bbi)

        COi = 0  # carbonyl atoms
        COi0_R_CLEAR()
        COi0_R_APPEND(COi)

        # NHi is 0 because
        # the first residue, remember the first residue as 'H1,2,3'
        # see bb_NH_builder definition
        NHi = 0
        NHi0_R_CLEAR()
        NHi0_R_APPEND(NHi)

        # residue integer number
        current_res_number = 0
        res_R_CLEAR()
        res_R_APPEND(current_res_number)

        backbone_done = False
        number_of_trials = 0
        # TODO: use or not to use number_of_trials2? To evaluate in future.
        #number_of_trials2 = 0
        # run this loop until a specific BREAK is triggered
        while 1:  # 1 is faster than True :-)

            # I decided to use an if-statement here instead of polymorph
            # the else clause to a `generative_function` variable because
            # the resulting overhead from the extra function call and
            # **kwargs handling was greater then the if-statement processing
            # https://pythonicthoughtssnippets.github.io/2020/10/21/PTS14-quick-in-if-vs-polymorphism.html
            if generative_function:
                agls = generative_function(
                    nres=RINT(1, 6),
                    cres=(bbi - 2) // 3,
                    )

            else:
                # following `aligndb` function,
                # `angls` will always be cyclic with:
                # phi - psi - omega - phi - psi - omega - (...)
                agls = angles[RC(slices), :].ravel()

            # index at the start of the current cycle
            try:
                for torsion in agls:
                    # bbi -2 makes the match, because despite the first atom
                    # being placed is a C, it is placed using the PHI angle of a
                    # residue. And PHI is the first angle of that residue angle
                    # set.
                    current_res_number = (bbi - 2) // 3
                    current_residue = input_seq[current_res_number]

                    # bbi is the same for bb_real and bb, but bb_real indexing
                    # is displaced by 1 unit from bb. So 2 in bb_read is the
                    # next atom of 2 in bb.
                    bb_real[bbi, :] = MAKE_COORD_Q_LOCAL(
                        bb[bbi - 2, :],
                        bb[bbi - 1, :],
                        bb[bbi, :],
                        next(bond_lens)[current_residue],
                        next(bond_bend)[current_residue],
                        torsion,
                        )
                    bbi += 1

            except IndexError:
                # IndexError happens when the backbone is complete
                # in this protocol the last atom build was a carbonyl C
                # bbi is the last index of bb + 1, and the last index of
                # bb_real + 2

                # activate flag to finish loop at the end
                backbone_done = True

                # add the carboxyls
                coords[[OXT1_index, OXT_index]] = MAKE_COORD_Q_COO_LOCAL(
                    bb[-2, :],
                    bb[-1, :],
                    )

            # builds carbonyl atoms. Two situations can happen here:
            # 1) the backbone is not complete - the last atom is CA
            # 2) the backbone is complete - the last atom is C
            # whichever the case, with `bbi - 1` applying on the range bellow
            # will build carbonyls for the new chain expect for the last C
            # if the chain is completed.
            # this is so because range(X, Y, Z) equals range(X, Y-1, Z)
            # if Y-1 is not in range(X, Y, Z). And, this is the case for N-CA
            # pair.
            for k in range(bbi0_register[-1], bbi - 1, 3):
                bb_CO[COi, :] = MAKE_COORD_Q_CO_LOCAL(
                    bb_real[k - 1, :],
                    bb_real[k, :],
                    bb_real[k + 1, :],
                    )
                COi += 1

            # Adds N-H Hydrogens
            for k in range(bbi0_register[-1] + 1, bbi, 3):

                if residue_labels_bb_simulating[k] == 'PRO':
                    continue

                # MAKE_COORD_Q_CO_LOCAL can be used for NH by giving
                # disntace and bend parameters
                bb_NH[NHi, :] = MAKE_COORD_Q_CO_LOCAL(
                    bb_real[k - 1, :],
                    bb_real[k, :],
                    bb_real[k + 1, :],
                    distance=DISTANCE_NH,
                    bend=BUILD_BEND_H_N_C,
                    )
                NHi += 1

            # Adds sidechain template structures
            if with_sidechains:
                for res_i in range(res_R[-1], current_res_number + backbone_done):  # noqa: E501

                    _sstemplate, _sidechain_idxs = \
                        SIDECHAIN_TEMPLATES[input_seq_3_letters[res_i]]  # 1 is faster than True :-)

                    sscoords = PLACE_SIDECHAIN_TEMPLATE(
                        bb_real[res_i * 3:res_i * 3 + 3, :],  # from N to C
                        _sstemplate,
                        )

                    ss_masks[res_i][1][:, :] = sscoords[_sidechain_idxs]

                # Transfers coords to the main coord array
                for _smask, _sidecoords in ss_masks[: current_res_number + backbone_done]:  # noqa: E501
                    coords[_smask] = _sidecoords

            coords[bb_mask] = bb_real  # do not consider the initial dummy atom
            coords[carbonyl_mask] = bb_CO
            coords[NHydrogen_mask] = bb_NH

            distances_ij = CALC_DISTS(coords)

            energy = njit_calc_LJ_energy(
                acoeff,
                bcoeff,
                distances_ij,
                bonds_ge_3_mask,
                )

            energy += NANSUM(charges_ij / distances_ij)
            #print(energy)

            if energy > 0:
                # reset coordinates to the original value
                # before the last chunk added

                # reset the same chunk maximum 5 times,
                # after that reset also the chunk before
                if number_of_trials > 5:
                    bbi0_R_POP()
                    COi0_R_POP()
                    NHi0_R_POP()
                    res_R_POP()
                    number_of_trials = 0
                    #number_of_trials2 += 1

                #if number_of_trials2 > 50:
                #    bbi0_R_POP()
                #    COi0_R_POP()
                #    NHi0_R_POP()
                #    res_R_POP()
                #    number_of_trials2 = 0

                try:
                    _bbi0 = bbi0_register[-1]
                    _COi0 = COi0_register[-1]
                    _NHi0 = NHi0_register[-1]
                    _resi0 = res_R[-1]
                except IndexError:
                    # if this point is reached,
                    # we erased until the beginning of the conformer
                    # discard conformer, something went really wrong
                    broke_on_start_attempt = True
                    break  # while loop

                # clean previously built protein chunk
                bb_real[_bbi0:bbi, :] = NAN
                bb_CO[_COi0:COi, :] = NAN
                bb_NH[_NHi0:NHi, :] = NAN

                # do the same for sidechains
                # ... it is not necessary because in the loop above
                # sidechain coordinates are replaces when repositioned

                # reset also indexes
                bbi = _bbi0
                COi = _COi0
                NHi = _NHi0
                current_res_number = _resi0

                # coords needs to be reset because size of protein next
                # chunks may not be equal
                coords[:, :] = NAN

                # prepares cycles for building process
                # this is required because the last chunk created may have been
                # the final part of the conformer
                if backbone_done:
                    bond_lens = cycle((
                        distances_CA_C,
                        distances_C_Np1,
                        distances_N_CA,
                        ))
                    bond_bend = cycle((
                        build_bend_angles_N_CA_C,
                        build_bend_angles_CA_C_Np1,
                        build_bend_angles_Cm1_N_CA,
                        ))

                # we do not know if the next chunk will finish the protein
                # or not
                backbone_done = False
                number_of_trials += 1
                continue  # send back to the while loop

            # if the conformer is valid
            number_of_trials = 0
            bbi0_R_APPEND(bbi)
            COi0_R_APPEND(COi)
            NHi0_R_APPEND(NHi)
            # the residue where the build process stopped
            res_R_APPEND(current_res_number)

            if backbone_done:
                # this point guarantees all protein atoms are built
                break  # while loop
        # END of while loop

        if broke_on_start_attempt:
            start_attempts += 1
            if start_attempts > max_start_attempts:
                log.error(
                    'Reached maximum amount of starts. Canceling... '
                    f'{end_conf - conf_n} of {end_conf - start_conf} '
                    'conformers were not built.'
                    )
                return
            broke_on_start_attempt = False
            continue

        # until sidechains are implemented this is needed
        # TODO: need to implement the N-terminal H1-3 atoms
        # lun 21 dic 2020 18:14:05 EST
        relevant = np.logical_not(np.isnan(coords[:, 0]))

        pdb_string = gen_PDB_from_conformer(
            input_seq,
            atom_labels[relevant],
            residue_numbers[relevant],
            ROUND(coords[relevant], decimals=3),
            )

        fname = f'{conformer_name}_{conf_n}.pdb'
        with open(fname, 'w') as fout:
            fout.write(pdb_string)

        conf_n += 1

    return


def gen_PDB_from_conformer(
        input_seq,
        atom_labels,
        residues,
        coords,
        ALF=atom_line_formatter,
        AA1TO3=aa1to3,
        ROUND=np.round,
        ):
    """."""
    lines = []
    LINES_APPEND = lines.append
    ALF_FORMAT = ALF.format
    resi = -1

    # this is possible ONLY because there are no DOUBLE CHARS atoms
    # in the atoms that constitute a protein chain
    ATOM_LABEL_FMT = ' {: <3}'.format

    for i in range(len(atom_labels)):

        if atom_labels[i] == 'N':
            resi += 1
            current_residue = input_seq[resi]
            current_resnum = residues[i]

        atm = atom_labels[i].strip()
        ele = atm.lstrip('123')[0]

        if len(atm) < 4:
            atm = ATOM_LABEL_FMT(atm)

        LINES_APPEND(ALF_FORMAT(
            'ATOM',
            i,
            atm,
            '',
            AA1TO3[current_residue],
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

    return '\n'.join(lines)


def generate_atom_labels(input_seq, atom_labels_dictionary):
    """."""
    labels = []
    LE = labels.extend

    first_residue_atoms = atom_labels_dictionary[input_seq[0]]

    # the first residue is a special case, we add here the three protons
    # for consistency with the forcefield
    # TODO: parametrize? multiple forcefields?
    for atom in first_residue_atoms:
        if atom == 'H':
            LE(('H1', 'H2', 'H3'))
        else:
            labels.append(atom)

    for residue in input_seq[1:]:
        LE(atom_labels_dictionary[residue])

    labels.append('OXT')

    assert Counter(labels)['N'] == len(input_seq)
    assert labels[-1] == 'OXT'
    assert 'H1' in labels
    assert 'H2' in labels
    assert 'H3' in labels
    return labels


def generate_residue_numbers(atom_labels, start=1):
    """
    Create a list of residue numbers based on atom labels.

    Considers `N` to be the first atom of the residue.
    If this is not the case, the output can be meaningless.

    Returns
    -------
    list of ints
    """
    if atom_labels[0] == 'N':
        start -= 1  # to compensate the +=1 implementation in the for loop

    residues = []
    RA = residues.append

    for al in atom_labels:
        if al == 'N':
            start += 1
        RA(start)

    return residues


def generate_residue_labels(
        input_seq,
        atom_labels,
        aa1to3=aa1to3,
        ):
    """Generate residue 3-letter labels."""
    assert atom_labels[0] == 'N'
    count_labels = Counter(atom_labels)['N']
    assert count_labels == len(input_seq), (count_labels, len(input_seq))

    residue_labels = []
    RL_APPEND = residue_labels.append

    counter = -1
    for atom in atom_labels:
        if atom == 'N':
            counter += 1
        RL_APPEND(aa1to3[input_seq[counter]])

    assert len(residue_labels) == len(atom_labels), 'Lengths differ.. strangely'
    return residue_labels



# TODO: abstract out this function
def populate_ff_parameters_in_structure(
        atom_labels,
        residue_numbers,
        residue_labels,
        force_field,
        ):
    """
    Creates a list of SIGMAS, EPSILONS, and CHARGES.

    Matches label information in the conformer to the respective values
    in the force field.
    """
    assert len(atom_labels) == len(residue_numbers) == len(residue_labels)

    sigmas_l = []
    epsilons_l = []
    charges_l = []

    sigmas_append = sigmas_l.append
    epsilons_append = epsilons_l.append
    charges_append = charges_l.append

    zipit = zip(atom_labels, residue_numbers, residue_labels)
    for atom_name, res_num, res_label in zipit:

        # adds C to the terminal residues
        if res_num == residue_numbers[-1]:
            res = 'C' + res_label
            was_in_C_terminal = True
            assert res.isupper() and len(res) == 4, res

        elif res_num == residue_numbers[0]:
            res = 'N' + res_label
            was_in_N_terminal = True
            assert res.isupper() and len(res) == 4, res

        else:
            res = res_label

        # TODO:
        # define protonation state in parameters
        if res_label.endswith('HIS'):
            res_label = res_label[:-3] + 'HIP'

        try:
            # force field atom type
            atype = force_field[res][atom_name]['type']

        # TODO:
        # try/catch is here to avoid problems with His...
        # for this purpose we are only using side-chains
        except KeyError:
            raise KeyError(tuple(force_field[res].keys()))

        ep = float(force_field[atype]['epsilon'])
        # epsilons by order of atom
        epsilons_append(ep)

        sig = float(force_field[atype]['sigma']) * 10  # conversion to Angstroms
        # sigmas by order of atom
        sigmas_append(sig)

        charge = float(force_field[res][atom_name]['charge'])
        charges_append(charge)

    assert len(epsilons_l) == len(sigmas_l) == len(charges_l), 'Lengths differ.'
    assert len(epsilons_l) == len(atom_labels), 'Lengths differ. There should be one epsilon per atom.'
    assert was_in_C_terminal, 'The C terminal residue was never computed. It should have been computed'
    assert was_in_N_terminal, 'The N terminal residue was never computed. It should have been computed'

    # brief theoretical review:
    # sigmas and epsilons are combined parameters self vs self (i vs i)
    # charges refer only to the self alone (i)
    sigmas_ii = np.array(sigmas_l)
    epsilons_ii = np.array(epsilons_l)
    charges_i = np.array(charges_l)

    assert epsilons_ii.shape == (len(epsilons_l),)
    assert sigmas_ii.shape == (len(sigmas_l),)
    assert charges_i.shape == (len(charges_l),)

    return sigmas_ii, epsilons_ii, charges_i


def gen_ij_pairs_upper_diagonal(data):
    """Generator"""
    for i in range(len(data) - 1):
        for j in range(i + 1, len(data)):
            yield (data[i], data[j])


def gen_atom_pair_connectivity_masks(
        res_names_ij,
        res_num_ij,
        atom_names_ij,
        connectivity_intra,
        connectivity_inter,
        ):
    """
    To generate bonds masks we need the residue numbers and the atom names.

    Depends
    -------
    `are_connected`
    """
    zipit = zip(res_names_ij, res_num_ij, atom_names_ij)
    counter = 0
    for (rn1, _), (n1, n2), (a1, a2) in zipit:

        found_connectivity = are_connected(
            int(n1),
            int(n2),
            rn1,
            a1,
            a2,
            connectivity_intra,
            connectivity_inter,
            )

        if found_connectivity:
            yield counter

        counter += 1


def are_connected(n1, n2, rn1, a1, a2, bonds_intra, bonds_inter):
    """
    Detect if a certain atom pair is bonded accordind to criteria.

    Considers only to the self residue and next residue
    """
    # requires
    assert isinstance(n1, int) and isinstance(n2, int), (type(n1), type(n2))
    assert all(isinstance(i, str) for i in (rn1, a1, a2)), \
        (type(i) for i in (rn1, a1, a2))
    assert all(isinstance(i, dict) for i in (bonds_intra, bonds_inter)), \
        (type(i) for i in (bonds_intra, bonds_inter))

    answer = (
        (n1 == n2 and a2 in bonds_intra[rn1][a1])
        or (
            n1 + 1 == n2
            and (
                a1 in bonds_inter  # this void KeyError
                and a2 in bonds_inter[a1]
                )
            )
        )

    assert isinstance(answer, bool)
    return answer


def calc_outer_multiplication_upper_diagonal_raw(data, result):
    """Calculate the upper diagonal multiplication with for loops."""
    #assert data.ndim == 1, \
    #    'Data has {data.ndim} dimensions, should have only 1.'

    #result = np.empty(data.size - 1, np.float64)
    c = 0
    len_ = len(data)
    for i in range(len_ - 1):
        for j in range(i + 1, len_):
            result[c] = data[i] * data[j]
            c += 1

    assert result.size == (data.size * data.size - data.size) // 2
    assert abs(result[0] - data[0] * data[1]) < 0.0000001
    assert abs(result[-1] - data[-2] * data[-1]) < 0.0000001
    #return result

njit_calc_multiplication_upper_diagonal_raw = \
    njit(calc_outer_multiplication_upper_diagonal_raw)


def calc_sum_upper_diagonal_raw(data, result):
    """
    Calculate outer sum for upper diagonal with for loops.

    Avoids creating temporary hiper large arrays
    """
    #assert data.ndim == 1, \
    #    'Data has {data.ndim} dimensions, should have only 1.'

    #result = np.empty(data.size - 1, np.float64)
    c = 0
    len_ = len(data)
    for i in range(len_ - 1):
        for j in range(i + 1, len_):
            result[c] = data[i] + data[j]
            c += 1

    assert result.size == (data.size * data.size - data.size) // 2
    assert abs(result[0] - (data[0] + data[1])) < 0.0000001
    assert abs(result[-1] - (data[-2] + data[-1])) < 0.0000001
    #return result


njit_calc_sum_upper_diagonal_raw = njit(calc_sum_upper_diagonal_raw)


def calc_LJ_energy(acoeff, bcoeff, dists_ij, to_eval_mask, NANSUM=np.nansum):
    """Calculates Lennard-Jones Energy."""
    #assert dists_ij.size == acoeff.size == bcoeff.size, (dists_ij.size, acoeff.size)
    #assert dists_ij.size == to_eval_mask.size

    ar = acoeff / (dists_ij ** 12)
    br = bcoeff / (dists_ij ** 6)
    energy_ij = ar - br

    # nansum is used because some values in dists_ij are expected to be nan
    #return energy_ij, NANSUM(energy_ij[to_eval_mask])
    return NANSUM(energy_ij[to_eval_mask])

njit_calc_LJ_energy = njit(calc_LJ_energy)


def create_sidechains_masks_per_residue(
        residue_numbers,
        atom_labels,
        backbone_atoms,
        ):
    """
    Create a map of numeric indexing masks pointing to side chains atoms.

    Create separate masks per residue.

    Parameters
    ----------
    residue_numbers : np.ndarray, shape (N,)
        The atom residue numbers of the protein.

    atom_labels : np.ndarray, shape (N,)
        The atom labels of the protein.

    backbone_atoms : list or tuple
        The labels of all possible backbone atoms.

    Returns
    -------
    list of tuples of length 2
        List indexes refer to protein residues, index 0 is residue 1.
        Per residue, a tuple of length 2 is given. Tuple index 0 are
        the indexes of that residue sidechain atoms mapped to an array
        of the `atom_labels` and `residue_numbers` characteristics.
        The tuple index 1 is an array of length M, where M is the number
        of sidechain atoms for that residue, defaults to np.nan.
    """
    assert residue_numbers.size == atom_labels.size
    assert type(backbone_atoms) in (list, tuple, np.ndarray)

    ss = []
    ssa = ss.append

    is_backbone = np.isin(atom_labels, backbone_atoms)
    is_sidechain = np.logical_not(is_backbone)

    float64 = np.float64
    full = np.full
    logical_and = np.logical_and
    nan = np.nan
    where = np.where

    # sorted here is mandatory because ss indexes must follow the residue
    # numbering
    for resnum in sorted(list(set(residue_numbers))):
        is_residue = residue_numbers == resnum
        bool_result = logical_and(is_residue, is_sidechain)
        ssa((
            where(bool_result)[0],
            full((np.sum(bool_result), 3), nan, dtype=float64)
            ))

    return ss


def create_bonds_apart_mask_for_ij_pairs(
        atom_labels,
        residue_numbers,
        residue_labels,
        bonds_intra,
        bonds_inter,
        base_bool=False,
        ):
    """."""
    atom_labels_ij_gen = gen_ij_pairs_upper_diagonal(atom_labels)
    residue_numbers_ij_gen = gen_ij_pairs_upper_diagonal(residue_numbers)
    residue_labels_ij_gen = gen_ij_pairs_upper_diagonal(residue_labels)

    bonds_indexes_gen = gen_atom_pair_connectivity_masks(
        residue_labels_ij_gen,
        residue_numbers_ij_gen,
        atom_labels_ij_gen,
        bonds_intra,
        bonds_inter,
        )

    num_ij_pairs = len(atom_labels) * (len(atom_labels) - 1) // 2
    other_bool = not base_bool
    bonds_mask = np.full(num_ij_pairs, base_bool)

    for idx in bonds_indexes_gen:
        bonds_mask[idx] = other_bool

    return bonds_mask


# NOT USED HERE
# kept to maintain compatibility with cli_validate.
# TODO: revisit cli_validate
def generate_vdW_data(
        atom_labels,
        residue_numbers,
        residue_labels,
        vdW_radii,
        bonds_apart=3,
        tolerance=0.4,
        ):
    # {{{
    """
    Generate van der Waals related data structures.

    Generated data structures are aligned to the output generated by
    scipy.spatial.distances.pdist used to compute all pairs distances.

    This function generates an (N,) array with the vdW thresholds aligned
    with the `pdist` output. And, a (N,) boolean array where `False`
    denotes pairs of covalently bonded atoms.

    Input should satisfy:

    >>> len(atom_labels) == len(residue_numbers) == len(residue_labels)

    Parameters
    ----------
    atom_labels : array or list
        The ordered list of atom labels of the target conformer upon
        which the returned values of this function will be used.

    residue_numbers : array or list
        The list of residue numbers corresponding to the atoms of the
        target conformer.

    residue_labels : array or list
        The list of residue 3-letter labels of each atom of the target
        conformer.

    vdW_radii : dict
        A dictionary containing the van der Waals radii for each atom
        type (element) present in `atom_labels`.

    bonds_apart : int
        The number of bonds apart to ignore vdW clash validation.
        For example, 3 means vdW validation will only be computed for
        atoms at least 4 bonds apart.

    tolerance : float
        The tolerance in Angstroms.

    Returns
    -------
    nd.array, dtype=np.float
        The vdW atom pairs thresholds. Equals to the sum of atom vdW
        radius.

    nd.array, dtype=boolean
        `True` where the distance between pairs must be considered.
        `False` where pairs are covalently bound.

    Note
    ----
    This function is slow, in the order of 1 or 2 seconds, but it is
    executed only once at the beginning of the building protocol.
    """
    # }}}
    print('here')
    assert len(atom_labels) == len(residue_numbers)
    assert len(atom_labels) == len(residue_labels)

    # we need to use re because hydrogen atoms start with integers some times
    atoms_char = re.compile(r'[CHNOPS]')
    findall = atoms_char.findall
    cov_topologies = generate_residue_template_topology()
    bond_structure_local = \
        expand_topology_bonds_apart(cov_topologies, bonds_apart)
    inter_connect_local = inter_residue_connectivities[bonds_apart]

    # adds OXT to the bonds connectivity, so it is included in the
    # final for loop creating masks
    add_OXT_to_residue(bond_structure_local[residue_labels[-1]])

    # the following arrays/lists prepare the conformer label data in agreement
    # with the function scipy.spatial.distances.pdist, in order for masks to be
    # used properly
    #
    # creates atom label pairs
    atoms = np.array([
        (a, b)
        for i, a in enumerate(atom_labels, start=1)
        for b in atom_labels[i:]
        ])

    # creates vdW radii pairs according to atom labels
    vdw_pairs = np.array([
        (vdW_radii[findall(a)[0]], vdW_radii[findall(b)[0]])
        for a, b in atoms
        ])

    # creats the vdW threshold according to vdW pairs
    # the threshold is the sum of the vdW radii pair
    vdW_sums = np.power(np.sum(vdw_pairs, axis=1) - tolerance, 2)

    # creates pairs for residue numbers, so we know from which residue number
    # is each atom of the confronted pair
    res_nums_pairs = np.array([
        (a, b)
        for i, a in enumerate(residue_numbers, start=1)
        for b in residue_numbers[i:]
        ])

    # does the same as above but for residue 3-letter labels
    res_labels_pairs = (
        (_a, _b)
        for i, _a in enumerate(residue_labels, start=1)
        for _b in residue_labels[i:]
        )

    # we want to compute clashes to all atoms that are not covalentely bond
    # that that have not certain connectivity between them
    # first we assum True to all cases, and we replace to False where we find
    # a bond connection
    vdW_non_bond = np.full(len(atoms), True)

    counter = -1
    zipit = zip(atoms, res_nums_pairs, res_labels_pairs)
    for (a1, a2), (res1, res2), (l1, _) in zipit:
        counter += 1
        if res1 == res2 and a2 in bond_structure_local[l1][a1]:
            # here we found a bond connection
            vdW_non_bond[counter] = False

        # handling the special case of a covalent bond between two atoms
        # not part of the same residue
        elif a1 in inter_connect_local \
                and res2 == res1 + 1 \
                and a2 in inter_connect_local[a1]:
            vdW_non_bond[counter] = False

    print('done')
    return vdW_sums, vdW_non_bond


if __name__ == "__main__":
    libcli.maincli(ap, main)
