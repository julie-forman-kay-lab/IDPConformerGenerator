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
    extract_interpairs_from_db,
    extract_intrapairs_from_db,
    )
from idpconfgen.libs.libmulticore import pool_function
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
            input_seq = list(input_seq.values())[0]
            seq_len = len(input_seq)
            log.info(S(f'input sequence: {input_seq}'))
    else:
        seq_len = len(input_seq)
        log.info(S(f'input sequence: {input_seq}'))
    
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
    log.info(T("Calculating all-contacts probability matrix"))
    matrix_consume = partial(contact_matrix, sequence=input_seq)
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
    contact_mtx = np.zeros((seq_len, seq_len))
    location_mtx = {}
    for id, hit, loc in matrix_execute_pool:
        contact_mtx = np.add(contact_mtx, hit)
        location_mtx[id] = loc

    norm_contact_mtx = (contact_mtx - np.min(contact_mtx)) / (np.max(contact_mtx) - np.min(contact_mtx))  # noqa: E501
    log.info(S('done'))
    
    log.info(T('saving heatmap distribution plot to disk'))
    if not (plot_name.endswith(".png") or plot_name.endswith(".jpg")):
        plot_path = Path(plot_name)
        plot_name = plot_path.with_suffix('.png')
    plot_contacts_matrix(norm_contact_mtx, input_seq, plot_name)
    log.info(S('done'))


if __name__ == "__main__":
    libcli.maincli(ap, main)
