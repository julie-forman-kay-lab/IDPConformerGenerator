"""# noqa: D205, D400, E501
Extracts secondary structure information from PDBs.

This client is upgraded to DSSP-PPII according to:

   Mansiaux Y, Joseph AP, Gelly J-C, de Brevern AG (2011) Assignment of
   PolyProline II conformation and analysis od sequence-structure
   relationship. PLoS One 6(3): e18401. doi:10.1371/journal.pone.0018401

Uses an external third party software.

USAGE:
    $ idpconfgen sscalc [PDBS]
    $ idpconfgen sscalc [PDBS] -o mysscalc.json -d splitsscalc.tar -cmd <DSSP EXEC> --plot -n
"""
import argparse
import os
import shutil
import traceback
from functools import partial

import numpy as np

from idpconfgen import Path, log
from idpconfgen.cli_dssppii import dssp_ppii_assignment
from idpconfgen.components.plots.plotfuncs import plot_fracSS
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_dictionary_from_disk,
    read_path_bundle,
    save_dictionary,
    save_pairs_to_disk,
    )
from idpconfgen.libs.libmulticore import (
    consume_iterable_in_list,
    pool_function_in_chunks,
    )
from idpconfgen.libs.libparse import split_pdb_by_dssp
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.idpconfgen_sscalc'
TMPDIR = '__tmpsscalc__'

_name = 'sscalc'
_help = 'Calculate secondary structure profile.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-cmd',
    '--dssp_cmd',
    help="Shell command for the DSSP executable. Defaults to `mkdssp`.",
    default="mkdssp",
    )

ap.add_argument(
    '-o',
    '--output',
    help=(
        "A path to the file where the PDBID secondary structure and FASTA"
        " information dictionary will be saved. "
        "Defaults to sscalc.json, requires \'.json\' extension."
        ),
    type=Path,
    default='sscalc.json',
    action=libcli.CheckExt({'.json'}),
    )

# can't use the libcli option because this one is different
ap.add_argument(
    '-d',
    '--destination',
    help=(
        'Destination folder where the split PDB files will be stored. '
        'Defaults to current working directory.'
        'Alternatively, you can provide a path to a .tar file '
        'where PDBs will be saved.'
        ),
    type=Path,
    default='sscalc_splitted.tar',
    )

ap.add_argument(
    '-u',
    '--update',
    help='Updates a previous SSCALC JSON file.',
    type=Path,
    default=None,
    )

ap.add_argument(
    '-tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )

libcli.add_argument_reduced(ap)
libcli.add_argument_minimum(ap)
libcli.add_argument_chunks(ap)
libcli.add_argument_ncores(ap)
libcli.add_argument_plot(ap)


def dssppi_helper(pdb_file, dssp_cmd, **kwargs):
    """."""
    pf = pdb_file.resolve()
    result = str.encode(os.linesep.join(dssp_ppii_assignment(str(pf), dssp_cmd)))  # noqa: E501
    yield from split_pdb_by_dssp(pf, result, **kwargs)


def main(
        pdb_files,
        dssp_cmd="dssp",
        chunks=1000,
        destination='sscalc_splitted.tar',
        func=None,  # here just to receive from main cli.py
        minimum=3,
        ncores=1,
        output='sscalc_output.json',
        reduced=False,
        tmpdir=TMPDIR,
        update=None,
        plot=False,
        plotvars=None,
        ):
    """# noqa: D205, D400, E501
    Run main cli logic.

    Parameters
    ----------
    ss_cmd : str or Path
        The command to run with subprocess module.

    pdb_files : str or Path, required
        Location for PDB files to operate on, can be within a folder or inside .TAR file

    output : string or Path, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.

    chunks : int, optional
        The number of items to process in memory before saving the
        results from the disk.

    reduced : Bool, optional
        Whether to reduce secondary structure information to L/H/E

    ncores : int
        The numbers of cores to use.

    plot : Bool, optional
        Whether to plot the fractional secondary structure information

    plotvars : dictionary, optional
        Parameters for creating the plot
    """
    log.info(T('Extracting Secondary structure information'))
    init_files(log, LOGFILESNAME)

    # update is inspected here because it can raise errors and it is better
    # that such errors are spotted before performing the expensive calculations
    if update:
        previous = read_dictionary_from_disk(update)

    # TODO: REFACTOR THIS BLOCK
    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))

    log.info(T('preparing task execution'))

    try:
        consume_func = partial(
            consume_iterable_in_list,
            dssppi_helper,
            dssp_cmd=dssp_cmd,
            reduced=reduced,
            minimum=minimum,
            )

        execute = partial(
            report_on_crash,
            consume_func,
            ROC_exception=Exception,
            ROC_prefix=_name,
            )

        # generator
        execute_pool = pool_function_in_chunks(
            execute,
            pdbs2operate,  # items to process
            ncores=ncores,
            chunks=chunks,
            )
        # this implementation is very specific for this case
        # this is why I haven spread it into functions for now
        dssp_data = {}  # stores DSSP data to save at the end
        pdb_data = {}  # temporarily stores data to be writen to the disk
        for chunk in execute_pool:
            for result in chunk:
                for fname, dsspdict, pdb_split in result:
                    assert fname not in dssp_data
                    dssp_data[fname] = dsspdict

                    # notice the copy, this is needed for the .clear()
                    # to work later on
                    pdb_data[f'{fname}.pdb'] = pdb_split
            save_pairs_to_disk(pdb_data.items(), destination=destination)
            pdb_data.clear()  # clears the dictionary to release memory

    except Exception as err:
        log.error('FAILED')
        log.debug(traceback.format_exc())
        raise err

    finally:
        if _istarfile:
            shutil.rmtree(tmpdir)

    if update:
        save_dictionary(previous.update(dssp_data), output)
    else:
        save_dictionary(dssp_data, output)

    if plot:
        log.info(T("Plotting fractional secondary structure information"))
        log.info(S("Tip: plotting is best for same protein-system ensembles."))

        max_residues = 0
        for key in dssp_data:
            if len(dssp_data[key]["dssp"]) > max_residues:
                max_residues = len(dssp_data[key]["dssp"])
        n_confs = len(dssp_data)

        plotvars = plotvars or dict()
        plt_default = {
            'sstype': 'DSSP',
            'filename': 'plot_dssp_fracSS.png',
            }
        plt_default.update(plotvars)

        if reduced:
            frac_dssp = {
                'L': np.zeros(max_residues),
                'H': np.zeros(max_residues),
                'E': np.zeros(max_residues),
                }
            p = np.ones(n_confs)
            p = p / len(p)

            c = 0
            for conf in dssp_data:
                r = 0
                for ss in dssp_data[conf]["dssp"]:
                    frac_dssp[ss][r] += p[c]
                    r += 1
                c += 1
        else:
            frac_dssp = {
                'H': np.zeros(max_residues),
                'B': np.zeros(max_residues),
                'E': np.zeros(max_residues),
                'G': np.zeros(max_residues),
                'I': np.zeros(max_residues),
                'T': np.zeros(max_residues),
                'S': np.zeros(max_residues),
                'P': np.zeros(max_residues),
                ' ': np.zeros(max_residues),
                }
            p = np.ones(n_confs)
            p = p / len(p)

            c = 0
            for conf in dssp_data:
                r = 0
                for ss in dssp_data[conf]["dssp"]:
                    frac_dssp[ss][r] += p[c]
                    r += 1
                c += 1

        errs = plot_fracSS(max_residues, frac_dssp, **plt_default)
        for e in errs:
            log.info(S(f'{e}'))

        log.info(S(f'saved plot: {plt_default["filename"]}'))

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
