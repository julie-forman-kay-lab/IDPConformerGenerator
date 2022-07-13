"""
Extracts backbone torsion angles from PDBs.

PROTOCOL:

1. Reads backbone coordinates (N, CA, C) from PDB/mmCIF files.
2. Calculates torsion angles from the backbone.
3. Saves results to a JSON dictionary where keys are the input file
    names and the value is a dictionary containing three lists: 'OMEGA',
    'PHI', and 'PSI'.
4. If `source` JSON file is given, updates that file with the new
    information. Preexisting keys are deleted.

CONTROLLED CHECKS:

For each PDB/mmCIF analyzed, fails if:

1. The number of N, CA, and C atoms differ.
2. There are unexpected inconsistencies in the PDB/mmCIF files.
3. Any of the consecutive atoms are more than 2.1A apart.

Failed PDBs are registered in `.rpr_on_crash` files and ignored.

PLOTTING:

The default torsion distribution plotted are "Phi" angles. To change this to
Omega or Psi, simply `--plots angtype=omega` or `--plots angtype=psi`.

USAGE:
    $ idpconfgen torsions [PDBS]
    $ idpconfgen torsions [PDBS] -sc file.json
    $ idpconfgen torsions [PDBS] -sc file.json -n
    $ idpconfgen torsions [PDBS] -sc file.json -o mytorsions.json -n
    $ idpconfgen torsions [PDBS] -sc file.json -o mytorsions.json -n -deg
    $ idpconfgen torsions [PDBS] -sc file.json -o mytorsions.json -n -deg --plot
"""
import argparse
from functools import partial
from math import radians

import numpy as np

from idpconfgen import Path, log
from idpconfgen.components.plots.plotfuncs import plot_fracSS, plot_torsions
from idpconfgen.libs import libcli
from idpconfgen.libs.libhigherlevel import cli_helper_calc_torsions
from idpconfgen.libs.libio import (
    FileReaderIterator,
    read_dictionary_from_disk,
    save_dict_to_json,
    )
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.libs.libparse import pop_difference_with_log, values_to_dict
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.idpconfgen_torsion'

_name = 'torsions'
_help = 'Calculate torsion angles for PDB files.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

libcli.add_argument_pdb_files(ap)
libcli.add_argument_source(ap)
libcli.add_argument_output(ap)
libcli.add_argument_degrees(ap)
libcli.add_argument_ncores(ap)
libcli.add_argument_plot(ap)

ap.add_argument(
    '--ramaplot',
    help=(
        "Change default parameters of fractional secondary structure plot "
        "using alpha, beta, and other regions of the ramachandran space. "
        "Optionally to be used with ``--plot``."
        "Example: --ramaplot filename=fracRama.png colors=r,g,b"
        ),
    nargs='*',
    )


def main(
        pdb_files,
        source=None,
        output=None,
        degrees=False,  # documentation for `add_argument_degrees()` specifies rad as default # noqa: E501
        func=None,
        ncores=1,
        plot=False,
        plotvars=None,
        ramaplot=None,
        ):
    """# noqa: D205, D400, E501
    Perform main script logic.

    Parameters
    ----------
    pdb_files : str or Path, required
        Location for PDB files to operate on, can be within a folder or inside .TAR file

    source : string or Path, optional
        If given, updates a preexisting torsions.JSON file.
        Defaults to `None`.

    output : string or Path, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.

    degrees : Bool, optional
        Whether to use degrees instead of radians.
        Defaults to False, use radians.

    ncores : int
        The numbers of cores to use.

    plot : Bool, optional
        Whether to plot the fractional secondary structure information.
        Defaults to False, don't plot.

    plotvars : dictionary, optional
        Parameters for creating the torsion angle distribution plot

    ramaplot : dictionary, optional
        Parameters for creating the ramachandran fraction secondary structure plot
    """
    # validates before performing time consuming calculations
    if source and not source.suffix == '.json':
        raise ValueError('Source file should have `.json` extension.')

    output = output or 'torsions.json'
    if not output.endswith('.json'):
        raise ValueError('Output file should have `.json` extension.')

    log.info(T('Extracting torsion angles'))
    init_files(log, LOGFILESNAME)

    if source:
        database_dict = read_dictionary_from_disk(source)

    log.info(T('reading input paths'))
    pdbs = FileReaderIterator(pdb_files, ext='.pdb')
    log.info(S('done'))

    consume = partial(starunpack,
                      cli_helper_calc_torsions,
                      degrees=degrees,
                      decimals=10
                      )

    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_prefix=_name,
        )

    execute_pool = pool_function(execute, pdbs, ncores=ncores)

    torsion_result = {
        Path(pdbid).stem: angles
        for pdbid, angles in execute_pool
        }

    if source:

        pop_difference_with_log(database_dict, torsion_result)

        for key, value in torsion_result.items():
            # where value is a dictionary {'chi':, 'phi':, 'omega':}
            database_dict[key].update(value)

        save_dict_to_json(database_dict, output=output)

    else:
        save_dict_to_json(torsion_result, output=output)

    if plot:
        # Plotting torsion angle distributions
        log.info(T("Plotting torsion angle distribution"))
        log.info(S("Reminder: PDBs must be conformers of the same protein-system."))  # noqa: E501

        plotvars = plotvars or dict()

        tor_defaults = {
            'angtype': 'phi',
            'filename': 'plot_torsions.png'
            }
        tor_defaults.update(plotvars)
        ttype = tor_defaults['angtype']

        max_residues = 0
        n_confs = len(torsion_result)

        for key in torsion_result:
            if len(torsion_result[key]['omega']) + 1 > max_residues:
                max_residues = len(torsion_result[key]['omega']) + 1

        angles = np.ndarray(shape=(n_confs, max_residues), dtype=float)
        j = 0
        for t in torsion_result:
            for i in range(1, max_residues):
                try:
                    angles[j:, i - 1] = torsion_result[t][ttype][i - 1]
                except Exception:
                    angles[j:, i - 1] = 0
            j += 1

        errs = plot_torsions(max_residues, angles, degrees, n_confs, **tor_defaults)  # noqa: E501
        for e in errs:
            log.info(S(f'{e}'))
        log.info(S(f'saved plot: {tor_defaults["filename"]}'))

        # Plotting ramachandran frac sec. str.
        log.info(T("Plotting ramachandran fractional secondary structure"))

        newfname = "_ramaSS.".join(tor_defaults['filename'].rsplit('.', 1))
        rama_defaults = {
            'sstype': 'Rama.',
            'filename': newfname,
            }

        if ramaplot:
            ramaplot = values_to_dict(ramaplot)
            rama_defaults.update(ramaplot)

        frac_rama = {
            'alpha': np.zeros(max_residues),
            'beta': np.zeros(max_residues),
            'other': np.zeros(max_residues),
            }
        p = np.ones(n_confs)
        p = p / len(p)
        c = 0
        
        deglimits = {
            '-180': -180.0,
            '-120': -120.0,
            '10': 10.0,
            '45': 45.0,
            '180': 180.0
            }
        
        if not degrees:
            for deg in deglimits:
                deglimits[deg] = radians(deglimits[deg])
            
        for conf in torsion_result:
            for res in range(max_residues - 1):
                try:
                    phi = torsion_result[conf]["phi"][res]
                    psi = torsion_result[conf]["psi"][res]
                except Exception:
                    continue
                if (deglimits['-180'] < phi and phi < deglimits['10']) and (deglimits['-120'] < psi and psi < deglimits['45']):  # noqa: E501
                    frac_rama['alpha'][res] += p[c]
                elif (deglimits['-180'] < phi and phi < 0.0) and ((deglimits['-180'] < psi and psi < deglimits['-120']) or (deglimits['45'] < psi and psi < deglimits['180'])):  # noqa: E501
                    frac_rama['beta'][res] += p[c]
                else:
                    frac_rama['other'][res] += p[c]
            c += 1
        errs = plot_fracSS(max_residues, frac_rama, **rama_defaults)
        for e in errs:
            log.info(S(f'{e}'))
        log.info(S(f'saved plot: {rama_defaults["filename"]}'))

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
