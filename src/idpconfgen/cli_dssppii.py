"""
Adds polyproline type-2 helix assignment to DSSP. For DSSP v2, v3 and v4.

USAGE:
    $ idpconfgen dssppii [--pdb-file]
    $ idpconfgen dssppii [--pdb-file] [--horiz] [--output]
"""
# Translated to Python3 by Nemo (Zi Hao @menoliu) Liu on Feb 3, 2022
# Upgrades: no longer have to specify where DSSP is installed, updated CLI
# based on python standards.
#
# REFERENCE
#   Please cite:
#   Mansiaux Y, Joseph AP, Gelly J-C, de Brevern AG (2011) Assignment of
#   PolyProline II conformation and analysis od sequence-structure
#   relationship. PLoS One 6(3): e18401. doi:10.1371/journal.pone.0018401
#
#   Chebrek R, Leonard S, de Brevern AG, Gelly J-C (2014)
#   PolyprOnline: polyproline helix II and secondary structure assignment
#   database. Database; Nov 7;2014 [pmid:25380779]
#
#
# Copyright Jean-Christophe Gelly (Jan 20 2012)
#
# jean-christophe.gelly@univ-paris-diderot.fr
#
# This software is a computer program whose purpose is to
# assign polyproline type II helix from dsspcmbi program output.
#
# This software is governed by the CeCILL-B license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-B
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-B license and that you accept its terms.
import argparse
import os
import re

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libdssp import (
    dssp_version_str,
    get_dssp_version,
    is_dssp4,
    remove_native_ppii,
    run_dssp,
    )
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = '.idpconfgen_dssppii'
_name = 'dssppii'
_help = 'Extracts DSSP-PPII secondary information from a PDB file.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-pdb',
    '--pdb_file',
    help="Path to the PDB file to operate on.",
    required=True,
    )

ap.add_argument(
    '-cmd',
    '--dssp_cmd',
    help=(
        "Command or path to the DSSP executable. "
        "DSSP v2, v3 and v4 are supported. Defaults to `mkdssp`."
        ),
    default="mkdssp",
    )

ap.add_argument(
    '-horiz',
    '--horizontal',
    help="Give output DSSP in 1D sequence fashion.",
    action='store_true',
    )

ap.add_argument(
    '-o',
    '--output',
    help="The output file containing DSSP with PPII added.",
    type=Path,
    default=None,
    )

libcli.add_argument_ncores(ap)


def dssp_ppii_assignment(
        pdb_file,
        dssp_cmd="mkdssp",
        dssp_version=None,
        native_ppii=False,
        ):
    """
    Run DSSP and include PPII as secondary structure in its output.

    Parameters
    ----------
    pdb_file : string
        Path to the PDB file to operate on as indicated by the user

    dssp_cmd : string
        The DSSP executable, or a command that runs it.
        Defaults to `mkdssp`.

    dssp_version : tuple of int, optional
        The DSSP version, as returned by
        :func:`idpconfgen.libs.libdssp.get_dssp_version`. If not given,
        it is read by running ``<dssp_cmd> --version``.

    native_ppii : bool, optional
        Only for DSSP v4, which assigns PPII helices itself. If ``True``,
        return the output of DSSP v4 as it is, without DSSP-PPII. If
        ``False`` (default), reset the PPII helices of DSSP v4 to loops
        and assign PPII with DSSP-PPII, as for DSSP v2 and v3.

    Returns
    -------
    tab_new_dssp : list
        List of every line of DSSP output, now with PPII identification
    """
    if dssp_version is None:
        dssp_version = get_dssp_version(dssp_cmd)

    dssp_output = run_dssp(pdb_file, dssp_cmd, dssp_version)
    tab_output = dssp_output.decode('utf-8', errors='replace').splitlines()

    if is_dssp4(dssp_version):
        if native_ppii:
            return tab_output
        tab_output = remove_native_ppii(tab_output)

    return assign_ppii(tab_output)


def assign_ppii(tab_output):
    """
    Parse DSSP output to include PPII as secondary structure.

    Parameters
    ----------
    tab_output : list of str
        List of every line of DSSP output, in the classic DSSP format.

    Returns
    -------
    tab_new_dssp : list
        List of every line of DSSP output, now with PPII identification
    """
    aa = ""
    ss = ""
    assign = 0
    i = 0
    line = ""

    chain = ""
    hash_chain = {}

    epsilon = 29
    cannonical_phi = -75
    cannonical_psi = 145

    sup_phi = cannonical_phi + epsilon
    inf_phi = cannonical_phi - epsilon

    sup_psi = cannonical_psi + epsilon
    inf_psi = cannonical_psi - epsilon

    tab_new_dssp = []

    # Parsing DSSP to detect Polyproline II Helix
    for line in tab_output:
        if re.search("#  RESIDUE AA STRUCTURE", line):
            assign = 1
        elif assign == 1:
            ss = line[16]
            aa = line[13]

            if "!" in aa:
                hash_chain[chain].insert(i, 0)
                i += 1
                continue

            index = line[0:5]

            chain = line[11]

            phi = line[103:109]
            psi = line[109:115]

            index = float(re.sub(r"\s", "", index))
            phi = float(re.sub(r"\s", "", phi))
            psi = float(re.sub(r"\s", "", psi))

            if chain == ' ':
                chain = "_"

            if chain not in hash_chain:
                i = 0
                hash_chain[chain] = []

            if \
                    phi <= sup_phi \
                    and phi >= inf_phi \
                    and psi <= sup_psi \
                    and psi >= inf_psi \
                    and ss == ' ':

                hash_chain[chain].insert(i, 1)
                if i != 0:
                    if hash_chain[chain][i - 1] >= 1:
                        hash_chain[chain][i] = 2
                        hash_chain[chain][i - 1] = 2
                else:
                    hash_chain[chain][i] = 0
            else:
                hash_chain[chain].insert(i, 0)

            i += 1

    # Print Modified Output
    i = 0
    assign = 0
    prev_chain = ''
    for line in tab_output:

        if re.search("#  RESIDUE AA STRUCTURE", line):
            assign = 1
            tab_new_dssp.append(line)
            i = 0

        elif assign == 1:
            aa = line[13]
            chain = line[11]

            if "!" in aa:
                i += 1
                tab_new_dssp.append(line)
                continue

            if chain != prev_chain:
                prev_chain = chain
                i = 0

            if hash_chain[chain][i] >= 2:
                line_list = list(line)
                line_list[16] = "P"
                line = "".join(line_list)

            tab_new_dssp.append(line)
            i += 1

        else:
            tab_new_dssp.append(line)

    return tab_new_dssp


def parsing_dssp(ref_tab_out_dssp):
    """
    Parse the DSSP-PPII information to convert into simplified 1D output.

    Parameters
    ----------
    ref_tab_out_dssp : list
        List of all the lines of DSSP-PPII output

    Returns
    -------
    tab_out : list
        Simplified DSSP-PPII output that only includes primary sequence
        sligned to secondary structure codes.
    """
    aa = ""
    AA = []

    tab_out = []
    ss = ""
    SS = []
    assign = 0

    dsspto3 = {
        "H": "H",
        "G": "G",
        "I": "I",
        "E": "E",
        "B": "B",
        "C": "C",
        "S": "S",
        "T": "T",
        " ": "-",
        "P": "P",
        }

    tab_out.append(">SEQ\n")

    for line in ref_tab_out_dssp:

        if re.search("#  RESIDUE AA STRUCTURE", line):
            assign = 1

        elif assign == 1:
            ss = line[16]
            aa = line[13]

            if "!" not in aa:
                SS.append(dsspto3[ss])
                AA.append(aa)

    tab_out.append(''.join(AA) + "\n")
    tab_out.append(">DSSPPII\n")
    tab_out.append(''.join(SS) + "\n")

    return tab_out


def main(
        pdb_file,
        output,
        dssp_cmd="mkdssp",
        horizontal=False,
        ncores=1,
        **kwargs,
        ):
    """
    Perform main logic of the script.

    Parameters
    ----------
    pdb_file : string, required
        A string to the path of a PDB to operate on

    output : string, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.

    dssp_cmd : string, optional
        If given runs dssp as per user specification, else runs dssp as
        default.  Defaults to `mkdssp`.

    horizontal : bool, optional
        If given the output will be a 1D line of DSSP-PPII codes, else
        output is full DSSP-PPII. Defaults to True.
    """
    init_files(log, LOGFILESNAME)
    log.info(T('reading and processing DSSP information...'))

    dssp_version = get_dssp_version(dssp_cmd)
    log.info(S(f'DSSP version: {dssp_version_str(dssp_version)}'))

    ref_tab_out_dssp = dssp_ppii_assignment(pdb_file, dssp_cmd, dssp_version)

    if horizontal:
        ref_tab_out_dssp_horiz = parsing_dssp(ref_tab_out_dssp)
        output_ = os.linesep.join(ref_tab_out_dssp_horiz)

    else:
        output_ = os.linesep.join(ref_tab_out_dssp)

    if output:
        log.info(S('saving DSSP-PPII output onto disk...'))
        with open(output, mode="w") as fout:
            fout.write(output_)
        log.info(S('done'))
    else:
        print(output_)  # noqa: T201

    log.info(S('done'))


if __name__ == '__main__':

    libcli.maincli(ap, main)
