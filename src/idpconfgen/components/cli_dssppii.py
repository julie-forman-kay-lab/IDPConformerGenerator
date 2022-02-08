"""
Adds polyproline type-2 helix assignment to DSSP.

USAGE:
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
#   PolyprOnline: polyproline helix II and secondary structure assignment database.
#   Database ; Nov 7;2014 [pmid:25380779]
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
from idpconfgen.libs.libio import save_pairs_to_disk
from idpconfgen.logger import T, init_files


LOGFILESNAME = 'idpconfgen_dssppii'
_name = 'dssppii'
_help = 'Extracts DSSP-PPII secondary information from a PDB file.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

#TODO: make this iterable with a folder/.tar of PDBs
ap.add_argument(
    '-pdb',
    '--pdb_file',
    help="Path to the PDB file to operate on.",
    required=True,
    )

ap.add_argument(
    '-cmd',
    '--dssp_cmd',
    help="Path the do DSSP executable. Defaults to `dssp`.",
    default="dssp",
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


def dssp_ppii_assignment(pdb_file, dssp_cmd="dssp"):
    """
    Parsing DSSP output to include PPII as secondary structure.

    Parameters
    ----------
    pdb_file : string
        Path to the PDB file to operate on as indicated by the user
        
    dssp_cmd : string
        Any custom dssp commands the user may want. Defaults to `dssp`

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

    #Launch DSSP
    run_dssp = os.popen(dssp_cmd + " -i " + pdb_file).read()
    tab_output = list(run_dssp.split("\n"))
    tab_output.pop()

    #Parsing DSSP to detect Polyproline II Helix
    for line in tab_output:
        if re.search("#  RESIDUE AA STRUCTURE", line):
            assign = 1
        elif assign == 1:
            ss = line[16]
            aa = line[13]

            if re.search("\!", aa):
                hash_chain[chain].insert(i, 0)
                i += 1
                continue

            index = line[0:5]

            chain = line[11]

            phi = line[103:109]
            psi = line[109:115]

            index = int(float(re.sub("\s", "", index)))
            phi = int(float(re.sub("\s", "", phi)))
            psi = int(float(re.sub("\s", "", psi)))

            if chain==' ':
                chain = "_"

            if  chain not in hash_chain:
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

    #Print Modified Output
    i = 0
    assign = 0
    for line in tab_output:

        if re.search("#  RESIDUE AA STRUCTURE", line):
            assign = 1
            tab_new_dssp.append(line)
            i = 0

        elif assign == 1:
            aa = line[13]
            chain = line[11]

            if re.search("\!", aa):
                i += 1
                tab_new_dssp.append(line)
                continue

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
        Simplified DSSP-PPII output that only includes primary sequence sligned to secondary structure codes
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

        if line.startswith("  #  RESIDUE AA STRUCTURE"):
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
        dssp_cmd="dssp",
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
        If given runs dssp as per user specification, else runs dssp as default.
        Defaults to `dssp`.
    
    horizontal : bool, optional
        If given the output will be a 1D line of DSSP-PPII codes, else output is full DSSP-PPII
        Defaults to True.
    """
    log.info(T('Reading and processing DSSP information...'))
    init_files(log, LOGFILESNAME)
    
    ref_tab_out_dssp = dssp_ppii_assignment(pdb_file, dssp_cmd)

    if not horizontal:
        log.info(T('Printing DSSP-PPII full output'))
        for line in ref_tab_out_dssp:
            print(line)
    else:
        log.info(T('Converting DSSP-PPII output to 1D simplified version'))
        ref_tab_out_dssp_horiz = parsing_dssp(ref_tab_out_dssp)
        for line in ref_tab_out_dssp_horiz:
            print(line)

    if output is not None:
        log.info(T('Saving DSSP-PPII output onto disk'))
        save_pairs_to_disk(ref_tab_out_dssp, output)


if __name__ == '__main__':

    libcli.maincli(ap, main)
