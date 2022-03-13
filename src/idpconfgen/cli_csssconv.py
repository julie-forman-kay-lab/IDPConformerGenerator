# CheSPI REFERENCE: https://github.com/protein-nmr/CheSPI
# Nielsen JT, Mulder FAA. CheSPI: chemical shift secondary structure population inference.
# J Biomol NMR. 2021 Jul;75(6-7):273-291. doi: 10.1007/s10858-021-00374-w.
# Epub 2021 Jun 19. PMID: 34146207.
#########################################################
# d2D REFERENCE: https://github.com/carlocamilloni/d2D
# Determination of Secondary Structure Populations in Disordered States of Proteins Using Nuclear Magnetic Resonance Chemical Shifts
# Carlo Camilloni, Alfonso De Simone, Wim F. Vranken, and Michele Vendruscolo
# Biochemistry 2012 51 (11), 2224-2231
# DOI: 10.1021/bi3001825
"""
Ability to parse probabilistic secondary-structure predictions based on user preference.

The output will be printed to the terminal window. To save the output to
a file use the `>` command. Note that output MUST have the extension .JSON.

USAGE:
    $ idpconfgen csssconv [--chespi_p8] [--delta2D] [--output]
    $ idpconfgen csssconv [--chespi_p8] [--delta2D] > [OUTPUT]
"""
import argparse
import re
import json

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = 'idpconfgen_csssconv'
_name = 'csssconv'
_help = 'Standardizes secondary-structure prediction output for CSSS.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-p8',
    '--chespi_p8',
    help="Path to the CheSPI probs8 file to operate on.",
    )

ap.add_argument(
    '-d2D',
    '--delta2D',
    help="Path to the delta2D output file to operate on."
)

ap.add_argument(
    '-f',
    '--full',
    help="Parses the secondary structure prediction files as is, without grouping DSSP.",
    action='store_true'
    )

libcli.add_argument_output(ap)

def chespi_probs8_convert_full(p8):
    """
    Parse the probs8_[ID].txt output from CheSPI as user configurable input file for CSSS.
    
    Parameters
    ----------
    p8 : string
        Path to the probs8_[ID].txt to operate on as indicated by the user
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates secondary structure (H/G/I/E/ /T/S/B) 
        as defined in DSSP. Values are their respective probabilities.
    """
    dict_out = {}
    dict_p8 = {}
    
    with open(p8) as reader:
        for line in reader:
            pline = line.split()
            pline.pop(0)
            data = [float(i) for i in pline]
            resid = int(data[0])
            for i in range(1, 9):
                prob = data[i]
                if i == 1: dict_p8["H"] = prob
                elif i == 2: dict_p8["G"] = prob
                elif i == 3: dict_p8["I"] = prob
                elif i == 4: dict_p8["E"] = prob
                elif i == 5: dict_p8[" "] = prob
                elif i == 6: dict_p8["T"] = prob
                elif i == 7: dict_p8["S"] = prob
                elif i == 8: dict_p8["B"] = prob
            dict_out[resid] = dict_p8
            dict_p8 = {}
    return dict_out

def chespi_probs8_convert_grouped(p8):
    """
    Parse the probs8_[ID].txt output from CheSPI as user configurable input file for CSSS.
    
    Groups together DSSP secondary structures as per idpconfgen definitions.
    If a residue has multiple SS probabilities, they are summative per L, H, E, G definition.
    
    Parameters
    ----------
    p8 : string
        Path to the probs8_[ID].txt to operate on as indicated by the user
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates grouped secondary structure (L, H, E, G) 
        as defined in idpconfgen. Values are their respective probabilities.
    """
    dict_out = {}
    dict_p8 = {}
    
    with open(p8) as reader:
        for line in reader:
            pline = line.split()
            pline.pop(0)
            data = [float(i) for i in pline]
            resid = int(data[0])
            
            dict_p8["L"] = round((data[3] + data[5] + data[6] + data[7] + data[8]), 4)
            dict_p8["H"] = round(data[1] + data[2], 4)
            dict_p8["E"] = data[4]
            dict_out[resid] = dict_p8
            dict_p8 = {}
            
    return dict_out

def d2D_convert_full(d2d):
    """
    Parse the .TXT output from delta2D as user configurable input file for CSSS.
    
    Parameters
    ----------
    d2d : string
        Path to the delta2D .TXT to operate on as indicated by the user
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates secondary structure (H/E/L/P) 
        as defined in DSSP. Values are their respective probabilities.
    """
    dict_out = {}
    dict_probs = {}
    
    with open(d2d) as reader:
        for line in reader:
            if "#" not in line and line != "\n":
                pline = line.split()
                pline.pop(1)
                pline.pop()
                data = [float(i) for i in pline]
                resid = int(data[0])
                dict_probs["H"] = data[1]
                dict_probs["E"] = data[2]
                dict_probs[" "] = data[3]
                dict_probs["P"] = data[4]
                dict_out[resid] = dict_probs
                dict_probs = {}
            elif re.match(r"\#+\d", line):
                sline = line.split()
                pline = sline[0].split("#")
                resid = int(pline[1].strip())
                dict_probs["H"] = 0.25
                dict_probs["E"] = 0.25
                dict_probs[" "] = 0.25
                dict_probs["P"] = 0.25
                dict_out[resid] = dict_probs
    return dict_out

def d2D_convert_grouped(d2d):
    """
    Parse the .TXT output from delta2D as user configurable input file for CSSS.
    
    Groups together DSSP secondary structures as per idpconfgen definitions.
    If a residue has multiple SS probabilities, they are summative per L, H, E definition.
    
    Parameters
    ----------
    d2d : string
        Path to the delta2D .TXT to operate on as indicated by the user
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates grouped secondary structure (L, H, E) 
        as defined in idpconfgen. Values are their respective probabilities.
    """
    dict_out = {}
    dict_probs = {}
    
    with open(d2d) as reader:
        for line in reader:
            if "#" not in line and line != "\n":
                pline = line.split()
                pline.pop(1)
                pline.pop()
                data = [float(i) for i in pline]
                resid = int(data[0])
                dict_probs["H"] = data[1]
                dict_probs["E"] = data[2]
                dict_probs["L"] = round((data[3]+data[4]), 3)
                dict_out[resid] = dict_probs
                dict_probs = {}
            elif re.match(r"\#+\d", line):
                sline = line.split()
                pline = sline[0].split("#")
                resid = int(pline[1].strip())
                dict_probs["H"] = 0.333
                dict_probs["E"] = 0.333
                dict_probs["L"] = 0.333
                dict_out[resid] = dict_probs
    return dict_out

def main(
    chespi_p8,
    delta2D,
    output,
    full=False,
    **kwargs,
        ):
    """
    Perform main logic of the script.
    
    Parameters
    ----------
    chespi_p8 : string
        A string to the path of probs8_[ID].txt output from CheSPI.
        
    delta2D : string
        A string to the path of the .TXT output from delta2D.
    
    output : string, optional
        If given, prints output to that file (must be .JSON), else prints to console.
        Defaults to `None`.
    
    full : boolean, optional
        If given, defaults to true and all of the DSSP codes are used.
        Defaults to False.
    """
    init_files(log, LOGFILESNAME)
    
    if chespi_p8:
        log.info(T('reading and processing CheSPI predictions...'))
        with open(chespi_p8) as reader:
            if not re.fullmatch(r"\s{3}[A-Z]{1}\d|\s|.{60,}\n", reader.readline()):
                log.info(S('Incorrect CheSPI input file. Please use probs8_[ID].txt'))
                return
            
        if full: converted_chespi = chespi_probs8_convert_full(chespi_p8)
        else: converted_chespi = chespi_probs8_convert_grouped(chespi_p8)
        
        _output = json.dumps(converted_chespi, indent=4)
        
    if delta2D:
        log.info(T('reading and processing delta2D predictions...'))
        if full: converted_delta2D = d2D_convert_full(delta2D)
        else: converted_delta2D = d2D_convert_grouped(delta2D)
        
        _output = json.dumps(converted_delta2D, indent=4)
    
    
    if output:
        log.info(S('saving converted secondary-predictions output onto disk...'))
        with open(output, mode="w") as fout:
            fout.write(_output)
    else:
        print(_output)
    
    log.info(S('done'))
    

if __name__ == '__main__':
    
    libcli.maincli(ap, main)