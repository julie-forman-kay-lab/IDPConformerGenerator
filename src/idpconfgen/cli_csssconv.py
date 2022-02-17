# CheSPI REFERENCE: https://github.com/protein-nmr/CheSPI
# Nielsen JT, Mulder FAA. CheSPI: chemical shift secondary structure population inference.
# J Biomol NMR. 2021 Jul;75(6-7):273-291. doi: 10.1007/s10858-021-00374-w.
# Epub 2021 Jun 19. PMID: 34146207.
"""
Parses probs8_[ID].txt output of CheSPI to generate a user-editable input
for the idpconfgen build -csss module.

The output will be printed to the terminal window. To save the output to
a file use the `>` command. Note that output MUST have the extension .JSON.

USAGE:
    $ idpconfgen csssconv [--chespi_p8]
    $ idpconfgen csssconv [--chespi_p8] > [OUTPUT]
"""
import argparse
import re
import json

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = 'idpconfgen_csssconv'
_name = 'csssconv'
_help = 'Standardizes CheSPI probs8 output for CSSS.'

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
    required=True,
    )

ap.add_argument(
    '-f',
    '--full',
    help="Parses the CheSPI probs8 file as is, without grouping DSSP.",
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
        and the second key later indicates secondary structure (H/G/I/E/-/T/S/B) 
        as defined in DSSP. Values are their respective probabilities.
    """
    dict_out = {}
    dict_p8 = {
        "H" : 0.0,
        "G" : 0.0,
        "I" : 0.0,
        "E" : 0.0,
        " " : 0.0,
        "T" : 0.0,
        "S" : 0.0,
        "B" : 0.0
    }
    
    with open(p8) as reader:
        for line in reader:
            data = line.split()
            resid = int(data[1])
            for i in range(2, 10):
                prob = float(data[i])
                if prob != 0:
                    idx = i - 2
                    if idx == 0: dict_p8["H"] = prob
                    elif idx == 1: dict_p8["G"] = prob
                    elif idx == 2: dict_p8["I"] = prob
                    elif idx == 3: dict_p8["E"] = prob
                    elif idx == 4: dict_p8[" "] = prob
                    elif idx == 5: dict_p8["T"] = prob
                    elif idx == 6: dict_p8["S"] = prob
                    elif idx == 7: dict_p8["B"] = prob
                else:
                    idx = i - 2
                    if idx == 0: dict_p8["H"] = 0.0
                    elif idx == 1: dict_p8["G"] = 0.0
                    elif idx == 2: dict_p8["I"] = 0.0
                    elif idx == 3: dict_p8["E"] = 0.0
                    elif idx == 4: dict_p8[" "] = 0.0
                    elif idx == 5: dict_p8["T"] = 0.0
                    elif idx == 6: dict_p8["S"] = 0.0
                    elif idx == 7: dict_p8["B"] = 0.0
            dict_out[resid] = dict_p8
            dict_p8 = {"H" : 0.0,"G" : 0.0,"I" : 0.0,"E" : 0.0," " : 0.0,"T" : 0.0,"S" : 0.0,"B" : 0.0}
    return dict_out

def chespi_probs8_convert_grouped(p8):
    """
    Parse the probs8_[ID].txt output from CheSPI as user configurable input file for CSSS.
    
    Groups together DSSP secondary structures as per idpconfgen definitions.
    If a residue has multiple SS probabilities, they are summative per L+, H+, E+ definition.
    
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
    dict_p8 = {
        "L+" : 0.0,
        "H+" : 0.0,
        "E+" : 0.0,
        "G+" : 0.0
    }
    
    with open(p8) as reader:
        for line in reader:
            data = line.split()
            resid = int(data[1])
            Lprob = round((float(data[4]) + float(data[6]) + float(data[7]) + float(data[8]) + float(data[9])), 4)
            Gprob = float(data[3])
            Hprob = float(data[2])
            Eprob = float(data[5])
            
            if Lprob != 0 : dict_p8["L+"] = Lprob 
            else: dict_p8["L+"] = 0.0
            if Hprob != 0 : dict_p8["H+"] = Hprob
            else: dict_p8["H+"] = 0.0
            if Eprob != 0 : dict_p8["E+"] = Eprob
            else: dict_p8["E+"] = 0.0
            if Gprob != 0 : dict_p8["G+"] = Gprob
            else: dict_p8["G+"] = 0.0
            dict_out[resid] = dict_p8
            dict_p8 = {"L+" : 0.0,"H+" : 0.0,"E+" : 0.0,"G+" : 0.0}
            
    return dict_out
        

def main(
    chespi_p8,
    output,
    full=False,
    **kwargs,
        ):
    """
    Perform main logic of the script.
    
    Parameters
    ----------
    chespi_p8 : string, required
        A string to the path of probs8_[ID].txt output from CheSPI.
    
    output : string, optional
        If given, prints output to that file (must be .JSON), else prints to console.
        Defaults to `None`.
    
    full : boolean, optional
        If given, defaults to true and all of the CheSPI DSSP codes are used.
        Defaults to False.
    """
    init_files(log, LOGFILESNAME)
    log.info(T('reading and processing CheSPI probs8 output...'))

    with open(chespi_p8) as reader:
        line = reader.readline()
        if not re.fullmatch(r"\s{3}[A-Z]{1}\d|\s|.{60,}\n", line):
            log.info(S('Incorrect CheSPI input file. Please use probs8_[ID].txt'))
            return
    
    if full:
        converted_chespi = chespi_probs8_convert_full(chespi_p8)
    else:
        converted_chespi = chespi_probs8_convert_grouped(chespi_p8)
    
    if output:
        log.info(S('saving converted CheSPI output onto disk...'))
        with open(output, mode="w") as fout:
            json.dump(converted_chespi, fout, indent=4)
    else:
        print(converted_chespi)
    
    log.info(S('done'))
    

if __name__ == '__main__':
    
    libcli.maincli(ap, main)