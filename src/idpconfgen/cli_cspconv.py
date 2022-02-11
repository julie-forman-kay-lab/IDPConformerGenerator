"""
Parses probs8_[ID].txt output of CheSPI to generate a user-editable input
for the idpconfgen build -csp module.

The output will be printed to the terminal window. To save the output to
a file use the `>` command.

USAGE:
    $ idpconfgen cspconv [--chespi_p8]
    $ idpconfgen cspconv [--chespi_p8] > [OUTPUT]
"""
import argparse
import re

from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = 'idpconfgen_cspconv'
_name = 'cspconv'
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
    '-g',
    '--group',
    help="Groups DSSP secondary structures into L+, H+, or E+.",
    action="store_true",
)

ap.add_argument(
    '-o',
    '--output',
    help="The user editable output file.",
    type=Path,
    default=None,
    )

#Formatting may change depending on how CSSS is implemented
def chespi_probs8_convert(p8):
    """
    Parse the probs8_[ID].txt output from CheSPI as user configurable input file for CSSS.
    
    Parameters
    ----------
    p8 : string
        Path to the probs8_[ID].txt to operate on as indicated by the user
    
    Returns
    -------
    header : string
        Header for the converted file containing the primary sequence and residue numbers horizontally
        
    dict_p8 : dictionary
        Dictionary of a list of lists (E.g. [AA, RESID, PROB]) for each
        possible secondary structure dssp code (H/G/I/E/-/T/S/B) as indicated
        in probs8_[ID].txt
    """
    primary = 'PRIMARY SEQ: '
    dict_p8 = {
        "H" : [],
        "G" : [],
        "I" : [],
        "E" : [],
        "-" : [],
        "T" : [],
        "S" : [],
        "B" : []
    }
    
    with open(p8) as reader:
        Lines = reader.readlines()
        for line in Lines:
            data = line.split()
            aa = data[0]
            resid = int(data[1])
            primary += aa
            for i in range(2, 10):
                prob = float(data[i])
                if prob != 0:
                    idx = i - 2
                    if idx == 0: dict_p8["H"].append([aa, resid, prob])
                    elif idx == 1: dict_p8["G"].append([aa, resid, prob])
                    elif idx == 2: dict_p8["I"].append([aa, resid, prob])
                    elif idx == 3: dict_p8["E"].append([aa, resid, prob])
                    elif idx == 4: dict_p8["-"].append([aa, resid, prob])
                    elif idx == 5: dict_p8["T"].append([aa, resid, prob])
                    elif idx == 6: dict_p8["S"].append([aa, resid, prob])
                    elif idx == 7: dict_p8["B"].append([aa, resid, prob])
                    
    return primary, dict_p8

def group_ss_structures(p8):
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
    header : string
        Header for the converted file containing the primary sequence and residue numbers horizontally
        
    dict_p8 : dictionary
        Dictionary of a list of lists (E.g. [AA, RESID, PROB]) for each
        grouped secondary structure (L+, H+, E+) as indicated in idpconfgen
    """
    primary = 'PRIMARY SEQ: '
    
    dict_p8 = {
        "L+" :[],
        "H+" : [],
        "E+" : []
    }
    
    with open(p8) as reader:
        Lines = reader.readlines()
        for line in Lines:
            data = line.split()
            aa = data[0]
            resid = int(data[1])
            Lprob = float(data[3]) + float(data[4]) + float(data[6]) + float(data[7]) + float(data[8]) + float(data[9])
            Hprob = float(data[2])
            Eprob = float(data[5])
            primary += aa
            if Lprob != 0 : dict_p8["L+"].append([aa, resid, Lprob])
            if Hprob != 0 : dict_p8["H+"].append([aa, resid, Hprob])
            if Eprob != 0 : dict_p8["E+"].append([aa, resid, Eprob])
            
    return primary, dict_p8
        

def main(
    chespi_p8,
    output,
    group=False,
    **kwargs,
        ):
    """
    Perform main logic of the script.
    
    Parameters
    ----------
    chespi_p8 : string, required
        A string to the path of probs8_[ID].txt output from CheSPI.
    
    output : string, optional
        If given prints output to that file, else prints to console.
        Defaults to `None`.
    """
    init_files(log, LOGFILESNAME)
    log.info(T('reading and processing CheSPI probs8 output...'))

    with open(chespi_p8) as reader:
        line = reader.readline()
        if not re.fullmatch(r"\s{3}[A-Z]{1}\d|\s|.{60,}\n", line):
            log.info(S('Incorrect CheSPI input file. Please use probs8_[ID].txt'))
            return
    
    if group:
        output_, converted_chespi = group_ss_structures(chespi_p8)
    else:
        output_, converted_chespi = chespi_probs8_convert(chespi_p8)

    for key in converted_chespi.keys():
        output_ += ("\n%s: " % key)
        for value in converted_chespi[key]:
            output_ += "(%s%d, %f)" % (value[0], value[1], value[2])
    
    if output:
        log.info(S('saving converted CheSPI output onto disk...'))
        with open(output, mode="w") as fout:
            fout.write(output_)
        log.info(S('done'))
    else:
        print(output_)
    
    log.info(S('done'))
    

if __name__ == '__main__':
    
    libcli.maincli(ap, main)