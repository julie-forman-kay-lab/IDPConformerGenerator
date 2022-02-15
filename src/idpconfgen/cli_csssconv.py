# CheSPI REFERENCE: https://github.com/protein-nmr/CheSPI
# Nielsen JT, Mulder FAA. CheSPI: chemical shift secondary structure population inference.
# J Biomol NMR. 2021 Jul;75(6-7):273-291. doi: 10.1007/s10858-021-00374-w.
# Epub 2021 Jun 19. PMID: 34146207.
"""
Parses probs8_[ID].txt output of CheSPI to generate a user-editable input
for the idpconfgen build -csss module.

The output will be printed to the terminal window. To save the output to
a file use the `>` command.

USAGE:
    $ idpconfgen csssconv [--chespi_p8]
    $ idpconfgen csssconv [--chespi_p8] > [OUTPUT]
"""
import argparse
import re

from idpconfgen import Path, log
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
    '-v',
    '--verbose',
    help="Parses the CheSPI probs8 file as is, without grouping DSSP.",
    action='store_true'
)

libcli.add_argument_output(ap)

def chespi_probs8_convert_verbose(p8):
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
        " " : [],
        "T" : [],
        "S" : [],
        "B" : []
    }
    
    with open(p8) as reader:
        for line in reader:
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
                    elif idx == 4: dict_p8[" "].append([aa, resid, prob])
                    elif idx == 5: dict_p8["T"].append([aa, resid, prob])
                    elif idx == 6: dict_p8["S"].append([aa, resid, prob])
                    elif idx == 7: dict_p8["B"].append([aa, resid, prob])
                    
    return primary, dict_p8

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
    header : string
        Header for the converted file containing the primary sequence and residue numbers horizontally
        
    dict_p8 : dictionary
        Dictionary of a list of lists (E.g. [AA, RESID, PROB]) for each
        grouped secondary structure (L, H, E, G) as indicated in idpconfgen
    """
    primary = 'PRIMARY SEQ: '
    
    dict_p8 = {
        "L+" :[],
        "H+" : [],
        "E+" : [],
        "G+" : []
    }
    
    with open(p8) as reader:
        for line in reader:
            data = line.split()
            aa = data[0]
            resid = int(data[1])
            Lprob = round((float(data[4]) + float(data[6]) + float(data[7]) + float(data[8]) + float(data[9])), 4)
            Gprob = float(data[3])
            Hprob = float(data[2])
            Eprob = float(data[5])
            primary += aa
            if Lprob != 0 : dict_p8["L+"].append([aa, resid, Lprob])
            if Hprob != 0 : dict_p8["H+"].append([aa, resid, Hprob])
            if Eprob != 0 : dict_p8["E+"].append([aa, resid, Eprob])
            if Gprob != 0 : dict_p8["G+"].append([aa, resid, Gprob])
            
    return primary, dict_p8
        

def main(
    chespi_p8,
    output,
    verbose=False,
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
    
    if verbose:
        output_, converted_chespi = chespi_probs8_convert_verbose(chespi_p8)
        output_ += "\nTYPE: CheSPI_full"
    else:
        output_, converted_chespi = chespi_probs8_convert_grouped(chespi_p8)
        output_ += "\nTYPE: CheSPI_rd"

    
    for key in converted_chespi.keys():
        output_ += (f"\n{key}: ")
        for value in converted_chespi[key]:
            output_ += f"({value[0]}{value[1]}, {value[2]})"
    
    if output:
        log.info(S('saving converted CheSPI output onto disk...'))
        with open(output, mode="w") as fout:
            fout.write(output_)
    else:
        print(output_)
    
    log.info(S('done'))
    

if __name__ == '__main__':
    
    libcli.maincli(ap, main)