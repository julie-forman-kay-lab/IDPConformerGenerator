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

#TODO: make an additional flag signaling if the user wants to group DSSP codes
# e.g. -g, --group to make L+, H+, E+ from the 8 different DSSP codes

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
    data = ''
    aa = ''
    resid = 0
    prob = 0.0
    idx = 0
    
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
    

def main(
    chespi_p8,
    output,
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
    #TODO: make a regex formatting matcher to see if user inputted the correct file
    
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