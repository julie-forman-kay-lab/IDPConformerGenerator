"""
Generates a CSSS.JSON file based on user specifications for secondary structure sampling
on a per residue basis.

The output will be printed to the terminal window. To save the output to
a file use the `>` command. Note that output MUST have the extension .JSON.

USAGE:
    $ idpconfgen makecsss [--custom-pattern]
    $ idpconfgen makecsss [--custom-pattern] > [OUTPUT]
"""
import argparse
import re
import json

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = 'idpconfgen_makecsss'
_name = 'makecsss'
_help = 'Creates a custom CSSS file per user specifications.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-cp',
    '--custom-pattern',
    help="Specifies the probability of secondary structure sampling on a per-residue basis.\n"+
        "E.g. -cp <DSSP_REGEX> <PROBABILITY> <RESIDUE-RANGE>|...",
    required=True,
    )


libcli.add_argument_output(ap)

def parse_custom_pattern(cp):
    """
    Parse the probs8_[ID].txt output from CheSPI as user configurable input file for CSSS.
    
    Parameters
    ----------
    cp : string
        User specified probabilities of secondary structures on a per residue basis.
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates secondary structure codes
        as defined in DSSP. Values are their respective probabilities.
    """
    dict_out = {}
    
    return dict_out

def main(
    custom_parameters,
    output,
    **kwargs,
        ):
    """
    Perform main logic of the script.
    
    Parameters
    ----------
    custom_parameters : string, required
        User specified probabilities of secondary structures on a per residue basis.
    
    output : string, optional
        If given, prints output to that file (must be .JSON), else prints to console.
        Defaults to `None`.
    """
    init_files(log, LOGFILESNAME)
    log.info(T('processing custom CSSS pattern...'))

    converted_csss = parse_custom_pattern(custom_parameters)
    
    _output = json.dumps(converted_csss, indent=4)
    
    if output:
        log.info(S('saving custom CSSS output onto disk...'))
        with open(output, mode="w") as fout:
            fout.write(_output)
    else:
        print(_output)
    
    log.info(S('done'))
    

if __name__ == '__main__':
    
    libcli.maincli(ap, main)