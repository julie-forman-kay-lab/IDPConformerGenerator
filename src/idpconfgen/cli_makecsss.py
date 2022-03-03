"""
Generates a CSSS.JSON file based on user specifications for secondary structure sampling on a per residue basis.

Has the ability to edit previously generated CSSS files as well with `--file-csss`.
The output will be printed to the terminal window. To save the output to
a file use the `>` or `--output` command. Note that output MUST have the extension .JSON.

USAGE:
    $ idpconfgen makecsss [--file-csss] [--custom-pattern] [--output]
    $ idpconfgen makecsss [--file-csss] [--custom-pattern] > [OUTPUT]
"""
import argparse
from collections import defaultdict
import re
import json

from idpconfgen.libs.libio import read_dict_from_json
from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = 'idpconfgen_makecsss'
_name = 'makecsss'
_help = 'Creates or edits a custom CSSS file per user specifications.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-fc',
    '--file-csss',
    help="Specifies the path of an existing CSSS.JSON file to edit.",
    )

ap.add_argument(
    '-cp',
    '--custom-pattern',
    help="Specifies the probability of secondary structure sampling on a per-residue basis.\n"
        "Please cover all residues in your sequence of interest.\n"
        "Overlapping residues with the same DSSP_REGEX will be overwritten by subsequent entries.\n"
        "E.g. -cp <RESIDUE or RANGE> <DSSP_REGEX> <PROBABILITY>|...",
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
    dict_out : dict-like
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates secondary structure codes
        as defined in DSSP. Values are their respective probabilities.
    """
    dict_out = defaultdict(dict)
    lst_cp = cp.split("|")
    for groups in lst_cp:
        split_gp = groups.split()
        if "-" in split_gp[0]:
            res_range = split_gp[0].split("-")
            start_res = int(res_range[0])
            end_res = int(res_range[1])
            for i in range(start_res, end_res+1):
                dict_out[i][split_gp[1]] = float(split_gp[2])
        else:
            num_res = int(split_gp[0])
            dict_out[num_res][split_gp[1]] = float(split_gp[2])
        
    return dict_out

def main(
    file_csss,
    custom_pattern,
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
    
    
    converted_csss = parse_custom_pattern(custom_pattern)
    
    if file_csss:
        dictCSSS = read_dict_from_json(file_csss)
        for resid in converted_csss:
            for ss_regex in converted_csss[resid]:
                try:
                    dictCSSS[str(resid)][ss_regex] = converted_csss[resid][ss_regex]
                except KeyError:
                    log.info(S(f'Warning: new residue {resid} has been appended'))
                    dictCSSS[str(resid)] = {ss_regex : converted_csss[resid][ss_regex]}
        _output = json.dumps(dictCSSS, indent=4)
    else:
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