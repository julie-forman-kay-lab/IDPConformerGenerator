# CheSPI REFERENCE: https://github.com/protein-nmr/CheSPI
# Nielsen JT, Mulder FAA. CheSPI: chemical shift secondary structure population inference. # noqa: E501
# J Biomol NMR. 2021 Jul;75(6-7):273-291. doi: 10.1007/s10858-021-00374-w.
# Epub 2021 Jun 19. PMID: 34146207.
#########################################################
# d2D REFERENCE: https://github.com/carlocamilloni/d2D
# Determination of Secondary Structure Populations in Disordered States of Proteins Using Nuclear Magnetic Resonance Chemical Shifts # noqa: E501
# Carlo Camilloni, Alfonso De Simone, Wim F. Vranken, and Michele Vendruscolo
# Biochemistry 2012 51 (11), 2224-2231
# DOI: 10.1021/bi3001825

"""# noqa: D400, D205, E501
Ability to parse probabilistic secondary-structure predictions based on user preference.

The output will be printed to the terminal window. To save the output to
a file use the `>` command. Note that output MUST have the extension .JSON.

USAGE:
    $ idpconfgen csssconv [--chespi_p8] [--delta2D] [--output]
    $ idpconfgen csssconv [--chespi_p8] [--delta2D] > [OUTPUT]
"""
import argparse
import json
import re

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = '.idpconfgen_csssconv'
_name = 'csssconv'
_help = ('Standardizes secondary-structure prediction output '
         'for custom secondary-structure sampling (CSSS).'
         )

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
    help="Path to the δ2D output file to operate on."
    )

ap.add_argument(
    '-f',
    '--full',
    help=('Parses the secondary structure prediction files '
          'as is, without grouping DSSP codes.'
          ),
    action='store_true'
    )

libcli.add_argument_output(ap)


def chespi_probs8_convert_full(p8):
    """# noqa: D205, D400, E501
    Parse the probs8_[ID].txt output from CheSPI as user
    configurable input file for CSSS.
    
    Parameters
    ----------
    p8 : string
        Path to the probs8_[ID].txt to operate on as indicated by the user.
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates secondary structure (H/G/I/E/ /T/S/B)
        as defined in DSSP. Values are their respective probabilities.
    """
    dict_out = {}
    dict_p8 = {}
    dssp_keys = ["H", "G", "I", "E", " ", "T", "S", "B"]
    
    with open(p8) as reader:
        for line in reader:
            pline = line.split()
            pline.pop(0)
            data = [float(i) for i in pline]
            resid = int(data[0])
            data.pop(0)
            for i, prob in enumerate(data):
                dict_p8[dssp_keys[i]] = prob
            dict_out[resid] = dict_p8
            dict_p8 = {}
    return dict_out


def chespi_probs8_convert_grouped(p8):
    """# noqa: D205, D400, E501
    Parse the probs8_[ID].txt output from CheSPI as user configurable input file for CSSS.
    
    Groups together DSSP secondary structures as per idpconfgen definitions.
    If a residue has multiple SS probabilities, they are summative per L, H, E definition.
    
    Parameters
    ----------
    p8 : string
        Path to the probs8_[ID].txt to operate on as indicated by the user
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates grouped secondary structure (L, H, E)
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
            # adds up the probabilities for SS codes
            # based on idpconfgen definitions
            dict_p8["L"] = round((data[3] + data[5] + data[6] + data[7] + data[8]), 4)  # noqa: E501
            dict_p8["H"] = round(data[1] + data[2], 4)
            dict_p8["E"] = data[4]
            dict_out[resid] = dict_p8
            dict_p8 = {}
            
    return dict_out


def d2D_convert_full(d2d):
    """# noqa: D205, D400, E501
    Parse the .TXT output from δ2D as user configurable input file for CSSS.
    
    Parameters
    ----------
    d2d : string
        Path to the δ2D .TXT to operate on as indicated by the user
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates secondary structure ( /H/E/P)
        as defined in DSSP. Values for the second key-layer are their
        respective probabilities.
    """
    dict_out = {}
    dict_probs = {}
    dssp_keys=["H", "E", " ", "P"]
    
    with open(d2d) as reader:
        for line in reader:
            if "#" not in line and line != "\n":
                pline = line.split()
                pline.pop(1)
                pline.pop()
                data = [float(i) for i in pline]
                resid = int(data[0])
                data.pop(0)
                for i, prob in enumerate(data):
                    dict_probs[dssp_keys[i]] = prob
                dict_out[resid] = dict_probs
                dict_probs = {}
            # if there aren't any predicted probabilities
            # for this residue, make them all equal
            elif re.match(r"\#+\d", line):
                sline = line.split()
                pline = sline[0].split("#")
                resid = int(pline[1].strip())
                for key in dssp_keys:
                    dict_probs[key] = 0.25
                dict_out[resid] = dict_probs
    return dict_out


def d2D_convert_grouped(d2d):
    """# noqa: D205, D400, E501
    Parse the .TXT output from δ2D as user configurable input file for CSSS.
    
    Groups together DSSP secondary structures as per idpconfgen definitions.
    If a residue has multiple SS probabilities, they are summative per L, H, E definition.
    
    Parameters
    ----------
    d2d : string
        Path to the δ2D .TXT to operate on as indicated by the user
    
    Returns
    -------
    dict_out : dictionary
        Nested dictionary where the first key layer indicates residue number
        and the second key later indicates grouped secondary structure (L, H, E)
        as defined in idpconfgen. Values for the second key-layer are their
        respective probabilities.
    """
    dict_out = {}
    dict_probs = {}
    dssp_rd_keys = ["H", "E", "L"]
    
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
                dict_probs["L"] = round((data[3] + data[4]), 3)
                dict_out[resid] = dict_probs
                dict_probs = {}
            # if there aren't any predicted probabilities
            # for this residue, make them all equal
            elif re.match(r"\#+\d", line):
                sline = line.split()
                pline = sline[0].split("#")
                resid = int(pline[1].strip())
                for key in dssp_rd_keys:
                    dict_probs[key] = 0.333
                dict_out[resid] = dict_probs
    return dict_out


def main(
    chespi_p8,
    delta2D,
    output,
    full=False,
    **kwargs,
        ):
    """# noqa: D205, D400, E501
    Perform main logic of the script.
    
    Parameters
    ----------
    chespi_p8 : string
        A string to the path of probs8_[ID].TXT output from CheSPI.
        
    delta2D : string
        A string to the path of the .TXT output from δ2D.
    
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
            # regex matching pattern for the format of CheSPI output files
            if not re.fullmatch(r"\s{3}[A-Z]{1}\d|\s|.{60,}\n", reader.readline()):     # noqa: E501
                log.info(S('Incorrect CheSPI input file. Please use probs8_[ID].txt'))  # noqa: E501
                return
            
        if full: converted_chespi = chespi_probs8_convert_full(chespi_p8)   # noqa: E701, E501
        else: converted_chespi = chespi_probs8_convert_grouped(chespi_p8)   # noqa: E701, E501
        
        _output = json.dumps(converted_chespi, indent=4)
        
    if delta2D:
        log.info(T('reading and processing delta2D predictions...'))
        # TODO: regex matching to see the format of alleged delta2D output files
        if full: converted_delta2D = d2D_convert_full(delta2D)   # noqa: E701
        else: converted_delta2D = d2D_convert_grouped(delta2D)   # noqa: E701
        
        _output = json.dumps(converted_delta2D, indent=4)

    if output:
        log.info(S('saving converted secondary-predictions output onto disk...'))   # noqa: E501
        with open(output, mode="w") as fout:
            fout.write(_output)
    else:
        print(_output)   # noqa: T201
    
    log.info(S('done'))
    

if __name__ == '__main__':
    
    libcli.maincli(ap, main)
