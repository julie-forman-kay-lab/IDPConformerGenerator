"""
Calculate fragment statistics for the input sequence in the database.

USAGE:
    $ idpconfgen stats -db <DATABASE> -seq <INPUTSEQ>
    $ idpconfgen stats -db <DATABASE> -seq <INPUTSEQ> -op <OUTPUT-PREFIX> -of <OUTPUT-FOLDER>
"""
import argparse, os, re, json
from pathlib import Path

from idpconfgen import log
from idpconfgen.components.residue_tolerance import add_substitution_groups
from idpconfgen.libs import libcli
from idpconfgen.libs.libbuild import build_regex_substitutions, make_combined_regex
from idpconfgen.libs.libfilter import regex_forward_no_overlap
from idpconfgen.libs.libio import read_dict_from_json
from idpconfgen.libs.libparse import get_mers
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = '.stats'

_name = 'stats'
_help = 'Gets statistics from DB and input sequence.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_db(ap)
libcli.add_argument_seq(ap)
libcli.add_argument_dloopoff(ap)
libcli.add_argument_dhelix(ap)
libcli.add_argument_dstrand(ap)
libcli.add_argument_dany(ap)
libcli.add_argument_duser(ap)
add_substitution_groups(ap)
libcli.add_argument_output_folder(ap)


OUTPUT_PREFIX = "count_stats"
ap.add_argument(
    '-op',
    '--output-prefix',
    help=(
        'The prefix for the output CSV files. One CSV file is generated '
        'for each fragment size. These files will have this prefix. '
        'Defaults to `count_stats`.'
        ),
    default=OUTPUT_PREFIX,
    type=str,
    )

# the func=None receives the `func` attribute from the main CLI interface
# defined at cli.py
def main(
        database,
        input_seq,
        dhelix=False,
        dstrand=False,
        dloop_off=False,
        duser=False,
        dany=False,
        residue_substitutions=None,
        output_prefix=OUTPUT_PREFIX,
        output_folder=None,
        func=None,
        **kwargs
    ):
    """
    Performs main client logic.

    The default function is to return the number/counts of hits
    for different fragment sequences within the IDPConfGen database.

    A secondary function denoted with `--search` can search through
    raw PDB files' headers for certain keywords and return the PDB IDs
    that had hits to the keywords.

    Parameters
    ----------
    database : str or Path
        Path to the database.JSON file for IDPConfGen.
        Required for default func.

    input_seq : str or Path
        Primary sequence or path to the FASTA file.
        Required for default func.

    output_prefix : str
        Prefix for the output files when comparing counts of fragment ID hits.
        Defaults to `count_stats`. Optional for default func.

    output_folder : str or Path
        Path of the output folder to store stats.CSV.
        Defaults to working directory. Optional for default func.
    """
    init_files(log, LOGFILESNAME)

    dloop = not dloop_off
    any_def_loops = any((dloop, dhelix, dstrand))
    non_overlapping_parameters = (any_def_loops, dany, duser)  # noqa: E501
    _sum = sum(map(bool, non_overlapping_parameters))

    if _sum > 1:
        emsg = (
            ' * ERROR * (dloop, dstrand, dhelix), dany, and duser '
            'are mutually exclusive.'
            )
        log.error(emsg)
        raise ValueError(emsg)
    elif _sum < 1:
        emsg = ' * ERROR * Give at least one sampling option.'
        log.error(emsg)
        raise ValueError(emsg)

    del _sum
    del non_overlapping_parameters

    log.info(T(f"reading database"))
    log.info(S(database))
    db = read_dict_from_json(database)
    log.info(S(f"The database as {len(db)} protein chains"))

    # joins all primary sequences together
    fastas = '|'.join(v["fasta"] for v in db.values())
    log.info(S(f"The total number of residues in the database is {len(fastas)}."))

    if any((dloop, dhelix, dstrand)):
        dssp_regexes = []
        if dloop: dssp_regexes.append("L")
        if dhelix: dssp_regexes.append("H")
        if dstrand: dssp_regexes.append("E")

        dssps = '|'.join(v["dssp"] for v in db.values())
        dsspre = make_combined_regex(dssp_regexes)
        log.info(T("Indexing individual secondary structure codes"))
        log.info(S(f"Searching for regex: {dsspre}"))
        ss_slices = regex_forward_no_overlap(dssps, dsspre)
        fastas = "|".join(fastas[s] for s in ss_slices)

    elif duser:
        log.info(T("Searching for user defined regex"))
        log.info(S(', '.join(duser)))
        dsspre = make_combined_regex(duser)
        # this is very advanced, users should know what they are doing :-)
        dssps = '|'.join(v["dssp"] for v in db.values())
        ss_slices = regex_forward_no_overlap(dssps, dsspre)
        fastas = "|".join(fastas[s] for s in ss_slices)

    else:
        if not dany:
            raise ValueError(
                "We expected `dany` to be true. "
                "Give at least one sampling option."
                )
        log.info(T("searching considering only primary structure"))


    # Creates the statistics of the database. How many fragments of each size does the database has.
    fasta_lists = [_f for _f in fastas.split("|") if not "X" in _f]
    mers = {}
    for i in range(1, 8):
        iset = mers.setdefault(i, set())
        for v in fasta_lists:
            cset = get_mers(v, i)
            iset.update(cset)


    # Note that `X` residue exists in the database.
    log.info(T("In the dataset there are"))
    for k, v in mers.items():
        log.info(S(f"{len(v)} different fragments of {k} residues, of {len(mers[1])**k}"))  # noqa: E501

    # Gets the input sequence different fragments
    log.info(T("loading input sequence:"))
    log.info(S(input_seq))
    imers = {}
    for i in range(1, 8):
        imers[i] = get_mers(input_seq, i)

    log.info(T("In the input sequence there are"))
    for k, v in imers.items():
        log.info(S(f"{len(v)} different fragments of {k} residues, of {20**k} possible"))

    # Searches the input sequence fragments in the database.
    log.info(T("Creating fragments statistics for input sequence"))
    if residue_substitutions:
        sub_map = json.dumps(residue_substitutions)[1:-1]
        log.info(S(f"residue replacement map: {sub_map}"))
    data_to_plot = {}

    for mer_len, set_of_mers in imers.items():
        labels, data = data_to_plot.setdefault(mer_len, ([], []))
        for mer in set_of_mers:
            if residue_substitutions:
                mer = build_regex_substitutions(mer, residue_substitutions)
            labels.append(mer)
            data.append(len(re.findall("(?=" + mer + ")", fastas)))

    log.info(S("done"))


    # Saves to files.
    log.info(T("saving CSV files to disk"))

    if output_folder is not None:
        pof = Path(output_folder)
        pof.mkdir(parents=True, exist_ok=True)
    else:
        pof = Path.cwd()

    for mer_len, labels_values in data_to_plot.items():
        to_save = sorted(zip(*labels_values), key=lambda x: x[1], reverse=True)
        fpath = Path(pof, f"{output_prefix}_{mer_len}.csv")
        with open(fpath, "w") as fout:
            fout.write(f"# fragment counts in the database for {input_seq}{os.linesep}")
            _ = os.linesep.join(','.join(map(str, pair)) for pair in to_save)
            fout.write(_)
        log.info(S(f"saved: {fpath}"))

    return

if __name__ == '__main__':
    libcli.maincli(ap, main)
