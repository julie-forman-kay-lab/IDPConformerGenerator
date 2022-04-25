"""
Calculate fragment statistics for the input sequence in the database.

USAGE:
    $ idpconfgen stats -db <DATABASE> -seq <INPUTSEQ>
"""
import argparse
import os
import re

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import read_dict_from_json
from idpconfgen.libs.libparse import get_mers
from idpconfgen.logger import S, T, init_files, pre_msg, report_on_crash


LOGFILESNAME = '.db_stats'

_name = 'stats'
_help = 'Gets statistics from DB and input sequence.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_idb(ap)
libcli.add_argument_seq(ap)


# the func=None receives the `func` attribute from the main CLI interface
# defined at cli.py
def main(database, input_seq, func=None):
    """Perform main logic."""
    log.info(T(f"reading database"))
    log.info(S(database))
    db = read_dict_from_json(database)
    log.info(S(f"The database as {len(db)} protein chains"))

    # joins all primary sequences together
    fastas = '|'.join(v["fasta"] for v in db.values())
    log.info(S(f"The total number of residues in the database is {len(fastas)}."))


    # Creates the statistics of the database. How many fragments of each size does the database has.
    mers = {}
    for i in range(1, 8):
        iset = mers.setdefault(i, set())
        for v in db.values():
            cset = get_mers(v["fasta"], i)
            # cull database usually have "X" residues
            # remove them.
            for xvalue in list(cset):
                if "X" in xvalue:
                    cset.remove(xvalue)
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
    data_to_plot = {}

    for mer_len, set_of_mers in imers.items():
        labels, data = data_to_plot.setdefault(mer_len, ([], []))
        for mer in set_of_mers:
            labels.append(mer)
            data.append(len(re.findall("(?=" + mer + ")", fastas)))

    log.info(S("done"))


    # Saves to files.
    log.info(T("saving CSV files to disk"))
    for mer_len, labels_values in data_to_plot.items():
        to_save = sorted(zip(*labels_values), key=lambda x: x[1], reverse=True)
        fpath = f"counts_in_database_{mer_len}.csv"
        with open(fpath, "w") as fout:
            fout.write(f"# fragment counts in the database for {input_seq}{os.linesep}")
            _ = os.linesep.join(','.join(map(str, pair)) for pair in to_save)
            fout.write(_)
        log.info(S(f"saved: {fpath}"))

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
