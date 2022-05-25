"""
Calculate number of hits of keywords in the header from raw PDB files fetched.

USAGE:
    $ idpconfgen search <PDB-FILES> -kw <KEYWORDS>
    $ idpconfgen search <PDB-FILES> -kw <KEYWORDS> -o <OUTPUT> -n <CORES>
"""
import argparse, os, json, shutil
from functools import partial
from pathlib import Path

from idpconfgen import log
from idpconfgen.libs.libio import (
    extract_from_tar,
    read_path_bundle,
    )
from idpconfgen.libs import libcli
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.logger import S, T, init_files, report_on_crash


LOGFILESNAME = '.idpconfgen_search'
_name = 'search'
_help = 'Searches through PDB headers for keywords of interest.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-kw',
    '--keywords',
    help=(
        "Keywords to look for in the header of PDB files. "
        "Delimited by commas, e.g. -kw 'pro,beta,elastin'"
        ),
    nargs='?',
    )

ap.add_argument(
    '-o',
    '--output',
    help=(
        "A path to the file where the keyword hits in the PDB headers "
        "dictionary will be saved. "
        "Defaults to dbHits.json, requires \'.json\' extension."
        ),
    type=Path,
    default='search_hits.json',
    action=libcli.CheckExt({'.json'}),
    )

TMPDIR = '__tmpsscalc__'
ap.add_argument(
    '-tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )

libcli.add_argument_ncores(ap)


def find_in_header(kwds, file_):
    with open(file_) as r:
        line = r.readline().lower()
        if "header" not in line:
            return None
        while "expdta" not in line:
            for word in kwds:
                if word in line:
                    return word, file_

            line=r.readline().lower()

    return None


def main(
        pdb_files,
        keywords,
        output="search_hits.json",
        tmpdir=TMPDIR,
        ncores=1,
        **kwargs,
        ):
    """
    Performs main client logic.

    Searches through raw PDB files' headers for certain
    keywords and return the PDB IDs that had hits to the keywords.

    Parameters
    ----------
    pdb_files : str or Path
        Path to a .TAR or folder of PDB files.
        Required for `--search`.

    keywords : str
        String of keywords delimited by commas to search for.
        Required for `--search`.

    output (str, optional): str or Path
        Path to the file to store the number of search hits.
        Defaults to "search_hits.json". Optional for `--search`.

    tmpdir : str or Path
        Path to the temporary directory if working with .TAR files.
        Defaults to TMPDIR. Optional for `--search`.

    ncores : int
        The number of cores to use.
        Defaults to 1. Optional for `--search`.
    """
    init_files(log, LOGFILESNAME)

    if not pdb_files or not keywords:
        log.info(S(
            'To use the search function you must provide PDB files '
            'and a list of keywords.'))
        return

    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, FileNotFoundError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))

    log.info(T('preparing task execution'))
    try:
        keywords = keywords.split(",")
    except:
        log.info(S('ERROR: keywords must be delimited by commas.'))
        return

    execute = partial(
        report_on_crash,
        find_in_header,
        keywords,
        )
    execute_pool=pool_function(execute, pdbs2operate, ncores=ncores)

    matches={}
    for word in keywords:
        matches[word]=[]

    for result in execute_pool:
        if result != None:
            wd = result[0]
            id = os.path.basename(result[1])
            matches[wd].append(id)

    for wd in keywords:
        log.info(S(f'total of {len(matches[wd])} instances of {wd}.'))

    log.info(S('done'))

    log.info(T('writing matches to disk'))
    _output = json.dumps(matches, indent = 4)
    with open(output, mode="w") as fout:
        fout.write(_output)
    log.info(S('done'))

    if _istarfile:
        shutil.rmtree(tmpdir)

    return


if __name__ == '__main__':
    libcli.maincli(ap, main)
