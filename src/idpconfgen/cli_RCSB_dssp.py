"""
Parses DSSP file from RCSB.

DSSP file is as provided at:
    https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz

USAGE:
    $ idpconfgen parse_dssp DSSP COMPARE -o OUTPUT
"""
import argparse
import gzip
from collections import defaultdict

from idpconfgen import Path, log
from idpconfgen.core.definitions import dssp_trans
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import save_dict_to_json
from idpconfgen.logger import T


_name = 'parse_dssp'
_help = 'Parse DSSP from RCSB.'
_prog, _des, _us = libcli.parse_doc_params(__doc__)


ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    'rcsb_dssp',
    help='The RCSB dssp file Gzipped.',
    )

libcli.add_argument_reduced(ap)

ap.add_argument(
    '-o',
    '--output',
    help=(
        "The output file containing the PDBID and "
        "respective DSSP and FASTA sequence information. "
        "Defaults to `allPDB.json`"
        ),
    type=Path,
    default='allPDB.json',
    action=libcli.CheckExt({'.json'}),
    )


def main(rcsb_dssp, func=None, output=None, reduced=False):
    """Perform main logic."""
    try:
        with gzip.open(rcsb_dssp, 'rb') as fin:
            ss = fin.read().decode('utf-8').split('>')
    except gzip.BadGzipFile:
        with open(rcsb_dssp, 'r') as fin:
            ss = fin.read().split('>')

    log.info(T('reading RCSB DSSP file'))
    allpdb = defaultdict(dict)
    for i, item in enumerate(ss[1:]):
        si = item.split('\n')
        pdbid, chain, _ = si[0].split(':')
        n = f'{pdbid}_{chain}'

        if i % 2 == 0:
            allpdb[n]['fasta'] = ''.join(si[1:])
        else:
            allpdb[n]['dssp'] = ''.join(si[1:])

    if reduced:
        _DT = dssp_trans
        for key in allpdb:
            allpdb[key]['dssp'] = allpdb[key]['dssp'].translate(_DT)

    save_dict_to_json(allpdb, output)
    return


if __name__ == '__main__':
    libcli.maincli()
