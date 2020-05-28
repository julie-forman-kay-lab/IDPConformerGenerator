"""
Parses DSSP file from RCSB.

DSSP file is as provided at: https://cdn.rcsb.org/etl/kabschSander/ss_dis.txt.gz

USAGE:
    $ idpconfgen parse_dssp DSSP COMPARE -o OUTPUT
"""
import gzip
import argparse


from idpconfgen.libs import libcli
from idpconfgen import Path, log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import concatenate_entries
from idpconfgen.libs.libpdb import PDBList
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = '.parse_dssp'

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

ap.add_argument(
    '-c',
    '--compare',
    help='A list of PDBID_CHAIN codes to compare with.',
    )

ap.add_argument(
    '-o',
    '--output',
    help='The output file.',
    default=None,
    )

def _load_args():
    cmd = ap.parse_args()
    return cmd


def main(
        rcsb_dssp,
        compare=None,
        output=None,
        func=None,
        ):

    #ss = Path(rcsb_dssp).read_bytes().split('>')
    with gzip.open(rcsb_dssp, 'rb') as fin:
        ss = fin.read().decode('utf-8').split('>')

    print(len(ss))

    if compare:
        pdbids_to_read = concatenate_entries([compare])
        pdblist = PDBList(pdbids_to_read)
        pdbids = set(str(id_) for id_ in pdblist)
        compare = True
        log.info(S('parsed compare'))

    log.info(T('reading RCSB DSSP file'))
    datass, fasta = {}, {}
    for i, item in enumerate(ss[1:]):
        si = item.split('\n')
        pdbid, chain, _ = si[0].split(':')
        n = f'{pdbid}_{chain.upper()}'

        if compare and n not in pdbids:
            continue

        if i % 2 == 0:
            fasta[n] = ''.join(si[1:])
        else:
            datass[n] = ''.join(si[1:])

    with open('allPDB.dssp', 'w') as fout:
        fout.write('\n'.join(f'{k}|{v}' for k, v in datass.items()))

    with open('allPDB.fasta', 'w') as fout:
        fout.write('\n'.join(f'{k}|{v}' for k, v in fasta.items()))

    return


def maincli():
    """Command-line interface entry point."""
    cmd = _load_args()
    main(**vars(cmd))


if __name__ == '__main__':
    maincli()
