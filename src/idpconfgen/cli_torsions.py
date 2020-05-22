"""
Extracts torsion angles from PDBs.

USAGE:
    $ idpconfgen torsions [PDBS]
"""
import argparse
import itertools as it
from multiprocessing.pool import Pool
from functools import partial

import numpy as np

from idpconfgen import Path, log
from idpconfgen.libs import libcli, libio, libpdb, libstructure, libcalc
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = 'idpconfgen_torsion'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
# https://stackoverflow.com/questions/24180527

ap.add_argument(
    'pdbs',
    help='PDB file list.',
    nargs='+',
    )

ap.add_argument(
    '-d',
    '--degrees',
    help='Write angles in degrees instead of radians. Defaults to False.',
    action='store_true',
    )


def _load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    """
    Execute main client function.

    Reads command line arguments and executes logic.
    """
    cmd = _load_args()
    main(**vars(cmd))


def main(pdbs, degrees=True, **kwargs):
    """Perform main script logic."""
    log.info(T('Extracting torsion angles'))
    init_files(log, LOGFILESNAME)

    log.info(T('reading input paths'))
    # tee is used to keep memory footprint low
    # though it would be faster to create a list from path_bundle
    pdbs = libio.read_path_bundle(pdbs, ext='.pdb')
    log.info(S('done'))

    with Pool() as pool:
        imuo = pool.imap_unordered(
            partial(get_torsions, degrees=degrees),
            pdbs,
            )

        for i in imuo:
            pass#log.info(f'Done {i}')

    return


def get_torsions(pdbfile, degrees=False):
    structure = libstructure.Structure(pdbfile)
    structure.build()
    structure.add_filter_record_name(('ATOM', 'HETATM'))
    structure.add_filter_backbone(minimal=True)

    data = structure.filtered_atoms

    # validates structure data
    # rare are the PDBs that produce errors, still they exist.
    # errors can be from a panoply of sources, that is why I decided
    # not to attempt correcting them and instead ignore and report.
    validation_error = libcalc.validate_array_for_torsion(data)
    if validation_error:
        errmsg = (
            f'Error found for pdb: \'{pdbfile}\'\n'
            f'{validation_error}\n'
            'ignoring this structure...\n\n'
            )
        log.error(errmsg)
        return f'ERROR for: {pdbfile}'

    coords = data[:, libstructure.cols_coords].astype(float)
    torsions = libcalc.calc_torsion_angles(coords)

    if degrees:
        torsions = np.degrees(torsions)

    CA_C = torsions[::3]
    C_N = torsions[1::3]
    N_CA = torsions[2::3]

    result = np.zeros((CA_C.size + 1, 4))
    result[:, 0] = np.char.strip(data[::3, libstructure.col_resSeq])
    result[1:, 1] = N_CA
    result[:-1, 2] = CA_C
    result[:-1, 3] = C_N

    np.savetxt(
        Path(pdbfile).with_suffix('.torsion'),
        result,
        fmt=['% 4d'] + ['% +20.12f'] * 3,
        header=f'resid,{"N_CA":>16},{"CA_C":>20},{"C_N":>20}',
        )

    return pdbfile


if __name__ == '__main__':
    maincli()
