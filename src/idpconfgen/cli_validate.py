"""
Conformer validator.

USAGE:
    $ idpconfgen validator FOLDER
"""
import argparse
from functools import partial

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import FileReaderIterator
from idpconfgen.libs.libstructure import Structure
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.logger import T, init_files, report_on_crash


# are they going to be here?
import numpy as np
from idpconfgen.libs.libstructure import col_name, col_resSeq
from idpconfgen.libs.libvalidate import vdW_clash_common_preparation, vdW_clash_calc, vdW_clash
from idpconfgen.core.definitions import vdW_radii_dict



LOGFILESNAME = '.validate'

_name = 'validate'
_help = 'Validate IDP conformers according to criteria.'

_prog, _des, _us = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_us,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-st',
    '--same-type',
    help=(
        'Whether all conformers are of the same type. '
        'If TRUE, consider atoms labels equal in all conformers. '
        'If FALSE, evaluates labels for every conformer separately. '
        'Defaults to TRUE - most common case.'
        ),
    action='store_false',
    )

ap.add_argument(
    '-ele',
    '--elements-to-consider',
    help='Which elements to consider in the validation process',
    default=None,
    type=tuple,
    )


ap.add_argument(
    '-vdw',
    '--vdW-radii',
    help='The van der Waals radii set.',
    choices=list(vdW_radii_dict.keys()),
    )

libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        elements_to_consider=('C'),
        func=None,
        ncores=1,
        residues_apart=3,
        same_type=True,
        vdW_radii='tsai1999',
        ):
    """Perform main logic."""
    log.info(T('Validating conformers'))
    init_files(log, LOGFILESNAME)

    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')

    execute_kwargs = {
        'elements_to_consider': elements_to_consider,
        'residues_apart': residues_apart,
        'vdW_radii': vdW_radii,
        }

    if same_type:
        execute = main_same_type(*next(pdbs2operate), **execute_kwargs)

    else:
        execute = partial(
            report_on_crash,
            validate_conformer,
            ROC_prefix=_name,
            **execute_kwargs,
            )

    execute_pool = pool_function(
        partial(starunpack, execute),
        pdbs2operate,
        ncores=ncores,
        )

    for _i in execute_pool:
        pass



def main_same_type(
        name,
        pdb_data,
        **preparation_kwargs,
        ):

    atc_mask, pure_radii_sum, distances_apart = preparation(
        name,
        pdb_data,
        **preparation_kwargs,
        )

    validate_conformer_from_preparation(
        name,
        pdb_data,
        atc_mask,
        pure_radii_sum,
        distances_apart,
        )
    # done

    execute = partial(
        report_on_crash,
        validate_conformer_from_preparation,
        atc_mask=atc_mask,
        pure_radii_sum=pure_radii_sum,
        distances_apart=distances_apart,
        ROC_prefix=_name,
        )

    return execute


def preparation(
        name,
        pdb_data,
        elements_to_consider,
        residues_apart,
        vdW_radii,
        ):
    """."""
    _, atom_elements, res_numbers = get_vdW_needs_from_structure(pdb_data)

    return vdW_clash_common_preparation(
        atom_elements,
        elements_to_consider=elements_to_consider,
        residue_numbers=res_numbers,
        residues_apart=residues_apart,
        vdW_radii=vdW_radii,
        )


def get_vdW_needs_from_structure(pdb_data):
    s = Structure(pdb_data)
    s.build()
    da = s.data_array
    atom_names = da[:, col_name]
    atom_elements = atom_names.astype('<U1')
    res_numbers = da[:, col_resSeq].astype(np.int)
    return s, atom_elements, res_numbers


def validate_conformer_from_preparation(
        name,
        pdb_data,
        **calc_kwargs,
        #atc_mask,
        #pure_radii_sum,
        #distances_apart,
        ):
    # here I don't use kwargs because of the need to use partial

    s = Structure(pdb_data)
    s.build()
    coords = s.coords

    return vdW_clash_calc(
        coords,
        **calc_kwargs,
        #atc_mask=atc_mask,
        #pure_radii_sum=pure_radii_sum,
        #distances_apart=distances_apart,
        )


def validate_conformer(
        name,
        pdb_data,
        elements_to_consider,
        residues_apart,
        vdW_radii,
        ):
    """."""
    s, atom_elements, res_numbers = get_vdW_needs_from_structure(pdb_data)

    coords = s.coords

    vdW_clash(
        coords,
        atom_elements,
        elements_to_consider,
        res_numbers,
        residues_apart,
        vdW_radii=vdW_radii,
        )


if __name__ == '__main__':
    libcli.maincli(ap, main)