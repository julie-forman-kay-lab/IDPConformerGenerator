"""
Conformer validator.

USAGE:
    $ idpconfgen validator FOLDER
"""
import argparse
from functools import partial

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libdownload import fetch_raw_PDBs
from idpconfgen.libs.libio import FileReaderIterator
from idpconfgen.libs.libhigherlevel import download_pipeline
from idpconfgen.libs.libstructure import Structure
from idpconfgen.libs.libvalidate import validate_conformer
from idpconfgen.libs.libmulticore import pool_function, starunpack
from idpconfgen.logger import T, init_files, report_on_crash


# are they going to be here?
import numpy as np
from idpconfgen.libs.libstructure import col_name, col_resSeq
from idpconfgen.libs.libvalidate import vdW_clash_common_preparation, vdW_clash_calc
from idpconfgen.core.definitions import vdW_radii_tsai_1999



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

libcli.add_argument_ncores(ap)


def main(
        pdb_files,
        func=None,
        ncores=1,
        residues_apart=3,):
    """Perform main logic."""
    log.info(T('Validating conformers'))
    init_files(log, LOGFILESNAME)

    pdbs2operate = FileReaderIterator(pdb_files, ext='.pdb')
    # structures = (pre_validation(name, pdb) for name, pdb_data in pdbs2operate)

    # #
    # do it for the first conformer
    name, pdb_data = next(pdbs2operate)

    atc_mask, pure_radii_sum, residue_distances = \
        preparation(name, pdb_data)

    validate_conformer_from_preparation(
        name,
        pdb_data,
        atc_mask,
        pure_radii_sum,
        residue_distances,
        residues_apart,
        )
    # done

    execute = partial(
        report_on_crash,
        validate_conformer_from_preparation,
        atc_mask=atc_mask,
        pure_radii_sum=pure_radii_sum,
        residues_apart=residues_apart,
        residue_distances=residue_distances,
        ROC_prefix=_name,
        )

    execute_pool = pool_function(
        partial(starunpack, execute),
        pdbs2operate,
        ncores=ncores,
        )

    for _i in execute_pool:
        pass


def preparation(
        name,
        pdb_data,
        ele_con=('C',),
        ):
    """."""
    s = Structure(pdb_data)
    s.build()
    da = s.data_array
    atom_names = da[:, col_name]
    atom_elements = atom_names.astype('<U1')
    res_numbers = da[:, col_resSeq].astype(np.int)

    return vdW_clash_common_preparation(
        atom_elements,
        elements_to_consider=ele_con,
        residue_numbers=res_numbers,
        )

def validate_conformer_from_preparation(
        name,
        pdb_data,
        atc_mask,
        pure_radii_sum,
        residue_distances,
        residues_apart,
        ):
    assert residue_distances.size

    s = Structure(pdb_data)
    s.build()
    coords = s.coords

    return vdW_clash_calc(
        coords,
        atc_mask=atc_mask,
        pure_radii_sum=pure_radii_sum,
        residue_distances=residue_distances,
        residues_apart=residues_apart)


if __name__ == '__main__':
    libcli.maincli(ap, main)
