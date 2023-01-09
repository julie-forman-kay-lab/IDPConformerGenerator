"""
Client for renaming 3 letter residues in a PDB file.

USAGE:
    $ idpconfgen resptm <PDB-FILES> -pt <PATTERN>
    $ idpconfgen resptm <PDB-FILES> -pt <PATTERN> -o <OUTPUT> -n <CORES>
"""
import argparse
import shutil
import re
from functools import partial
from pathlib import Path

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import extract_from_tar, read_path_bundle
from idpconfgen.libs.libmulticore import pool_function
from idpconfgen.logger import S, T, init_files, report_on_crash

from idpconfgen.libs.libio import make_folder_or_cwd
from idpconfgen.libs.libstructure import (
    Structure,
    structure_to_pdb,
    col_resName,
    col_resSeq,
)

LOGFILESNAME = '.idpconfgen_resptm'
TMPDIR = '__tmpresptm__'

_name = 'resptm'
_help = 'Renames residues of interest within a PDB file for post-translational modifications.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)

ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_argument_pdb_files(ap)

ap.add_argument(
    '-pt',
    '--pattern',
    help=(
        "Pattern is denoted by colons for residue number and desired "
        "new residue name. Additional residues are delimited by commas "
        "and pattern must end at a comma. For e.g. -pt 12:TPO,14:SEP,"
        ),
    nargs='?',
    )

libcli.add_argument_output_folder(ap)
libcli.add_argument_ncores(ap)

ap.add_argument(
    '-tmpdir',
    help=(
        'Temporary directory to store data during calculation '
        'if needed.'
        ),
    type=Path,
    default=TMPDIR,
    )


def rename_residues(target_struc, pattern):
    """
    Renames 3 letter residue codes in a PDB based on the pattern.

    Parameters
    ----------
    target_struc : Path
        Path to PDB file to edit.
        
    pattern : dict of str
        Dictionary of patterns where the first key-layer is the residue number
        and values are 3-letter code to be changed.
    
    Returns
    -------
    new_struc_arr : np.ndarray (N x 16)
        Array format of new re-named PDB file.
    """
    struc = Structure(target_struc)
    struc.build()
    
    new_struc_array = struc.data_array
    for target_res in pattern:
        for i, res in enumerate(new_struc_array[:, col_resSeq]):
            if res == target_res:
                new_struc_array[:, col_resName][i] = pattern[target_res]
            elif res == str(int(target_res) + 1):
                break
            
    return new_struc_array


def main(
        pdb_files,
        pattern,
        output_folder=None,
        ncores=1,
        tmpdir=TMPDIR,
        **kwargs,
        ):
    output_folder = make_folder_or_cwd(output_folder)
    init_files(log, Path(output_folder, LOGFILESNAME))
    
    REGEX_PATTERN = re.compile(r"[0-9]*:[A-Z]{3},")
    RES_COMBO = {}
    
    log.info(T('reading input paths'))
    try:
        pdbs2operate = extract_from_tar(pdb_files, output=tmpdir, ext='.pdb')
        _istarfile = True
    except (OSError, TypeError):
        pdbs2operate = list(read_path_bundle(pdb_files, ext='pdb'))
        _istarfile = False
    log.info(S('done'))
    
    log.info(T('reading input pattern'))
    if REGEX_PATTERN.match(pattern):
        residues = pattern.split(',')
        residues.pop()  # last element should be empty
        for res in residues:
            tmp = res.split(':')
            RES_COMBO[tmp[0]] = tmp[1]
            log.info(S(f'Residue {tmp[0]} will be changed to {tmp[1]}.'))
    else:
        log.info(S('Incorrect pattern input.'))
        log.info(S('Pattern is as follows: 1:ABC,...'))
        return
    log.info(S('done'))
    
    consume = partial(
        rename_residues,
        pattern=RES_COMBO,
    )
    execute = partial(
        report_on_crash,
        consume,
        ROC_exception=Exception,
        ROC_folder=output_folder,
        ROC_prefix=_name
    )
    execute_pool = pool_function(execute, pdbs2operate, ncores=ncores)
    
    for i, conf in enumerate(execute_pool):
        struc = structure_to_pdb(conf)
        output = output_folder.joinpath(f"resptm_conf_{i + 1}.pdb")
        with open(output, 'w') as f:
            for line in struc:
                f.write(line + "\n")
    
    if _istarfile:
        shutil.rmtree(tmpdir)


if __name__ == '__main__':
    libcli.maincli(ap, main)
