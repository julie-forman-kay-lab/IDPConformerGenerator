"""Functions and variables to download files and data."""
import tarfile
import time
import urllib.request
from urllib.error import URLError

from idpconfgen import Path, log
from idpconfgen.libs.libio import make_destination_folder, save_file_to_tar
from idpconfgen.libs.libmulticore import pool_function_in_chunks
from idpconfgen.libs.libstructure import save_structure_chains_and_segments
from idpconfgen.logger import S


PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]


def download_dispacher(func, destination, *args, **kwargs):
    """Dispaches the appropriate download env based on `destination`."""
    DOWNLOADER = {
        True: download_pdbs_to_folder,
        destination.suffix == '.tar': download_pdbs_to_tar,
        }
    return DOWNLOADER[True](destination, *args, func=func, **kwargs)


def download_pdbs_to_folder(destination, items, func=None, **kwargs):
    """
    Download PDBs to folder.

    Uses :func:`idpconfgen.libs.libmulticore.pool_function_in_chunks`
    """
    dest = make_destination_folder(destination)
    for result_dict in pool_function_in_chunks(func, items, **kwargs):
        for fname, data in result_dict.items():
            with open(Path(dest, fname), 'w') as fout:
                fout.write(data)


def download_pdbs_to_tar(destination, items, func=None, **kwargs):
    """
    Download PDBs to tarfile.

    Uses :func:`idpconfgen.libs.libmulticore.pool_function_in_chunks`
    """
    # append 'a' here combos with the
    # read_PDBID_from_source(destination), in cli_pdbdownloader
    _exists = {True: 'a', False: 'w'}
    dests = destination.str()
    dest = tarfile.open(dests, mode=_exists[destination.exists()])
    dest.close()

    for result_dict in pool_function_in_chunks(func, items, **kwargs):
        with tarfile.open(dests, mode='a:') as tar:
            for fout, _data in sorted(result_dict.items()):
                try:
                    save_file_to_tar(tar, fout, _data)
                except Exception:
                    log.error(f'failed for {fout}')


def download_structure(pdbid, **kwargs):
    """
    Download a PDB/CIF structure chains.

    Parameters
    ----------
    pdbid : tuple of 2-elements
        0-indexed, the structure ID at RCSB.org;
        1-indexed, a list of the chains to download.

    **kwargs : as for :func:`save_structure_chains_and_segments`.
    """
    pdbname = pdbid[0]
    chains = pdbid[1]

    try:
        downloaded_data = fetch_pdb_id_from_RCSB(pdbname)
    except IOError:
        log.error(f'Complete download failure for {pdbname}')
        return
    save_structure_chains_and_segments(
        downloaded_data,
        pdbname,
        chains=chains,
        **kwargs,
        )


def fetch_pdb_id_from_RCSB(pdbid):
    """Fetch PDBID from RCSB."""
    possible_links = (l.format(pdbid) for l in POSSIBLELINKS)

    attempts = 0
    while attempts < 10:
        try:
            for weblink in possible_links:
                try:
                    response = urllib.request.urlopen(weblink)
                    return response.read()
                except urllib.error.HTTPError:
                    continue
                except (AttributeError, UnboundLocalError):  # response is None
                    log.debug(S(f'Download {weblink} failed.'))
                    continue
            else:
                raise IOError(f'Failed to download {pdbid}')
        except (TimeoutError, URLError):
            log.error(
                f'failed download for {pdbid} because of TimeoutError. '
                'Retrying...'
                )
            time.sleep(15)
            attempts += 1
    else:
        raise IOError(f'Failed to download {pdbid} - attempts exhausted')


def fetch_raw_PDBs(pdbid, mdict=None, **kwargs):
    """Download raw PDBs without any filtering."""
    pdbname = pdbid[0]
    try:
        downloaded_data = fetch_pdb_id_from_RCSB(pdbname)
    except IOError:
        log.error(f'Complete download failure for {pdbname}')
    mdict[f'{pdbname}.pdb'] = downloaded_data.decode('utf-8')
