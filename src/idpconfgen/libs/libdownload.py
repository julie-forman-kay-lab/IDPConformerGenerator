"""Functions and variables to download files and data."""
import tarfile
import time
import urllib.request
import traceback
from urllib.error import URLError

from idpconfgen import Path, log
from idpconfgen.libs.libio import make_destination_folder, save_file_to_tar
from idpconfgen.libs.libmulticore import pool_function_in_chunks, pool_function, consume_iterable_in_list
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
    for chunk in pool_function_in_chunks(consume_iterable_in_list, items, func, **kwargs):
        for result in chunk:
            for fname, data in result:
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
    dests = str(destination)
    with tarfile.open(dests, mode=_exists[destination.exists()]) as tar:
        for chunk in pool_function_in_chunks(consume_iterable_in_list, items, func, **kwargs):
            for result in chunk:
                for fout, data in result:
                    save_file_to_tar(tar, fout, data)


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
    #try:
    yield from save_structure_chains_and_segments(
            downloaded_data,
            pdbname,
            chains=chains,
            **kwargs,
            )
    #except Exception as err:
    #    log.debug(traceback.format_exc())
    #    log.error(err)


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
        raise
    yield f'{pdbname}.pdb', downloaded_data.decode('utf-8')
