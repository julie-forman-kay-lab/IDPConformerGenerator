"""Functions and variables to download files and data."""
import time
import urllib.request
from functools import partial
from urllib.error import URLError

from idpconfgen import log
from idpconfgen.core.exceptions import DownloadFailedError
from idpconfgen.libs.libstructure import save_structure_by_chains
from idpconfgen.logger import S


PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]


def download_structure(pdbid, mmcif=False, **kwargs):
    """
    Download a PDB/CIF structure chains.

    Parameters
    ----------
    pdbid : tuple of 2-elements
        0-indexed, the structure ID at RCSB.org;
        1-indexed, a list of the chains to download.

    **kwargs : as for :func:`save_structure_by_chains`.
    """
    pdbname = pdbid[0]
    chains = pdbid[1]

    downloaded_data, _ = fetch_pdb_id_from_RCSB(pdbname, mmcif=mmcif)

    # save_structure_by_chains is a specific function that always saves
    # in PDB format.
    yield from save_structure_by_chains(
        downloaded_data,
        pdbname,
        chains=chains,
        **kwargs,
        )


def fetch_pdb_id_from_RCSB(pdbid, mmcif=False):
    """
    Fetch PDBID from RCSB.

    Returns
    -------
    tuple (str, str)
        urllib.request.urlopen.response.read()
        PDB file extension (.pdb, .cif, ...)
    """
    if mmcif:
        POSSIBLELINKS.reverse()

    possible_links = (link.format(pdbid) for link in POSSIBLELINKS)

    attempts = 0
    while attempts < 10:
        try:
            for weblink in possible_links:
                try:
                    response = urllib.request.urlopen(weblink)
                    return response.read(), weblink.rpartition(".")[-1]
                except urllib.error.HTTPError:
                    continue
                except (AttributeError, UnboundLocalError):  # response is None
                    log.debug(S(f'Download {weblink} failed.'))
                    continue
            else:
                break
        except (TimeoutError, URLError) as err:
            log.error(
                f'failed download for {pdbid} because of {repr(err)}. '
                'Retrying...'
                )
            time.sleep(15)
            attempts += 1
    else:
        raise DownloadFailedError(f'Failed to download {pdbid}')


def fetch_raw_structure(pdbid, ext, **kwargs):
    """Download raw structure from RCSB without any filtering."""
    pdbname = pdbid[0]
    mmcif = True if ext == "cif" else False
    downloaded_data, ext = fetch_pdb_id_from_RCSB(pdbname, mmcif)
    yield f'{pdbname}.{ext}', downloaded_data.decode('utf-8')


fetch_raw_PDBs = partial(fetch_raw_structure, ext='pdb')
fetch_raw_CIFs = partial(fetch_raw_structure, ext='cif')
