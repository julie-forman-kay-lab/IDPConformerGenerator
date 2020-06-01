"""Functions and variables to download files and data."""
import traceback
import urllib.request
import urllib.error.URLError as URLError
from contextlib import contextmanager

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libstructure import Structure, col_altLoc, write_PDB, structure_to_pdb, col_resSeq, col_resName, col_element
from idpconfgen.libs.libparse import save_structure_chains_and_segments
from idpconfgen.logger import S


PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]


def download_raw_PDBS(pdbid, folder=''):
    """
    Download raw PDBs without any filtering.
    """
    downloaded_data = fetch_pdb_id_from_RCSB(pdbid.name)
    s = Structure(downloaded_data)
    s.build()
    s.add_filter_chain(pdbid.chain)
    fout = Path(folder, f'{pdbid}.pdb')
    fout.write_bytes(downloaded_data)


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

    downloaded_data = fetch_pdb_id_from_RCSB(pdbname)
    save_structure_chains_and_segments(
        downloaded_data,
        pdbname,
        chains=chains,
        **kwargs,
        )



def fetch_pdb_id_from_RCSB(pdbid):
    possible_links = (l.format(pdbid) for l in POSSIBLELINKS)

    while attempts < 10:
        try:
            for weblink in possible_links:
                try:
                    response = urllib.request.urlopen(weblink)
                    return response.read()
                except urllib.error.HTTPError:
                    continue
                except (AttributeError, UnboundLocalError):  # response is None
                    #log.error(S(f'Download {pdbid} failed to read data.'))
                    continue
            else:
                raise IOError(f'Failed to download {pdbid}')
        except (TimeoutError, URLError):
            log.error(f'failed downlaod for {pdbid} because of TimeoutError. Retrying...')
            time.sleep(15)
            attempts += 1


