"""Functions and variables to download files and data."""
import time
import traceback
import urllib.request
from urllib.error import URLError
from contextlib import contextmanager
from os import SEEK_END
from io import BytesIO

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
                    #log.error(S(f'Download {pdbid} failed to read data.'))
                    continue
            else:
                raise IOError(f'Failed to download {pdbid}')
        except (TimeoutError, URLError):
            log.error(f'failed download for {pdbid} because of TimeoutError. Retrying...')
            time.sleep(15)
            attempts += 1
    else:
        raise IOError(f'Failed to download {pdbid} - too much attempts')



def get_pdbs_downloader(destination):
    """Get proper function to download PDBs based on the destination type."""
    # classes to manage the download action
    DOWNLOADER = {
        True: download_pdbs_to_folder,  # this is the equivalent to the else statement
        destination.suffix == '.tar': download_pdbs_to_tar,
        }
    return DOWNLOADER[True]



def download_pdbs_to_folder(destination, items, **kwargs):
    """
    Download PDBs to folder.

    Uses :func:`idpconfgen.libs.libmulticore.pool_function_in_chunks`
    """
    dest = make_destination_folder(destination)

    for result_dict in pool_function_in_chunks(
            download_structure,
            items,
            **kwargs,
            ):

        for fname, data in sorted(result_dict.items()):
            with open(Path(self.dest, fname), 'w') as fout:
                fout.write('\n'.join(data))



def download_pdbs_to_tar(destination, items, **kwargs):
    """
    Download PDBs to tarfile.

    Uses :func:`idpconfgen.libs.libmulticore.pool_function_in_chunks`
    """
    _exists = {True: 'a', False: 'w'}
    dests = destination.str()
    dest = tarfile.open(dests, mode=_exists[destination.exists()])
    dest.close()

    for result_dict in pool_function_in_chunks(
            download_structure,
            items,
            **kwargs,
            ):

        tar = tarfile.open(dests, mode='a:')

        for fout, _data in sorted(result_dict.items()):
            try:
                sIO = BytesIO()
                sIO.write('\n'.join(_data).encode())
                info = tarfile.TarInfo(name=fout)
                info.size=sIO.seek(0, SEEK_END)
                sIO.seek(0)
                tar.addfile(tarinfo=info, fileobj=sIO)
            except Exception:
                log.error(f'failed for {fout}')
        tar.close()


