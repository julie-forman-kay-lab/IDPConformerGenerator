"""Functions and variables to download files and data."""
import traceback
import urllib.request
from contextlib import contextmanager

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libstructure import Structure, col_altLoc
from idpconfgen.logger import S


PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]


def download_structure(pdbid, folder='', record_name='ATOM'):
    """
    Download a PDB/CIF structure.

    Parameters
    ----------
    pdbid : tuple of 2-elements
        0-indexed, the structure ID at RCSB.org;
        1-indexed, a list of the chains to download.

    folder : str or Path, optional
        The destination folder. Default to current working directory.

    record_name : str or tuple
        The RECORD name lines to save.
    """
    pdbname = pdbid[0]
    chains = pdbid[1]

    possible_links = (l.format(pdbname) for l in POSSIBLELINKS)

    for weblink in possible_links:
        try:
            response = urllib.request.urlopen(weblink)
            downloaded_data = response.read()
            break
        except urllib.error.HTTPError:
            continue
        except (AttributeError, UnboundLocalError):  # response is None
            log.error(S(f'Download {pdbname} failed to read data.'))
            continue
    else:
        log.error(S(f'Failed to download {pdbname}'))
        return

    pdbdata = Structure(downloaded_data)
    pdbdata.build()

    chains = chains or pdbdata.chain_set

    pdbdata.add_filter_record_name(record_name)
    pdbdata.add_filter(
        lambda x: x[col_altLoc] in ('A', '', ' ')
        )

    for chain in chains:
        pdbdata.add_filter_chain(chain)
        fout = Path(folder, f'{pdbname}_{chain}.pdb')
        with try_to_write(downloaded_data, fout):
            pdbdata.write_PDB(fout)
        pdbdata.pop_last_filter()


@contextmanager
def try_to_write(data, fout):
    """Context to download."""
    try:
        yield
    except EXCPTS.IDPConfGenException as err:
        log.debug(traceback.format_exc())
        log.error(S('error found for {}: {}', fout, repr(err)))
        erp = Path(fout.myparents(), 'errors')
        erp.mkdir(parents=True, exist_ok=True)
        p = Path(erp, fout.stem).with_suffix('.structure')
        try:
            p.write_bytes(data)
        except TypeError:
            p.write_text(data)
