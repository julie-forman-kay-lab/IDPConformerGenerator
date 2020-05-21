"""Functions and variables to download files and data."""
import traceback
import urllib.request
from contextlib import contextmanager

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libstructure import Structure, col_altLoc, write_PDB, structure_to_pdb, col_resSeq, col_resName, col_element
from idpconfgen.logger import S


PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]


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


# considers solvent and DNA/RNA
# http://www.wwpdb.org/documentation/file-format-content/format33/sect4.html#HET
_discarded_residues = ('HOH', 'SOL', 'I', 'C', 'G', 'A',
    'U', 'I', 'DC', 'DG', 'DA', 'DU', 'DT', 'DI', 'N')

_allowed_elements = ('C', 'O', 'N', 'H', 'S', 'Se', 'D')


def save_structure_chains_and_segments(
        downloaded_data,
        pdbname,
        chains=None,
        record_name='ATOM',
        altlocs=('A', '', ' '),
        folder='',
        ):

    _DR = _discarded_residues
    _AE = _allowed_elements

    pdbdata = Structure(downloaded_data)
    pdbdata.build()

    chains = chains or pdbdata.chain_set

    pdbdata.add_filter_record_name(record_name)
    pdbdata.add_filter(lambda x: x[col_resName] not in _DR)
    pdbdata.add_filter(lambda x: x[col_element] in _AE)
    pdbdata.add_filter(lambda x: x[col_altLoc] in altlocs)

    for chain in chains:
        pdbdata.add_filter_chain(chain)

        # passar esto a una function
        pdbsegs = pdbdata.residue_segments

        # more than one residue
        valid_segments = filter(
            lambda x: set(line[col_resSeq] for line in x),
            pdbsegs,
            )

        if len(pdbsegs) > 1:
            for i, segment in enumerate(valid_segments):
                fout_seg = Path(folder, f'{pdbname}_{chain}_seg{i}.pdb')
                with try_to_write(downloaded_data, fout_seg):
                    # because segments are pure arrays
                    # and not Structure objects
                    write_PDB(structure_to_pdb(segment), fout_seg)
        #

        else:
            fout = Path(folder, f'{pdbname}_{chain}.pdb')
            with try_to_write(downloaded_data, fout):
                pdbdata.write_PDB(fout)

        pdbdata.pop_last_filter()


def fetch_pdb_id_from_RCSB(pdbid):
    possible_links = (l.format(pdbid) for l in POSSIBLELINKS)

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
        if p.exists():
            return
        else:
            try:
                p.write_bytes(data)
            except TypeError:
                p.write_text(data)
