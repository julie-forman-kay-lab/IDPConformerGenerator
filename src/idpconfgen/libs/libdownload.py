import contextlib
import functools
import re
import string
import traceback
import urllib.request
from abc import ABC, abstractmethod
from collections import namedtuple
from multiprocessing.pool import ThreadPool

import numpy as np

from idpconfgen import Path, log
from idpconfgen.core import count_string_formatters
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs import libtimer
from idpconfgen.logger import S, T
from idpconfgen.libs.libstructure import Structure
from idpconfgen.libs import libpdb


PDB_WEB_LINK = "https://files.rcsb.org/download/{}.pdb"
CIF_WEB_LINK = "https://files.rcsb.org/download/{}.cif"
POSSIBLELINKS = [
    PDB_WEB_LINK,
    CIF_WEB_LINK,
    ]


class PDBDownloader:
    """
    Control PDB downloading operations.

    Given a list of :class:`PDBID` downloads those PDB files.

    Parameters
    ----------
    pdb_list : list
        List of PDBIDs.

    destination : str or Path
        Destination folder.

    record_name : tuple
        A tuple containing the atom record names to store from the
        downloaded PDBs. Record names are as described for the PDB
        format v3. Usually 'ATOM' and 'HETATM' can be used.
        Defaults to ('ATOM',)

    ncores : int, optional
        Number of cores to use during the download phase.
        All downloads operations are stored in a pool, each core
        performs a download operation per time grabing those operations
        from the pool. The more cores the more operations are performed
        simultaneoursly.
        Defaults to 1.
    """

    def __init__(
            self,
            pdb_list,
            destination,
            ncores=1,
            record_name=('ATOM',),
            **kwargs,
            ):

        self.pdb_list = pdb_list

        self.destination = Path(destination)
        self.destination.mkdir(parents=True, exist_ok=True)

        self.ncores = ncores

        self.record_name = record_name

        self.kwargs = kwargs

        self.prepare_pdb_list()

        return

    def prepare_pdb_list(self):
        """Prepare a list with the PDBs to download."""
        self.pdbs_to_download = {}
        for pdbid in self.pdb_list:
            pdbentry = self.pdbs_to_download.setdefault(pdbid.name, [])
            pdbentry.append(pdbid.chain)

    @libtimer.record_time()
    def run(self):
        """Run download operation."""
        log.info(T('starting raw PDB download'))
        log.info(S(
            f'{len(self.pdbs_to_download)} PDBs will be downloaded '
            f'and at least {len(self.pdb_list)} chains will be saved'
            ))

        results = ThreadPool(self.ncores).imap_unordered(
            self._download_single_pdb,
            self.pdbs_to_download.items(),
            )

        for _r in results:
            continue

    def _download_single_pdb(self, pdbid_and_chains_tuple):

        pdbname = pdbid_and_chains_tuple[0]
        chains = pdbid_and_chains_tuple[1]

        possible_links = [l.format(pdbname) for l in POSSIBLELINKS]

        with self._attempt_download(pdbname):
            response = self._download_data(possible_links)

        try:
            downloaded_data = response.read()
        except (AttributeError, UnboundLocalError):  # response is None
            return

        pdbdata = Structure(downloaded_data)

        pdbdata.build()

        if chains[0] is None:
            chains = pdbdata.chain_set

        for chain in chains:

            pdbdata.add_filter_record_name(self.record_name)
            pdbdata.add_filter_chain(chain)
            pdbdata.add_filter(lambda x: x[libpdb.atom_altLoc.col] in ('A', ''))
            destination = Path(self.destination, f'{pdbname}_{chain}.pdb')
            try:
                pdbdata.write_PDB(destination)
            except EXCPTS.EmptyFilterError:
                log.debug(traceback.format_exc())
                log.error(S(f'Empty Filter for:'))
                log.error(S(f'{destination}'))
                log.error(S(f'record_name: {self.record_name}'))
                log.error(S(f'chain filter: {chain}'))
            pdbdata.clear_filters()

    def _download_data(self, possible_links):
        for weblink in possible_links:
            try:
                response = urllib.request.urlopen(weblink)
            except urllib.error.HTTPError:
                log.error(S(f'failed from {weblink}'))
                continue
            else:
                log.info(S(f'completed from {weblink}'))
                return response
        else:
            raise EXCPTS.DownloadFailedError

    @contextlib.contextmanager
    def _attempt_download(self, pdbname):
        try:
            yield
        except EXCPTS.DownloadFailedError as e:
            log.error(S(f'{repr(e)}: FAILED {pdbname}'))
            return

