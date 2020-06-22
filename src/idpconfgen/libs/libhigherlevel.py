"""
Higher level functions.

Function which operate with several libraries
and are defined here to avoid circular imports.
"""
from collections import defaultdict
from contextlib import suppress
from functools import reduce, partial

from idpconfgen import Path, log
from idpconfgen.core.definitions import blocked_ids
from idpconfgen.libs.libio import concatenate_entries, read_PDBID_from_source, save_pairs_dispacher
from idpconfgen.libs.libparse import group_runs, group_by
from idpconfgen.libs.libpdb import PDBList, PDBIDFactory, atom_resSeq, atom_name, atom_resName, atom_resSeq
from idpconfgen.libs.libstructure import Structure, col_resSeq
from idpconfgen.libs.libtimer import record_time
from idpconfgen.logger import S, T, init_files
from idpconfgen.libs.libmulticore import pool_function_in_chunks, consume_iterable_in_list


# USED OKAY!
def download_pipeline(func, logfilename='.download'):
    """
    Context pipeline to download PDB/mmCIF files.

    Exists because fetching and download filtered PDBs shared the same
    operational contextm, only differing on the function which
    orchestrates the download process.

    Parameters
    ----------
    func : function
        The actual function that orchestrates the download operation.

    logfilename : str, optional
        The common stem name of the log files.
    """
    LOGFILESNAME = logfilename
    def main(
            pdbids,
            chunks=5_000,
            destination=None,
            ncores=1,
            update=False,
            ):
        """Run main script logic."""
        init_files(log, LOGFILESNAME)

        #
        log.info(T('reading input PDB list'))

        pdblist = PDBList(concatenate_entries(pdbids))

        log.info(
            f"{S(str(pdblist))}\n"
            f"{S('done')}\n"
            )

        #
        log.info(T('Filtering input'))
        destination = destination or Path.cwd()
        log.info(
            f"{S(f'from destination: {destination}')}\n"
            f"{S('and other sources...')}"
            )

        # comparison block
        def diff(first, other):
            return first.difference(other)

        remove_from_input = [
            read_PDBID_from_source(destination),
            PDBList(blocked_ids),
            ]

        # yes, there are just two items in remove_from_input, why use reduce?
        # what if more are added in the future? :-P the engine is already created
        pdblist_comparison = reduce(diff, remove_from_input, pdblist)
        log.info(S(f'Found {str(pdblist_comparison)} to download'))
        #

        something_to_download = len(pdblist_comparison) > 0
        if something_to_download and update:

            execute = partial(
                pool_function_in_chunks,
                consume_iterable_in_list,
                list(pdblist_comparison.name_chains_dict.items()),
                func,
                ncores=ncores,
                chunks=chunks,
                )

            save_pairs = save_pairs_dispacher[destination]

            for chunk in execute():
                # flatted to avoid open/close tar file every loop iteration
                flatted = (item for result in chunk for item in result)
                save_pairs(flatted, destination)

            log.info(T('Reading UPDATED destination'))
            pdblist_updated = read_PDBID_from_source(destination)
            pdblist_up_comparison = pdblist.difference(pdblist_updated)
            log.info(S(f'{str(pdblist_up_comparison)}'))
            if len(pdblist_up_comparison) > 0:
                log.info(S(
                    'There are PDBIDs not downloaded\n.'
                    'Those IDs have been registered in the '
                    f'{LOGFILESNAME}.debug file.'
                    ))
                log.debug('\n'.join(str(_id) for _id in pdblist_up_comparison))

        elif not something_to_download and update:
            log.info(S('There is nothing to download.'))
            log.info(S('All requested IDs are already at the destination folder.'))

        log.info(T('PDB Downloader finished'))
        return

    return main


# USED OKAY
def extract_secondary_structure(
        pdbid,
        ssdata,
        atoms='all',
        minimum=0,
        structure='all',
        ):
    """
    Extract secondary structure elements from PDB data.

    Parameters
    ----------
    pdbid : tuple
        Where index 0 is the PDB id code, for example `12AS.pdb`, or
        `12AS_A`, or `12AS_A_seg1.pdb`.
        And, index 1 is the PDB data itself in bytes.

    ssdata : dict
        Dictionary containing the DSSP information. Must contain a key
        equal to `Path(pdbid).stem, where `dssp` key contains the DSSP
        information for that PDB.

    atoms : str or list of str or bytes, optional
        The atom names to keep.
        Defaults to `all`.

    minimum : int, optional
        The minimum size a segment must have in order to be considered.
        Defaults to 0.

    structure : str or list of chars
        The secondary structure character to separate.
        Multiple can be given in the form of a list.
    """
    pdbname = Path(pdbid[0]).stem
    pdbdata = pdbid[1].split(b'\n')
    pdbdd = ssdata[pdbname]

    ss_identified = set(pdbdd['dssp'])
    if structure == 'all':
        ss_to_isolate = ss_identified
    else:
        ss_to_isolate = set(s for s in ss_identified if s in structure)

    # general lines filters
    line_filters = []
    LF_append = line_filters.append
    LF_pop = line_filters.pop

    # in atoms filter
    if atoms != 'all':
        with suppress(AttributeError):  # it is more common to receive str
            atoms = [c.encode() for c in atoms]
        line_filters.append(lambda x: x[atom_name].strip() in atoms)

    dssp_slices = group_by(pdbdd['dssp'])
    # DR stands for dssp residues
    DR = [c.encode() for c in pdbdd['resids'].split(',')]

    for ss in ss_to_isolate:

        ssfilter = (slice_ for char, slice_ in dssp_slices if char == ss)
        minimum_size = (s for s in ssfilter if s.stop - s.start >= minimum)

        for counter, seg_slice in enumerate(minimum_size):

            LF_append(lambda x: x[atom_resSeq].strip() in DR[seg_slice])
            pdb = b'\n'.join(
                l for l in pdbdata if all(f(l) for f in line_filters)
                )
            LF_pop()

            yield f'{pdbname}_{ss}_{counter}.pdb', pdb

            counter += 1
