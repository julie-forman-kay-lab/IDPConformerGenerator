"""Use FASPR to build sidechains."""
from idpconfgen import Path, source_folder


faspr_dun2010bbdep_path = Path(
    source_folder,
    'components',
    'sidechain_packing',
    'dun2010bbdep.bin',
    )


def init_faspr_sidechains(
        input_seq,
        template_masks,
        all_atom_masks,
        **ignore,
        ):
    """
    Instantiate dedicated function environment for FASPR sidehchain calculation.

    Examples
    --------
    >>> calc_faspr = init_fastpr_sidechains('MASFRTPKKLCVAGG', ...)
    >>> # a (N, 3) array with the N,CA,C,O coordinates
    >>> coords = np.array( ... )
    >>> calc_faspr(coords)

    Parameters
    ----------
    input_seq : str
        The FASTA sequence of the protein for which this function will
        be used.

    template_masks : `libbuil.ConfMasks` object
        Related to the building template in `cli_build`.

    all_atom_masks : `libbuil.ConfMasks` object
        Related to the all atoom coordinate system in `cli_build`.

    Returns
    -------
    np.ndarray, dtype=bool, (M, 3)
        Array mask for the builder.

    np.ndarray (M, 3)
        Heavy atom coordinates of the protein sequence.
    """
    # we import idpcpp here to avoid problems with tox
    import idpcpp
    faspr_func = idpcpp.faspr_sidechains
    faspr_dun2010_bbdep_str = str(faspr_dun2010bbdep_path)

    def compute_faspr_sidechains(coords):
        """Do calculation."""
        _coords = faspr_func(
            coords[template_masks.bb4],
            input_seq,
            faspr_dun2010_bbdep_str,
            )

        return all_atom_masks.non_Hs_non_OXT, _coords

    return compute_faspr_sidechains
