"""Test libfilter functions."""

import numpy as np
import pytest

from idpconfgen.libs.libfilter import aligndb

@pytest.mark.parametrize(
    'in1,out',
    [
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5",
                    },
                },
            (
                {"pdb1": slice(0, 6, None)},
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL",
                "MMMMMM",
                ),
            ),
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "HHHHHH",
                    "fasta": "AAAAAA",
                    "omega": [5, 5, 5, 5, 5],
                    "psi": [4, 4, 4, 4, 4],
                    "phi": [6, 6, 6, 6, 6],
                    "resids": "7,8,9,10,11,12",
                    },
                },
            (
                {
                    "pdb1": slice(0, 6, None),
                    "pdb2": slice(7, 13, None),
                    },
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    [np.nan, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL|HHHHHH",
                "MMMMMM|AAAAAA",
                ),
            ),
        (
            {
                "pdb1": {
                    "dssp": "LLLLLL",
                    "fasta": "MMMMMM",
                    "omega": [2, 2, 2, 2, 2],
                    "psi": [1, 1, 1, 1, 1],
                    "phi": [3, 3, 3, 3, 3],
                    "resids": "1,2,3,4,5,6",
                    },
                "pdb2": {
                    "dssp": "HHHHHH",
                    "fasta": "AAAAAA",
                    "omega": [5, 5, 5, 5, 5],
                    "psi": [4, 4, 4, 4, 4],
                    "phi": [6, 6, 6, 6, 6],
                    "resids": "7,8,9,10,11,12",
                    },
                "pdb3": {
                    "dssp": "EEEEEE",
                    "fasta": "WWWWWW",
                    "omega": [8, 8, 8, 8, 8],
                    "psi": [7, 7, 7, 7, 7],
                    "phi": [9, 9, 9, 9, 9],
                    "resids": "13,14,15,16,17,18",
                    },
                },
            (
                {
                    "pdb1": slice(0, 6, None),
                    "pdb2": slice(7, 13, None),
                    "pdb3": slice(14, 20, None),
                    },
                np.array([
                    [np.nan, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, 1, 2],
                    [3, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    [np.nan, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, 4, 5],
                    [6, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    [np.nan, 7, 8],
                    [9, 7, 8],
                    [9, 7, 8],
                    [9, 7, 8],
                    [9, 7, 8],
                    [9, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                    ]),
                "LLLLLL|HHHHHH|EEEEEE",
                "MMMMMM|AAAAAA|WWWWWW",
                ),
            ),
        ],
    )
def test_aligndb(in1, out):
    """."""
    assert isinstance(in1, dict)
    pdbs, angles, dssp, fasta = aligndb(in1)
    assert pdbs == out[0]
    assert angles.shape == out[1].shape
    assert np.array_equal(angles, out[1], equal_nan=True)
    assert dssp == out[2]
    assert fasta == out[3]
