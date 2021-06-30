"Test cli_build."""
import numpy as np
import pytest

from idpconfgen.libs.libbuild import *


def test_conf_masks():
    """Test confmask creation."""

    atom_labels = np.array([
        #0     1     2
        'H1', 'H2', 'H3',
        #3    4     5    6
        'N', 'CA', 'C', 'O',
        #7     8      9     10   11
        'CB', 'H', 'HB1', 'HB2', 'HB3',
        #12   13   14    15
        'N', 'CA', 'C', 'O',
        #16     17   18
        'H', 'HA2', 'HA3',
        #19
        'OXT',
        ])

    assert len(atom_labels) == 20

    confmask = init_confmasks(atom_labels)

    assert np.array_equal(confmask.bb3, [3, 4, 5, 12, 13, 14])
    assert np.array_equal(confmask.bb4, [3, 4, 5, 6, 12, 13, 14, 15])
    assert np.array_equal(confmask.NHs, [8, 16])
    assert np.array_equal(confmask.COs, [6])
    assert confmask.OXT1 == 15
    assert confmask.OXT2 == 19
    assert np.array_equal(confmask.non_Hs, [3, 4, 5, 6, 7, 12, 13, 14, 15, 19])
    assert np.array_equal(confmask.non_Hs_non_OXT, [3, 4, 5, 6, 7, 12, 13, 14, 15])
    assert np.array_equal(confmask.H1_N_CA_CB, [0, 3, 4, 7])



@pytest.mark.parametrize(
    'in1,options,expected',
    [
        ('ASD', {'S': 'SE'}, 'A[SE]D'),
        ('S', {'S': 'SE'}, '[SE]'),
        ('ASDS', {'S': 'SE'}, 'A[SE]D[SE]'),
        ('ASDSE', {'S': 'SE', 'E': 'DE'}, 'A[SE]D[SE][DE]'),
        ('ASD', {'D': 'DEWQ'}, 'AS[DEWQ]'),
        ('ASD', {'R': 'REWQ'}, 'ASD'),
        ('ASD', {}, 'ASD'),
        ]
    )
def test_add_regex_substitutions(in1, options, expected):
    """Test add regex substitutions."""
    result = build_regex_substitutions(in1, options)
    assert result == expected


@pytest.mark.parametrize(
    'primary,input_seq,mers_size,res_tolerance,expected',
    [
        (
            'AADDAADDAADD',
            'AD',
            (2,),
            None,
            {2: {"AD": [slice(1, 3), slice(5, 7), slice(9, 11)]}},
            ),
        ],
    )
def test_prepare_slice_dict(
        primary,
        input_seq,
        mers_size,
        res_tolerance,
        expected):
    """Test prepare slice dict."""
    result = prepare_slice_dict(primary, input_seq, mers_size, res_tolerance)
    assert result.keys() == expected.keys()
    assert result == expected
