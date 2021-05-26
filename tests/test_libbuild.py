"Test cli_build."""
import numpy as np
import pytest

from idpconfgen.libs.libbuild import init_confmasks, build_regex_substitutions


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
        ('ASDS', {'S': 'SE'}, 'A[SE]D[SE]'),
        ('ASD', {'D': 'DEWQ'}, 'AS[DEWQ]'),
        ('ASD', {'R': 'REWQ'}, 'ASD'),
        ]
    )
def test_add_regex_substitutions(in1, options, expected):
    """Test add regex substitutions."""
    result = build_regex_substitutions(in1, options)
    assert result == expected
