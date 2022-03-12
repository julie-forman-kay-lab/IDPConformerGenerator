"""Test cli build."""
import pytest

from idpconfgen.cli_build import ap


@pytest.fixture(
    params=  # noqa: E251, E131
        [  # noqa: E251, E131
            (1, 5),
            (2, 10),
            (5, 40),
            ]
    )
def dtuples_good(request):
    """Tuples that should work in dloop/dhelix/dstrand parameters."""
    return request.param


@pytest.fixture(
    params=  # noqa: E251, E131
        [  # noqa: E251, E131
            (1, -5),
            (0, 5),
            (-2, 10),
            (-5, -40),
            ('a', -40),
            (1, 'b'),
            ('a', 'b'),
            (1.3, 4.6),
            ]
    )
def dtuples_bad(request):
    """Tuples that should break in dloop/dhelix/dstrand parameters."""
    return request.param


@pytest.mark.parametrize(
    'command',
    ['dloop', 'dhelix', 'dstrand'],
    )
def test_ap_dsecondary_structure(command, dtuples_good):
    """Test dloop argument."""
    one, two = dtuples_good
    cmd = ap.parse_args(f'-db dummy.json -seq AAAAA --{command} {one} {two}'.split())  # noqa: E501
    d = vars(cmd)
    assert d[command] == [one, two]


@pytest.mark.parametrize(
    'command',
    ['dloop', 'dhelix', 'dstrand'],
    )
def test_ap_dSS_with_error(command, dtuples_bad):
    """Test dloop argument error on negative numbers."""
    one, two = dtuples_bad
    with pytest.raises(SystemExit) as err:
        ap.parse_args(f'-db dummy.json -seq AAAAA --{command} {one} {two}'.split())  # noqa: E501
    assert err.value.code == 2


@pytest.mark.parametrize(
    'command,one,two',
    [
        ('duser', "L+", "H+"),
        ],
    )
def test_duser_argument(command, one, two):
    """Test dloop argument."""
    cmd = ap.parse_args(f'-db dummy.json -seq AAAAA --{command} {one} {two}'.split())  # noqa: E501
    d = vars(cmd)
    assert d[command] == [one, two]
