"""Test cli build."""
import pytest

from idpconfgen.cli_build import ap


@pytest.mark.parametrize(
    'command',
    ['dhelix', 'dstrand', 'dany'],
    )
def test_ap_dsecondary_structure_true(command):
    """Test dloop argument."""
    cmd = ap.parse_args(f'-db dummy.json -seq AAAAA --{command}'.split())
    d = vars(cmd)
    assert d[command] is True


@pytest.mark.parametrize(
    'command',
    ['dhelix', 'dstrand', 'dany'],
    )
def test_ap_dsecondary_structure_false(command):
    """Test dloop argument."""
    cmd = ap.parse_args('-db dummy.json -seq AAAAA'.split())
    d = vars(cmd)
    assert d[command] is False


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
