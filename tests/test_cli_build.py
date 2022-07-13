"""Test cli build."""
import shutil

import pytest

from idpconfgen.cli_build import ap, main
from idpconfgen.libs.libio import extract_from_tar

from .tcommons import idp_db_test


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


@pytest.fixture(params=["fixed", "sampling"])
def bgeo_strategies(request):
    return request.param


@pytest.fixture(params=[True, False])
def disable_sidechains(request):
    return request.param


class TestBuildOptions:

    db_path = extract_from_tar(idp_db_test, ext='.json')[0]
    test_folder = 'test_build_conformers'

    def test_build_option_1(self, bgeo_strategies, disable_sidechains):
        """Test a direct command to build conformers."""
        assert self.db_path.exists()
        main(
            "EGAAGAASS",
            self.db_path,
            custom_sampling=None,
            bgeo_strategy=bgeo_strategies,
            nconfs=10,
            output_folder=self.test_folder,
            disable_sidechain=disable_sidechains,
            )

    def test_remove_folder(self):
        # in case the above test fails
        shutil.rmtree(self.test_folder, ignore_errors=True)
        self.db_path.unlink()
