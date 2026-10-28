"""Test cli build."""
import shutil

import numpy as np
import pytest

from idpconfgen.cli_build import ap, main
from idpconfgen.libs.libio import extract_from_tar
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resSeq,
    cols_coords,
    )

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


def _angle(a, b, c):
    u, v = a - b, c - b
    return np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))


def test_build_long(tmp_path):
    """Test --long joins its fragments with peptide-bond geometry."""
    db_path = extract_from_tar(idp_db_test, ext='.json', output=tmp_path)[0]
    seq = 'EGAAGAASS' * 12
    main(
        seq,
        db_path,
        custom_sampling=None,
        long=True,
        force_long=True,
        dany=True,
        dloop_off=True,
        nconfs=2,
        output_folder=tmp_path / 'long',
        disable_sidechains=True,
        )
    files = sorted((tmp_path / 'long').glob('conformer_*.pdb'))
    assert len(files) == 2
    for path in files:
        structure = Structure(path)
        structure.build()
        rows = structure.data_array
        resseq = rows[:, col_resSeq].astype(int)
        assert sorted(set(resseq)) == list(range(1, len(seq) + 1))

        def xyz(res, name):
            sel = (resseq == res) & (rows[:, col_name] == name)
            return rows[sel][0, cols_coords].astype(float)

        # fragments of 50 residues meet at 49-50 and 99-100
        for i in (49, 99):
            c, n = xyz(i, 'C'), xyz(i + 1, 'N')
            assert abs(np.linalg.norm(c - n) - 1.33) < 0.05
            assert abs(_angle(xyz(i, 'CA'), c, n) - 116.2) < 10
            assert abs(_angle(xyz(i, 'O'), c, n) - 122.7) < 10
            assert abs(_angle(c, n, xyz(i + 1, 'CA')) - 121.7) < 10
