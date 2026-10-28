"""Test the Linker-IDR loop-closure options of cli_ldrs."""
import math
import os

import numpy as np
import pytest

from idpconfgen import cli_ldrs
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.ldrs_helper import split_long_idr
from idpconfgen.libs import libccd
from idpconfgen.libs.libstructure import (
    Structure,
    col_name,
    col_resSeq,
    cols_coords,
    )

from .tcommons import data_folder


BASE = '-db db.json -seq EGAAGAASS -fld fld.pdb'


def test_loop_closure_defaults():
    """Test ccd is the default method."""
    args = vars(cli_ldrs.ap.parse_args(BASE.split()))
    assert args['loop_closure'] == 'ccd'
    assert args['ccd_max_trials'] == 100
    assert cli_ldrs.LOOP_CLOSURE_METHODS[0] == 'ccd'


@pytest.mark.parametrize('method', cli_ldrs.LOOP_CLOSURE_METHODS)
def test_loop_closure_arguments(method):
    """Test both methods and the ccd budget parse."""
    cmd = f'{BASE} --loop-closure {method} --ccd-max-trials 7'
    args = vars(cli_ldrs.ap.parse_args(cmd.split()))
    assert args['loop_closure'] == method
    assert args['ccd_max_trials'] == 7


@pytest.mark.parametrize(
    'option',
    ['--loop-closure kic', '--ccd-solver dls'],
    )
def test_loop_closure_invalid_choice(option):
    """Test unknown methods and the removed solver option are refused."""
    with pytest.raises(SystemExit):
        cli_ldrs.ap.parse_args(f'{BASE} {option}'.split())


def test_main_refuses_unknown_methods():
    """Test main validates the method before doing any work."""
    with pytest.raises(ValueError):
        cli_ldrs.main('seq.fasta', 'db.json', folded_structure='f.pdb', loop_closure='kic')  # noqa: E501


def test_unchecked_pairs():
    """Test next_seeker checks every pair of conformers exactly once."""
    pairs = set()
    checked_C, checked_N = [], []
    C_files, N_files = [], []
    for round_ in range(3):
        C_files = C_files + [f'C{round_}{i}' for i in range(2)]
        N_files = N_files + [f'N{round_}{i}' for i in range(3)]
        groups = cli_ldrs._unchecked_pairs(
            C_files,
            N_files,
            checked_C,
            checked_N,
            )
        for tocheck_C, tocheck_N in groups:
            for pair in ((c, n) for c in tocheck_C for n in tocheck_N):
                assert pair not in pairs
                pairs.add(pair)
        checked_C, checked_N = C_files, N_files
        assert len(pairs) == len(C_files) * len(N_files)
    assert cli_ldrs._unchecked_pairs(C_files, N_files, C_files, N_files) == []


def test_half_seeds():
    """Test every round and side of next_seeker has seeds of its own."""
    assert cli_ldrs._half_seeds(7, 0, 0, 3) == [7, 8, 9]
    seeds = [
        seed
        for round_ in range(4)
        for side in (0, 1)
        for seed in cli_ldrs._half_seeds(7, round_, side, 3)
        ]
    assert len(set(seeds)) == len(seeds) == 24


def test_quiet_stdout(capfd):
    """Test writes to the stdout descriptor are dropped, then restored."""
    with cli_ldrs._quiet_stdout():
        os.write(1, b'hidden\n')
    os.write(1, b'shown\n')
    assert capfd.readouterr().out == 'shown\n'


def _helix(size):
    """N, CA, C, O of an ideal helix (NeRF, as in test_libccd)."""
    def place(a, b, c, bond, angle, torsion):
        bc = (c - b) / np.linalg.norm(c - b)
        n = np.cross(b - a, bc)
        n /= np.linalg.norm(n)
        m = np.cross(n, bc)
        return c + bond * (
            -math.cos(angle) * bc
            + math.sin(angle) * (math.cos(torsion) * m + math.sin(torsion) * n)  # noqa: E501
            )
    phi, psi, omega = (math.radians(x) for x in (-60.0, -45.0, 180.0))
    X = np.zeros((size, 4, 3))
    X[0, 1] = (1.458, 0, 0)
    angle = math.radians(111.2)
    X[0, 2] = X[0, 1] + 1.525 * np.array([-math.cos(angle), math.sin(angle), 0])  # noqa: E501
    for i in range(1, size):
        X[i, 0] = place(*X[i - 1, :3], 1.329, math.radians(116.2), psi)
        X[i, 1] = place(X[i - 1, 1], X[i - 1, 2], X[i, 0], 1.458, math.radians(121.7), omega)  # noqa: E501
        X[i, 2] = place(X[i - 1, 2], X[i, 0], X[i, 1], 1.525, math.radians(111.2), phi)  # noqa: E501
    for i in range(size):
        X[i, 3] = place(*X[i, :3], 1.231, math.radians(120.8), psi + math.pi)  # noqa: E501
    return X


def test_faspr_in_context():
    """Test FASPR packs a loop next to template residues it then drops."""
    pytest.importorskip('idpcpp')
    from scipy.spatial import cKDTree

    from idpconfgen.core.build_definitions import forcefields
    from idpconfgen.libs.libbuild import init_conflabels, init_confmasks
    loop_seq = 'AKLWEA'
    chain = _helix(12)
    loop = chain[:6]
    context_bb = chain[8:]  # a later part of the helix as neighbours
    packing = {
        'sequence': 'FYLR',
        'backbone': context_bb,
        'atoms': context_bb.reshape(-1, 3),
        'owner': np.repeat(np.arange(4), 4),
        }
    packing['tree'] = cKDTree(packing['atoms'])
    topology = forcefields['Amberff14SB'](add_OXT=True, add_Nterminal_H=True)
    labels = init_conflabels(loop_seq, topology.atom_names)
    masks = init_confmasks(labels.atom_labels)
    n_heavy = len(masks.non_Hs_non_OXT)
    packed = cli_ldrs._faspr_in_context(loop, loop_seq, packing, 2, 3, n_heavy)  # noqa: E501
    assert packed.shape == (n_heavy, 3)
    assert np.all(np.isfinite(packed))
    heavy_labels = labels.atom_labels[masks.non_Hs_non_OXT]
    heavy_res = labels.res_nums[masks.non_Hs_non_OXT]
    # the backbone given is kept and side chains hang on their residues
    for k, (name, res) in enumerate(zip(heavy_labels, heavy_res)):
        if name in ('N', 'CA', 'C', 'O'):
            atom = ('N', 'CA', 'C', 'O').index(name)
            assert np.allclose(packed[k], loop[res - 1, atom], atol=1e-3)
        elif name == 'CB':
            ca = loop[res - 1, 1]
            assert 1.4 < np.linalg.norm(packed[k] - ca) < 1.65
    # without neighbours in reach, the loop is packed alone
    far = dict(packing, tree=cKDTree(packing['atoms'] + 1000.0))
    alone = cli_ldrs._faspr_in_context(loop, loop_seq, far, 2, 3, n_heavy)
    assert alone.shape == (n_heavy, 3)


def _write_loop(path, offset):
    """Write an 8-residue C-alpha trace whose IDR residues move by `offset`."""
    lines = []
    for r in range(1, 9):
        x = 3.8 * r + (offset if 3 <= r <= 6 else 0.0)
        lines.append(
            f'ATOM  {r:5d}  CA  ALA A{r:4d}    '
            f'{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C'
            )
    path.write_text('\n'.join(lines) + '\n')


def _close_with(plan, tmp_path, monkeypatch, nconfs=3):
    """Run `_close_linker_ccd` with loops given by `plan(round)`."""
    rounds = []

    def build(task, *, output_folder, conformer_name, **kwargs):
        worker, quota, _seed, round_ = task
        offsets = plan(round_)[:quota]
        for n, offset in enumerate(offsets):
            name = f'{conformer_name}_{round_}_{worker}_{n}.pdb'
            _write_loop(output_folder / name, offset)
        rounds.append(round_)
        return {'built': len(offsets), 'trials': 1, 'reasons': {}}

    monkeypatch.setattr(cli_ldrs, '_build_ccd_loops', build)
    cli_ldrs._close_linker_ccd(
        'GGAAAAGG',
        index=0,
        ccd_context={},
        ccd_params=libccd.CCDParameters(),
        rama=None,
        output_folder=tmp_path,
        conformer_name='loop',
        nconfs=nconfs,
        ncores=1,
        random_seed=0,
        max_trials=1,
        max_clash=40,
        tolerance=0.4,
        report_folder=tmp_path,
        )
    return rounds, sorted(path.name for path in tmp_path.glob('*.pdb'))


def test_close_linker_ccd_keeps_closing(tmp_path, monkeypatch):
    """Test rounds go on after one closes nothing, until all are built."""
    plan = {0: [0.0], 1: [], 2: [5.0, 10.0]}
    rounds, files = _close_with(plan.get, tmp_path, monkeypatch)
    assert rounds == [0, 1, 2]
    assert len(files) == 3


def test_close_linker_ccd_repeats(tmp_path, monkeypatch):
    """Test repeated loops are deleted, and ldrs stops if only repeats come."""
    rounds, files = _close_with(lambda r: [0.0, 0.1, 0.2], tmp_path, monkeypatch)  # noqa: E501
    assert rounds == list(range(cli_ldrs.LOOP_STALL_ROUNDS + 1))
    assert files == ['loop_0_0_0.pdb']


def test_close_linker_ccd_stalls(tmp_path, monkeypatch):
    """Test rounds that close nothing new end the search with a warning."""
    rounds, files = _close_with(
        lambda r: [0.0] if r == 0 else [],
        tmp_path,
        monkeypatch,
        )
    assert rounds == list(range(cli_ldrs.LOOP_STALL_ROUNDS + 1))
    assert files == ['loop_0_0_0.pdb']


def test_close_linker_ccd_first_round(tmp_path, monkeypatch):
    """Test a first round without any closed loop stops ldrs."""
    with pytest.raises(IDPConfGenException):
        _close_with(lambda r: [], tmp_path, monkeypatch)


class _NoLog:
    def save(self, *args):
        pass


@pytest.mark.parametrize('case', ['N-IDR', 'C-IDR'])
def test_build_long_conformers(tmp_path, monkeypatch, case):
    """Test an IDR built in fragments comes out as if built in one piece."""
    structure = Structure(data_folder / 'EXPL_A.pdb')
    structure.build()
    rows = structure.data_array
    resseq = rows[:, col_resSeq].astype(int)
    sequence = 'AYIAKQRQISFVKSHF'  # EXPL_A, residues 4 to 19
    # IDR of 10 residues with the 2 folded residues it overlaps
    if case == 'N-IDR':
        folded = 14, 19
        idr = 4, 15
        anchor = ((14, 'C'), (15, 'N'), (15, 'CA'))
    else:
        folded = 4, 9
        idr = 8, 19
        anchor = ((8, 'C'), (9, 'N'), (9, 'CA'))
    fld = Structure(data_folder / 'EXPL_A.pdb')
    fld.build()
    fld._data_array = rows[(resseq >= folded[0]) & (resseq <= folded[1])]
    fld_xyz = np.array([
        rows[(resseq == r) & (rows[:, col_name] == n)][0, cols_coords]
        for r, n in anchor
        ]).astype(float)
    seq = sequence[idr[0] - 4:idr[1] - 3]
    parts = split_long_idr(seq, case, size=5)
    assert len(parts) == 2

    def fake_generator(*, index=None, input_seq=None, random_seed=0, **kw):
        first = 4 + sequence.index(input_seq)
        sel = (resseq >= first) & (resseq < first + len(input_seq))
        piece = rows[sel]
        numbers = resseq[sel] - first + 1
        yield piece[:, col_name], numbers, None
        rng = np.random.default_rng(random_seed)
        xyz = piece[:, cols_coords].astype(float)
        while True:
            rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            rotation *= np.sign(np.linalg.det(rotation))
            yield 0.0, xyz @ rotation.T + rng.normal(scale=10, size=3)

    monkeypatch.setattr(cli_ldrs, 'conformer_generator', fake_generator)
    monkeypatch.setattr(cli_ldrs, 'ENERGYLOGSAVER', _NoLog())
    monkeypatch.setattr(cli_ldrs, 'FRAGMENT_GLOBALS', {0: {}, 1: {}})
    cli_ldrs.RANDOMSEEDS.put(3)
    cli_ldrs.CONF_NUMBER.put(1)
    cli_ldrs._build_long_conformers(
        fragments=[(part, key) for key, part in enumerate(parts)],
        fld_xyz=fld_xyz,
        fld_struc=fld,
        disorder_case=case,
        max_clash=10**6,
        conformer_name='conformer',
        output_folder=tmp_path,
        )
    built = Structure(tmp_path / 'conformer_1.pdb')
    built.build()
    out = built.data_array
    native = rows[(resseq >= idr[0]) & (resseq <= idr[1])]
    assert list(out[:, col_name]) == list(native[:, col_name])
    assert sorted(set(out[:, col_resSeq].astype(int))) == list(range(1, 13))
    assert np.allclose(out[:, cols_coords].astype(float), native[:, cols_coords].astype(float), atol=0.05)  # noqa: E501
