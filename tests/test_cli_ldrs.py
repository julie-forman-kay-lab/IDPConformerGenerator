"""Test the Linker-IDR loop-closure options of cli_ldrs."""
import math
import os

import numpy as np
import pytest

from idpconfgen import cli_ldrs
from idpconfgen.libs import libccd


BASE = '-db db.json -seq EGAAGAASS -fld fld.pdb'


def test_loop_closure_defaults():
    """Test next_seeker stays the default method."""
    args = vars(cli_ldrs.ap.parse_args(BASE.split()))
    assert args['loop_closure'] == 'next_seeker'
    assert args['ccd_solver'] == 'dls'
    assert args['ccd_max_trials'] == 100
    assert cli_ldrs.LOOP_CLOSURE_METHODS[0] == 'next_seeker'


@pytest.mark.parametrize('solver', libccd.SOLVERS)
def test_loop_closure_ccd_arguments(solver):
    """Test the ccd method and its options parse."""
    cmd = f'{BASE} --loop-closure ccd --ccd-solver {solver} --ccd-max-trials 7'
    args = vars(cli_ldrs.ap.parse_args(cmd.split()))
    assert args['loop_closure'] == 'ccd'
    assert args['ccd_solver'] == solver
    assert args['ccd_max_trials'] == 7


@pytest.mark.parametrize(
    'option',
    ['--loop-closure kic', '--ccd-solver newton'],
    )
def test_loop_closure_invalid_choice(option):
    """Test unknown methods are refused by the parser."""
    with pytest.raises(SystemExit):
        cli_ldrs.ap.parse_args(f'{BASE} {option}'.split())


@pytest.mark.parametrize(
    'kwargs',
    [{'loop_closure': 'kic'}, {'ccd_solver': 'newton'}],
    )
def test_main_refuses_unknown_methods(kwargs):
    """Test main validates the methods before doing any work."""
    with pytest.raises(ValueError):
        cli_ldrs.main('seq.fasta', 'db.json', folded_structure='f.pdb', **kwargs)  # noqa: E501


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
