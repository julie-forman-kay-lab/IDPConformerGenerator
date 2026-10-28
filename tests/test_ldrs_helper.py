"""Test LDRS helpers."""
import os
from itertools import islice, product

import numpy as np
import pytest

from idpconfgen import ldrs_helper
from idpconfgen.core.build_definitions import build_bend_CA_C_O, distance_C_O
from idpconfgen.core.definitions import vdW_radii_tsai_1999
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.ldrs_helper import (
    _assemble_combinations,
    _random_combinations,
    align_coords,
    atom_vdw_radii,
    attachment_anchor,
    attachment_overlaps,
    backbone_anchor,
    build_in_fragments,
    combinations_clash_check,
    count_clashes,
    create_all_combinations,
    create_combinations,
    find_junctions,
    graft_fragment,
    half_backbone,
    inter_chain_cc,
    is_hydrogen,
    join_overlaps,
    junction_runs,
    linker_numbered,
    linker_reach,
    linker_stop,
    next_seeker,
    packing_context,
    prepare_ccd_context,
    prune_duplicate_loops,
    psurgeon,
    reachable_residues,
    splice_folded_rows,
    split_long_idr,
    unreachable_end,
    )
from idpconfgen.libs.libcalc import calc_torsion_angles
from idpconfgen.libs.libstructure import (
    Structure,
    col_chainID,
    col_name,
    col_resName,
    col_resSeq,
    col_segid,
    col_serial,
    cols_coords,
    structure_to_pdb,
    write_PDB,
    )

from .tcommons import data_folder


@pytest.fixture
def expl_rows():
    """Atoms of EXPL_A (residues 4 to 19)."""
    structure = Structure(data_folder / 'EXPL_A.pdb')
    structure.build()
    return structure.data_array


def test_create_all_combinations_linker_order(tmp_path, monkeypatch):
    """Test Linker-IDR folders are read in order, each with its matches."""
    linkers = tmp_path / 'chainA' / 'Linker-IDR'
    for num in (10, 2, 0, 1):
        match = linkers / f'linker_{num}' / f'{num}_match'
        match.mkdir(parents=True)
        (match / f'{num}_conf.pdb').write_text('')
    seen = []

    def fake_combinations(lst, num_combinations, ncores=1, **kwargs):
        seen.append(lst)
        return [[files[0] for files in lst]]

    monkeypatch.setattr(
        'idpconfgen.ldrs_helper.create_combinations',
        fake_combinations,
        )
    combos = create_all_combinations(tmp_path, ['A'], nconfs=1)
    names = [[path.name for path in files] for files in seen[0]]
    assert names == [['0_conf.pdb'], ['1_conf.pdb'], ['2_conf.pdb'], ['10_conf.pdb']]  # noqa: E501
    assert [path.name for path in combos['A'][0]] == [
        '0_conf.pdb', '1_conf.pdb', '2_conf.pdb', '10_conf.pdb',
        ]


def _write_backbone_fragment(rows, path, shift_n=()):
    """Write backbone rows renumbered from 1; optionally displace some N."""
    rows = np.array(rows, copy=True)
    resseq = rows[:, col_resSeq].astype(int)
    rows[:, col_resSeq] = (resseq - resseq.min() + 1).astype(str)
    for number in shift_n:
        atom = (rows[:, col_resSeq].astype(int) == number) & (rows[:, col_name] == 'N')  # noqa: E501
        rows[atom, cols_coords[0]] = str(float(rows[atom, cols_coords[0]][0]) + 5.0)  # noqa: E501
    write_PDB(structure_to_pdb(rows), path)


def test_next_seeker_writes_every_junction(expl_rows, tmp_path):
    """Test every junction of every pair is joined from the full halves."""
    resseq = expl_rows[:, col_resSeq].astype(int)
    backbone = expl_rows[
        (resseq >= 4) & (resseq <= 11)
        & np.isin(expl_rows[:, col_name], ('N', 'CA', 'C', 'O'))
        ]
    lower = tmp_path / 'lower.pdb'
    upper_a = tmp_path / 'upper_a.pdb'
    upper_b = tmp_path / 'upper_b.pdb'
    _write_backbone_fragment(backbone, lower)
    _write_backbone_fragment(backbone, upper_a)
    # no junction before residues 2 to 4 in the second N-side conformer
    _write_backbone_fragment(backbone, upper_b, shift_n=(2, 3, 4))
    matches = tmp_path / 'matches'
    matches.mkdir()
    found = next_seeker(lower, [upper_a, upper_b], 40, 0.4, matches)
    expected = (
        [f'lower+upper_a_j{r}.pdb' for r in range(1, 8)]
        + [f'lower+upper_b_j{r}.pdb' for r in range(4, 8)]
        )
    assert found == len(expected)
    assert sorted(path.name for path in matches.iterdir()) == sorted(expected)
    for name in expected:
        joined = Structure(matches / name)
        joined.build()
        rows = joined.data_array
        numbers = rows[:, col_resSeq].astype(int)
        assert sorted(set(numbers)) == list(range(1, 9))
        assert len(rows) == len(backbone)
        # the carbonyl O of the junction lies in the new peptide plane
        term = int(name[:-4].rsplit('_j', 1)[1])

        def xyz(residue, atom):
            row = (numbers == residue) & (rows[:, col_name] == atom)  # noqa: B023, E501
            return rows[row][0, cols_coords].astype(float)  # noqa: B023

        c, o = xyz(term, 'C'), xyz(term, 'O')
        ca, n = xyz(term, 'CA'), xyz(term + 1, 'N')
        assert abs(np.linalg.norm(o - c) - distance_C_O) < 2e-3
        cos = np.dot(ca - c, o - c) / np.linalg.norm(ca - c) / np.linalg.norm(o - c)  # noqa: E501
        assert abs(np.arccos(cos) - 2 * build_bend_CA_C_O) < 2e-3
        normal = np.cross(ca - c, n - c)
        assert abs(np.dot(normal / np.linalg.norm(normal), o - c)) < 5e-3


def test_next_seeker_trimmed_halves(expl_rows, tmp_path):
    """Test halves without some residues are joined by residue number."""
    resseq = expl_rows[:, col_resSeq].astype(int)
    backbone = expl_rows[
        (resseq >= 4) & (resseq <= 11)
        & np.isin(expl_rows[:, col_name], ('N', 'CA', 'C', 'O'))
        ]
    lower = tmp_path / 'lower.pdb'
    upper = tmp_path / 'upper.pdb'
    _write_backbone_fragment(backbone[backbone[:, col_resSeq].astype(int) <= 8], lower)  # noqa: E501
    rows = np.array(backbone, copy=True)
    rows[:, col_resSeq] = (rows[:, col_resSeq].astype(int) - 3).astype(str)
    write_PDB(structure_to_pdb(rows[rows[:, col_resSeq].astype(int) >= 4]), upper)  # noqa: E501
    matches = tmp_path / 'matches'
    matches.mkdir()
    # C-side residues 1 to 5, N-side residues 4 to 8: junctions at 3 to 5
    assert next_seeker(lower, [upper], 40, 0.4, matches) == 3
    for term in (3, 4, 5):
        joined = Structure(matches / f'lower+upper_j{term}.pdb')
        joined.build()
        numbers = joined.data_array[:, col_resSeq].astype(int)
        assert sorted(set(numbers)) == list(range(1, 9))


def _backbone_array(expl_rows, first, last):
    rows = _chain_rows(expl_rows, first, last)
    return half_backbone(rows, last - first + 1)


def test_find_junctions(expl_rows):
    """Test junctions are found as the scalar next_seeker check finds them."""
    native = _backbone_array(expl_rows, 4, 13)
    rng = np.random.default_rng(5)
    c_side = native[None] + rng.normal(0, 0.06, size=(3,) + native.shape)
    n_side = native[None] + rng.normal(0, 0.06, size=(4,) + native.shape)
    n_side[1, 6:] += 3.0  # no junction from residue 6 on
    c_side[2, 4:] = np.nan  # a C-side half cut after residue 4
    found = find_junctions(c_side, n_side)
    expected = []
    for c, n, residue in product(range(3), range(4), range(1, 10)):
        ca, cc = c_side[c, residue - 1, 1], c_side[c, residue - 1, 2]
        nn, ca2 = n_side[n, residue, 0], n_side[n, residue, 1]
        if not np.all(np.isfinite([ca, cc, nn, ca2])):
            continue
        omega = calc_torsion_angles(np.array([ca, cc, nn, ca2]))
        if (
                1.32 <= ldrs_helper.calculate_distance(cc, nn) <= 1.56
                and 1.91 <= ldrs_helper.calculate_angle(ca, cc, nn) <= 2.15
                and 2.2 <= ldrs_helper.calculate_distance(cc, ca2) <= 2.7
                and np.abs(omega)[0] >= 2.61
                ):
            expected.append((c, n, residue))
    assert len(expected) > 20
    assert [tuple(row) for row in found.tolist()] == expected
    assert find_junctions(c_side[:0], n_side).shape == (0, 3)


def test_junction_runs():
    """Test junctions are grouped by pair in runs of consecutive residues."""
    found = np.array([
        (0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 1, 7), (0, 2, 5), (3, 1, 6),
        ])
    assert junction_runs(found) == {
        (0, 1): [[2, 3, 4], [7]],
        (0, 2): [[5]],
        (3, 1): [[6]],
        }


def _straight_half(length, direction=1.0):
    """Rows of an extended chain along x, residues numbered 1 to `length`."""
    rows = []
    for residue in range(1, length + 1):
        x = direction * 3.6 * (residue - 1)
        for k, (name, element) in enumerate((('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O'))):  # noqa: E501
            rows.append([
                'ATOM', str(4 * residue + k), name, '', 'ALA', 'A',
                str(residue), '', f'{x + 0.5 * k:.3f}', '0.000', '0.000',
                '1.00', '0.00', '', element, '',
                ])
    return np.array(rows)


@pytest.mark.parametrize('case', ['C-IDR', 'N-IDR'])
def test_reachable_residues(case):
    """Test a half is cut where its residues can no longer reach the far end."""
    length = 12
    rows = _straight_half(length)
    xyz = rows[:, cols_coords].astype(float)
    numbers = rows[:, col_resSeq].astype(int)
    # the far end next to residue 1 (C-side) or to residue 12 (N-side)
    far = xyz[numbers == (1 if case == 'C-IDR' else length)].mean(axis=0) + (0, 3.0, 0)  # noqa: E501
    atom = 'C' if case == 'C-IDR' else 'N'
    pick = rows[:, col_name] == atom
    distance = np.linalg.norm(xyz[pick] - far, axis=1)
    reach = distance <= linker_reach(case, numbers[pick], length)
    kept = reachable_residues(rows, case, far, length)
    kept_numbers = sorted(set(kept[:, col_resSeq].astype(int)))
    if case == 'C-IDR':
        assert not reach.all()
        assert kept_numbers == list(range(1, int(np.argmin(reach)) + 1))
    else:
        first = numbers[pick][~reach].max() + 1
        assert kept_numbers == list(range(first, length + 1))
    # the bound is monotonic: past the first residue out of reach, none
    # is within it
    assert np.all(np.diff(linker_reach(case, np.arange(1, length + 1), length)) * (1 if case == 'N-IDR' else -1) > 0)  # noqa: E501


def test_reachable_residues_drops_incomplete_end():
    """Test a residue without its C at the free end of a C-side half goes."""
    rows = _straight_half(6)
    numbers = rows[:, col_resSeq].astype(int)
    rows = rows[~((numbers == 6) & np.isin(rows[:, col_name], ('CA', 'C', 'O')))]  # noqa: E501
    kept = reachable_residues(rows, 'C-IDR', np.zeros(3), 12)
    assert sorted(set(kept[:, col_resSeq].astype(int))) == [1, 2, 3, 4, 5]


def test_linker_stop(expl_rows):
    """Test the builder stop matches the reach of the aligned residues."""
    length = 16
    rows = _chain_rows(expl_rows, 4, 15)
    rows = rows[np.isin(rows[:, col_name], ('N', 'CA', 'C', 'O'))]
    target = np.array([[1.0, 2.0, 3.0], [2.2, 2.6, 3.0], [3.4, 2.0, 3.6]])
    aligned = align_coords(np.array(rows, copy=True), target, 'C-IDR')
    c_atoms = aligned[aligned[:, col_name] == 'C'][:, cols_coords].astype(float)  # noqa: E501
    limit = linker_reach('C-IDR', np.arange(1, 13), length)
    far = c_atoms[7] + (0.0, 0.0, limit[7] + 1.0)
    out_of_reach = np.linalg.norm(c_atoms - far, axis=1) > limit
    assert out_of_reach.any() and not out_of_reach.all()
    # the builder backbone: N, CA, C of each residue in another frame
    xyz = rows[np.isin(rows[:, col_name], ('N', 'CA', 'C'))][:, cols_coords].astype(float)  # noqa: E501
    rotation, _ = np.linalg.qr(np.random.default_rng(3).normal(size=(3, 3)))
    rotation *= np.sign(np.linalg.det(rotation))
    xyz = xyz @ rotation.T + (5.0, -7.0, 2.0)
    stop = linker_stop(target, far, length)
    for built in range(1, 13):
        backbone = np.full_like(xyz, np.nan)
        backbone[:3 * built] = xyz[:3 * built]
        expected = built >= 2 and bool(out_of_reach[built - 1])
        assert stop(backbone) == expected


def test_unreachable_end_and_numbering():
    """Test the free end of halves built from their folded end."""
    length = 12
    rows = _straight_half(5)
    xyz = rows[:, cols_coords].astype(float)
    # an N-side half of 5 residues is the end of the Linker-IDR
    renumbered = linker_numbered(rows, 'N-IDR', length)
    assert sorted(set(renumbered[:, col_resSeq].astype(int))) == list(range(8, 13))  # noqa: E501
    assert linker_numbered(rows, 'C-IDR', length) is rows
    n_first = xyz[0]
    limit = linker_reach('N-IDR', 8, length)
    assert not unreachable_end(rows, 'N-IDR', n_first + (limit - 0.1, 0, 0), length)  # noqa: E501
    assert unreachable_end(rows, 'N-IDR', n_first + (limit + 0.1, 0, 0), length)  # noqa: E501
    c_last = xyz[(rows[:, col_resSeq] == '5') & (rows[:, col_name] == 'C')][0]
    limit = linker_reach('C-IDR', 5, length)
    assert not unreachable_end(rows, 'C-IDR', c_last + (0, limit - 0.1, 0), length)  # noqa: E501
    assert unreachable_end(rows, 'C-IDR', c_last + (0, limit + 0.1, 0), length)  # noqa: E501


def test_atom_vdw_radii_and_hydrogens():
    """Test radii follow count_clashes and fall back to atom names."""
    radii = atom_vdw_radii(['C', '', 'XX'], ['CA', '1HB', 'ZN9'], default=1.5)
    assert radii[0] == vdW_radii_tsai_1999['C']
    assert radii[1] == vdW_radii_tsai_1999['H']
    assert radii[2] == 1.5
    assert list(is_hydrogen(['H', '', 'C', ''], ['H', 'HA', 'CB', 'N'])) == [
        True, True, False, False,
        ]


def _anchor_xyz(rows, atoms):
    """Coordinates of the (residue, atom name) pairs of `rows`."""
    resseq = rows[:, col_resSeq].astype(int)
    return np.array([
        rows[(resseq == res) & (rows[:, col_name] == name)][0, cols_coords]
        for res, name in atoms
        ]).astype(float)


def _handedness(rows):
    """Sign of (N - CA) x (C - CA) . (CB - CA) of the residues with a CB."""
    resseq = rows[:, col_resSeq].astype(int)
    signs = []
    for res in sorted(set(resseq)):
        names = rows[resseq == res][:, col_name]
        if 'CB' not in names:
            continue
        n, ca, c, cb = _anchor_xyz(rows, [(res, a) for a in ('N', 'CA', 'C', 'CB')])  # noqa: E501
        signs.append(np.sign(np.dot(np.cross(n - ca, c - ca), cb - ca)))
    return np.array(signs)


def _fragment_and_anchors(expl_rows, case):
    rows = np.array(expl_rows, copy=True)
    resseq = rows[:, col_resSeq].astype(int)
    rows[:, col_resSeq] = (resseq - resseq.min() + 1).astype(str)
    last = int(rows[:, col_resSeq].astype(int).max())
    if case == 'N-IDR':
        anchors = [(last - 1, 'C'), (last, 'N'), (last, 'CA')]
    else:
        anchors = [(1, 'C'), (2, 'N'), (2, 'CA')]
    return rows, anchors


@pytest.mark.parametrize('case', ['N-IDR', 'C-IDR'])
def test_align_coords(expl_rows, case):
    rows, anchors = _fragment_and_anchors(expl_rows, case)
    original = rows[:, cols_coords].astype(float)
    target = _anchor_xyz(rows, anchors)
    handedness = _handedness(rows)
    assert len(handedness) and np.all(handedness == handedness[0])

    for seed in range(100):
        # move the fragment away by a random rotation and translation
        rng = np.random.default_rng(seed)
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        rotation *= np.sign(np.linalg.det(rotation))
        moved = original @ rotation.T + rng.normal(scale=20, size=3)
        rows[:, cols_coords] = np.round(moved, 3).astype(str)

        aligned = align_coords(rows, target, case)
        assert np.allclose(_anchor_xyz(aligned, anchors), target, atol=0.005)
        assert np.allclose(aligned[:, cols_coords].astype(float), original, atol=0.05)  # noqa: E501
        assert np.array_equal(_handedness(aligned), handedness)


def test_backbone_anchor(expl_rows):
    """Test N, CA, C, O of one residue, NaN when missing."""
    rows = expl_rows[expl_rows[:, col_resSeq].astype(int) == 5]
    anchor = backbone_anchor(rows)
    names = list(rows[:, col_name])
    assert np.allclose(anchor[1], rows[names.index('CA'), cols_coords].astype(float))  # noqa: E501
    no_oxygen = rows[rows[:, col_name] != 'O']
    assert np.all(np.isnan(backbone_anchor(no_oxygen)[3]))


def test_prepare_ccd_context(expl_rows):
    """Test anchors and environment masks for a gap between 8 and 12."""
    resseq = expl_rows[:, col_resSeq].astype(int)
    chain = expl_rows[(resseq <= 8) | (resseq >= 12)]
    context = prepare_ccd_context(chain, chain, 8, 12)
    assert context['lower_anchor'].shape == (4, 3)
    assert np.all(np.isfinite(context['lower_anchor']))
    assert len(context['env_coords']) == len(chain)
    lower = chain[:, col_resSeq].astype(int) == 8
    upper = chain[:, col_resSeq].astype(int) == 12
    assert np.array_equal(context['env_lower'], lower)
    assert np.array_equal(context['env_upper'], upper)
    assert len(context['lower_rows']) == lower.sum()
    assert not context['env_hydrogen'].any()
    outer = chain[chain[:, col_resSeq].astype(int) == 7]
    assert np.allclose(context['lower_outer'], backbone_anchor(outer))
    assert np.all(np.isfinite(context['upper_outer']))
    # the four overlap residues 7, 8, 12 and 13 are not packing context
    packing = context['packing']
    assert len(packing['sequence']) == len(set(resseq)) - 3 - 4
    assert packing['backbone'].shape == (len(packing['sequence']), 4, 3)
    with pytest.raises(IDPConfGenException):
        prepare_ccd_context(chain, chain, 9, 12)
    # no residue before or after the anchors
    edge = prepare_ccd_context(chain, chain, 4, 12)
    assert np.isnan(edge['lower_outer']).all()


def test_packing_context(expl_rows):
    """Test which template residues side chains are packed next to."""
    rows = np.array(expl_rows, copy=True)
    resseq = rows[:, col_resSeq].astype(int)
    rows[resseq == 5, col_resName] = 'HIP'
    rows[resseq == 6, col_resName] = 'HOH'
    rows = rows[~((resseq == 9) & (rows[:, col_name] == 'O'))]
    resseq = rows[:, col_resSeq].astype(int)
    exclude = resseq == 10
    context = packing_context(rows, exclude)
    kept = [r for r in sorted(set(resseq)) if r not in (6, 9, 10)]
    assert len(context['sequence']) == len(kept)
    assert context['sequence'][1] == 'H'  # HIP of residue 5
    heavy = sum(int(np.sum(resseq == r)) for r in kept)
    assert len(context['atoms']) == len(context['owner']) == heavy
    assert np.array_equal(np.unique(context['owner']), np.arange(len(kept)))
    empty = packing_context(rows[:0])
    assert empty['sequence'] == '' and empty['atoms'].shape == (0, 3)


def test_splice_folded_rows(expl_rows):
    """Test junction residues are taken from the template."""
    resseq = expl_rows[:, col_resSeq].astype(int)
    # a fake loop fragment numbered 1..8 made of residues 4..11
    fragment = np.array(expl_rows[(resseq >= 4) & (resseq <= 11)], copy=True)
    fragment[:, col_resSeq] = (fragment[:, col_resSeq].astype(int) - 3).astype(str)  # noqa: E501
    fragment[:, col_chainID] = 'B'
    lower = expl_rows[resseq == 15]
    upper = expl_rows[resseq == 16]
    out = splice_folded_rows(fragment, lower, upper)
    out_res = out[:, col_resSeq].astype(int)
    assert sorted(set(out_res)) == list(range(1, 9))
    assert np.array_equal(out[out_res == 2][:, cols_coords], lower[:, cols_coords])  # noqa: E501
    assert np.array_equal(out[out_res == 7][:, cols_coords], upper[:, cols_coords])  # noqa: E501
    assert set(out[:, col_chainID]) == {'B'}
    assert list(out[:, col_serial]) == [str(i) for i in range(1, len(out) + 1)]


def _write_residues(rows, path, first, last, shift=0.0, renumber=True):
    """Write residues `first` to `last`, numbered from 1, moved in x."""
    resseq = rows[:, col_resSeq].astype(int)
    out = np.array(rows[(resseq >= first) & (resseq <= last)], copy=True)
    if renumber:
        out[:, col_resSeq] = (out[:, col_resSeq].astype(int) - first + 1).astype(str)  # noqa: E501
    x = out[:, cols_coords[0]].astype(float) + shift
    out[:, cols_coords[0]] = [f'{value:.3f}' for value in x]
    write_PDB(structure_to_pdb(out), path)
    return path


@pytest.fixture
def fragments(expl_rows, tmp_path):
    """IDR fragments 100 A apart; fragment 'X' sits on fragment 'a'."""
    paths = {}
    for i, name in enumerate('abcd'):
        paths[name] = _write_residues(
            expl_rows, tmp_path / f'{name}.pdb', 4, 11, shift=100.0 * i,
            )
    paths['X'] = _write_residues(expl_rows, tmp_path / 'X.pdb', 4, 11)
    return paths


@pytest.mark.parametrize(
    'combo,clash',
    [
        (('a', 'b'), False),
        (('X', 'a'), True),  # N-IDR + C-IDR
        (('X', ['b', 'a']), True),  # N-IDR + Linker-IDRs
        ((['b', 'X'], 'a'), True),  # Linker-IDRs + C-IDR
        (('b', ['X'], 'a'), True),  # N-IDR + one Linker-IDR + C-IDR
        (('X', ['b', 'c'], 'a'), True),  # N-IDR against C-IDR
        (('X', 'b', 'a'), True),  # three Linker-IDRs
        (('a', ['b', 'c'], 'd'), False),
        ],
    )
def test_combinations_clash_check(fragments, combo, clash):
    """Test every pair of IDRs of a combination is clash-checked."""
    selected = tuple(
        [fragments[name] for name in item] if isinstance(item, list)
        else fragments[item]
        for item in combo
        )
    result = combinations_clash_check(selected)
    if clash:
        assert result is False
    else:
        assert result == selected


def test_combinations_clash_check_carried_residues(expl_rows, tmp_path):
    """Test folded residues carried by close IDRs are not compared."""
    # N-IDR 4-8 carries 9-10, the Linker-IDR 11-14 carries 9-10, 15-16
    nidr = _write_residues(expl_rows, tmp_path / 'n.pdb', 4, 10)
    lidr = _write_residues(expl_rows, tmp_path / 'l.pdb', 9, 16)
    selected = (nidr, [lidr])
    assert combinations_clash_check(selected) is False
    cases = ['N-IDR', 'Linker-IDR']
    assert combinations_clash_check(selected, cases) == selected


def test_random_combinations():
    """Test combinations come once each, also from large products."""
    lists = [[1, 2, 3], ['a', 'b'], [None]]
    combos = list(_random_combinations(lists))
    assert sorted(combos, key=str) == sorted(product(*lists), key=str)
    first = list(islice(_random_combinations([list(range(1000))] * 3), 50))
    assert len(set(first)) == 50
    assert all(len(combo) == 3 for combo in first)


def test_create_combinations(fragments):
    """Test only combinations without clashes are returned."""
    nidr = [fragments['X'], fragments['b']]
    combos = create_combinations([nidr, [fragments['a']]], 2)
    assert combos == [[fragments['b'], fragments['a']]]
    linkers = [[fragments['c'], fragments['a']]]
    combos = create_combinations([nidr, linkers], 5)
    assert combos == [[fragments['b'], [fragments['c'], fragments['a']]]]
    with pytest.raises(IDPConfGenException):
        create_combinations([[fragments['X']], [fragments['a']]], 2)


def test_assemble_combinations_one_to_one():
    """Test every item is used once when all combinations pass."""
    chosen, tried = _assemble_combinations([6, 6, 6], 6, lambda idx: True, 100)
    assert len(chosen) == 6 and tried == 6
    for position in range(3):
        assert sorted(idx[position] for idx in chosen) == list(range(6))


def test_assemble_combinations_short_list():
    """Test a short list is used evenly and the others once each."""
    chosen, tried = _assemble_combinations([2, 6, 6], 6, lambda idx: True, 100)  # noqa: E501
    assert len(chosen) == len(set(chosen)) == 6 and tried == 6
    first = [idx[0] for idx in chosen]
    assert first.count(0) == first.count(1) == 3
    for position in (1, 2):
        assert sorted(idx[position] for idx in chosen) == list(range(6))


def test_assemble_combinations_swaps(monkeypatch):
    """Test items left over are placed by an exchange, not reused."""
    # without shuffling, (0, 0) is tried first and fails
    monkeypatch.setattr(
        'idpconfgen.ldrs_helper.random.shuffle',
        lambda items: None,
        )
    chosen, _ = _assemble_combinations([3, 3], 3, lambda idx: idx != (0, 0), 100)  # noqa: E501
    assert len(chosen) == 3 and (0, 0) not in chosen
    for position in range(2):
        assert sorted(idx[position] for idx in chosen) == [0, 1, 2]


def test_assemble_combinations_reuse():
    """Test rejected items are paired again and reuse is a last resort."""
    # item 0 of the first list never passes
    chosen, _ = _assemble_combinations([4, 4], 4, lambda idx: idx[0] != 0, 100)  # noqa: E501
    assert len(chosen) == len(set(chosen)) == 4
    assert all(idx[0] != 0 for idx in chosen)
    # the second list is used once each; one item of the first is reused
    assert sorted(idx[1] for idx in chosen) == [0, 1, 2, 3]
    first = [idx[0] for idx in chosen]
    assert max(first.count(i) for i in (1, 2, 3)) == 2
    # all four distinct combinations of two lists of two
    chosen, tried = _assemble_combinations([2, 2], 10, lambda idx: True, 100)
    assert sorted(chosen) == [(0, 0), (0, 1), (1, 0), (1, 1)] and tried == 4
    chosen, tried = _assemble_combinations([2, 2], 10, lambda idx: False, 100)
    assert chosen == [] and tried == 4
    chosen, tried = _assemble_combinations([5, 5], 10, lambda idx: False, 3)
    assert chosen == [] and tried == 3


def test_create_all_combinations_uses_each_conformer_once(expl_rows, tmp_path):  # noqa: E501
    """Test a chain with N-IDR, two Linker-IDRs and C-IDR."""
    chain = tmp_path / 'chainA'
    folders = {
        'N': chain / 'N-IDR',
        'L0': chain / 'Linker-IDR' / 'linker_0' / '0_match',
        'L1': chain / 'Linker-IDR' / 'linker_1' / '1_match',
        'C': chain / 'C-IDR',
        }
    files = {}
    for k, (name, folder) in enumerate(folders.items()):
        folder.mkdir(parents=True)
        files[name] = [
            _write_residues(
                expl_rows,
                folder / f'{name}_{i}.pdb',
                4,
                11,
                shift=100.0 * (3 * k + i),
                )
            for i in range(3)
            ]
    combos = create_all_combinations(tmp_path, ['A'], nconfs=3)['A']
    assert len(combos) == 3
    for combo in combos:
        assert len(combo) == 3 and len(combo[1]) == 2
    used = {
        'N': [combo[0] for combo in combos],
        'L0': [combo[1][0] for combo in combos],
        'L1': [combo[1][1] for combo in combos],
        'C': [combo[2] for combo in combos],
        }
    for name, paths in used.items():
        assert sorted(p.name for p in paths) == sorted(f.name for f in files[name])  # noqa: E501


def test_inter_chain_cc(fragments):
    """Test chains with nested Linker-IDRs and fewer entries than nconfs."""
    nested = [fragments['a'], [fragments['b']]]
    all_files = {'A': [nested], 'B': [fragments['c'], fragments['X']]}
    combos = inter_chain_cc(all_files, 3)
    assert combos == [{'A': nested, 'B': [fragments['c']]}]
    with pytest.raises(IDPConfGenException):
        inter_chain_cc({'A': [fragments['a']], 'B': [fragments['X']]}, 2)


def _template(rows, missing):
    """Return a built structure of `rows` without residues `missing`."""
    resseq = rows[:, col_resSeq].astype(int)
    template = Structure('\n'.join(structure_to_pdb(
        rows[~np.isin(resseq, missing)]
        )))
    template.build()
    return template


def _assert_same_residues(out, rows):
    assert list(out[:, col_resSeq].astype(int)) == list(rows[:, col_resSeq].astype(int))  # noqa: E501
    assert list(out[:, col_name]) == list(rows[:, col_name])
    assert list(out[:, col_resName]) == list(rows[:, col_resName])
    assert np.allclose(
        out[:, cols_coords].astype(float),
        rows[:, cols_coords].astype(float),
        )


@pytest.mark.parametrize('nested', [False, True])
def test_psurgeon_linkers_only(expl_rows, tmp_path, nested):
    """Test each Linker-IDR goes into its own chain break."""
    # sequence position 0 is residue 4; Linker-IDRs at 9-10 and 15-16
    linkers = [
        _write_residues(expl_rows, tmp_path / 'l0.pdb', 7, 12),
        _write_residues(expl_rows, tmp_path / 'l1.pdb', 13, 18),
        ]
    out = psurgeon(
        {'A': [linkers] if nested else linkers},
        _template(expl_rows, [9, 10, 15, 16]),
        {'A': ['Linker-IDR', 'Linker-IDR']},
        {'A': [(5, 7), (11, 13)]},
        )
    _assert_same_residues(out, expl_rows)


def test_psurgeon_n_linker_c(expl_rows, tmp_path):
    """Test N-IDR, Linker-IDR and C-IDR keep the template numbering."""
    # sequence position 0 is residue 4; the template starts at residue 7
    nidr = _write_residues(expl_rows, tmp_path / 'n.pdb', 4, 8)
    lidr = _write_residues(expl_rows, tmp_path / 'l.pdb', 10, 15)
    cidr = _write_residues(expl_rows, tmp_path / 'c.pdb', 16, 19)
    out = psurgeon(
        {'A': [nidr, [lidr], cidr]},
        _template(expl_rows, [4, 5, 6, 12, 13, 18, 19]),
        {'A': ['N-IDR', 'Linker-IDR', 'C-IDR']},
        {'A': [(0, 3), (8, 10), (14, 16)]},
        )
    _assert_same_residues(out, expl_rows)


def test_prune_duplicate_loops(expl_rows, tmp_path):
    """Test loops repeating another loop of the folder are deleted."""
    loops = tmp_path / 'matches'
    loops.mkdir()
    # 8-residue fragments: residues 3 to 6 are the IDR
    _write_residues(expl_rows, loops / 'a.pdb', 4, 11)
    _write_residues(expl_rows, loops / 'b.pdb', 4, 11, shift=0.3)
    _write_residues(expl_rows, loops / 'c.pdb', 4, 11, shift=3.0)
    # different folded residues only: still a repeat of 'a'
    resseq = expl_rows[:, col_resSeq].astype(int)
    rows = np.array(expl_rows, copy=True)
    moved = (resseq == 4) | (resseq == 11)
    x = rows[moved, cols_coords[0]].astype(float) + 5.0
    rows[moved, cols_coords[0]] = [f'{value:.3f}' for value in x]
    _write_residues(rows, loops / 'd.pdb', 4, 11)
    assert prune_duplicate_loops(loops) == 2
    assert sorted(p.name for p in loops.iterdir()) == ['a.pdb', 'c.pdb']
    assert prune_duplicate_loops(loops, min_rmsd=5.0) == 1


def test_prune_duplicate_loops_keeps_earlier(expl_rows, tmp_path):
    """Test a new loop is deleted, not the earlier one it repeats."""
    loops = tmp_path / 'matches'
    loops.mkdir()
    old = _write_residues(expl_rows, loops / 'z.pdb', 4, 11)
    new = _write_residues(expl_rows, loops / 'a.pdb', 4, 11, shift=0.3)
    os.utime(old, (1000, 1000))
    os.utime(new, (2000, 2000))
    assert prune_duplicate_loops(loops) == 1
    assert sorted(p.name for p in loops.iterdir()) == ['z.pdb']


@pytest.mark.parametrize('case', ['N-IDR', 'C-IDR'])
def test_split_long_idr(case):
    """Test long IDRs split into nearly equal, overlapping fragments."""
    folded = 'QK'
    idr = ''.join('ACDEFGHIKLMNPQRSTVWY'[i % 20] for i in range(310))
    seq = idr + folded if case == 'N-IDR' else folded + idr
    assert split_long_idr(seq[:152] if case == 'C-IDR' else seq[-152:], case) == [seq[:152] if case == 'C-IDR' else seq[-152:]]  # noqa: E501
    fragments = split_long_idr(seq, case)
    assert len(fragments) == 3
    assert [len(f) - 2 for f in fragments] == [104, 103, 103]
    # the first fragment built holds the folded residues; each next one
    # shares two residues with the one before
    if case == 'C-IDR':
        assert fragments[0].startswith(folded)
        for before, after in zip(fragments, fragments[1:]):
            assert after[:2] == before[-2:]
        joined = fragments[0] + ''.join(f[2:] for f in fragments[1:])
    else:
        assert fragments[0].endswith(folded)
        for before, after in zip(fragments, fragments[1:]):
            assert after[-2:] == before[:2]
        joined = ''.join(f[:-2] for f in reversed(fragments[1:])) + fragments[0]  # noqa: E501
    assert joined == seq


@pytest.mark.parametrize('case', ['N-IDR', 'C-IDR'])
def test_split_long_idr_linker(case):
    """Test the folded residues past a Linker-IDR are not IDR residues."""
    idr = ''.join('ACDEFGHIKLMNPQRSTVWY'[i % 20] for i in range(151))
    seq = 'QK' + idr[:150] + 'RL'
    assert split_long_idr(seq, case, extra=2) == [seq]
    seq = 'QK' + idr + 'RL'
    fragments = split_long_idr(seq, case, extra=2)
    # 76 and 75 IDR residues; the last fragment built also holds the two
    # folded residues at the other end
    assert [len(f) - 2 for f in fragments] == [76, 77]
    if case == 'C-IDR':
        assert fragments[0].startswith('QK') and fragments[1].endswith('RL')
        assert fragments[0] + fragments[1][2:] == seq
    else:
        assert fragments[0].endswith('RL') and fragments[1].startswith('QK')
        assert fragments[1][:-2] + fragments[0] == seq


def _chain_rows(expl_rows, first, last):
    """Residues first..last of EXPL_A renumbered from 1."""
    resseq = expl_rows[:, col_resSeq].astype(int)
    rows = np.array(expl_rows[(resseq >= first) & (resseq <= last)], copy=True)
    numbers = rows[:, col_resSeq].astype(int)
    rows[:, col_resSeq] = (numbers - first + 1).astype(str)
    return rows


def _moved(rows, seed):
    rows = np.array(rows, copy=True)
    rng = np.random.default_rng(seed)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rotation *= np.sign(np.linalg.det(rotation))
    xyz = rows[:, cols_coords].astype(float) @ rotation.T + rng.normal(scale=15, size=3)  # noqa: E501
    rows[:, cols_coords] = np.round(xyz, 3).astype(str)
    return rows


@pytest.mark.parametrize('case', ['N-IDR', 'C-IDR'])
def test_graft_fragment_rebuilds_the_chain(expl_rows, case):
    """Test a chain cut in two overlapping pieces is put back together."""
    first = int(expl_rows[:, col_resSeq].astype(int).min())
    last = int(expl_rows[:, col_resSeq].astype(int).max())
    whole = _chain_rows(expl_rows, first, last)
    if case == 'C-IDR':
        built, piece = (first, first + 8), (first + 7, last)
    else:
        built, piece = (first + 7, last), (first, first + 8)
    chain = _chain_rows(expl_rows, *built)
    fragment = _moved(_chain_rows(expl_rows, *piece), 1)
    target = attachment_anchor(chain, case)
    aligned = align_coords(fragment, target, case)
    assert attachment_overlaps(aligned, chain, case, target) == 0
    joined = graft_fragment(chain, aligned, case)
    assert len(joined) == len(whole)
    assert list(joined[:, col_name]) == list(whole[:, col_name])
    assert list(joined[:, col_resSeq]) == list(whole[:, col_resSeq])
    assert np.allclose(joined[:, cols_coords].astype(float), whole[:, cols_coords].astype(float), atol=0.05)  # noqa: E501


@pytest.mark.parametrize('case', ['N-IDR', 'C-IDR'])
def test_join_overlaps(expl_rows, case):
    """Test atoms across a join are checked, the bonded residues not."""
    first = int(expl_rows[:, col_resSeq].astype(int).min())
    last = int(expl_rows[:, col_resSeq].astype(int).max())
    # the chain has 9 residues, 8 of them kept; the join is between the
    # residues 8 and 9 of the result
    if case == 'C-IDR':
        built, piece, old, new = (first, first + 8), (first + 7, last), 8, 10  # noqa: E501
    else:
        built, piece, old, new = (first + 7, last), (first, first + 8), 9, 7  # noqa: E501
    chain = _chain_rows(expl_rows, *built)
    target = attachment_anchor(chain, case)
    aligned = align_coords(_moved(_chain_rows(expl_rows, *piece), 3), target, case)  # noqa: E501
    joined = graft_fragment(chain, aligned, case)
    # the peptide bond at the join is no overlap
    assert join_overlaps(joined, 8, case) == 0
    numbers = joined[:, col_resSeq].astype(int)
    names = joined[:, col_name]
    pushed = np.array(joined, copy=True)
    atom = np.flatnonzero((numbers == new) & (names == 'CA'))[0]
    pushed[atom, cols_coords] = joined[(numbers == old) & (names == 'CA')][0, cols_coords]  # noqa: E501
    assert join_overlaps(pushed, 8, case) > 0


def test_build_in_fragments_join_overlaps(expl_rows, monkeypatch):
    """Test a fragment that overlaps the chain across the join is redrawn."""
    first = int(expl_rows[:, col_resSeq].astype(int).min())
    last = int(expl_rows[:, col_resSeq].astype(int).max())
    pieces = [(first, first + 8), (first + 7, last)]
    calls = []

    def draw(k):
        calls.append(k)
        return _moved(_chain_rows(expl_rows, *pieces[k]), len(calls)), 0.0

    found = iter([1, 0])
    monkeypatch.setattr(ldrs_helper, 'join_overlaps', lambda *args: next(found))  # noqa: E501
    rows, _, rebuilt = build_in_fragments(draw, 2, 'C-IDR')
    assert calls == [0, 1, 1]
    assert rebuilt == 0
    assert len(rows) == len(_chain_rows(expl_rows, first, last))


def test_build_in_fragments_stop(expl_rows):
    """Test no fragment is added once `stop` tells the chain is done."""
    first = int(expl_rows[:, col_resSeq].astype(int).min())
    pieces = [(first, first + 5), (first + 4, first + 10), (first + 9, first + 15)]  # noqa: E501
    calls = []

    def draw(k):
        calls.append(k)
        return _moved(_chain_rows(expl_rows, *pieces[k]), k + 1), float(k)

    def stop(rows):
        return rows[:, col_resSeq].astype(int).max() >= 10

    rows, energy, _ = build_in_fragments(draw, 3, 'C-IDR', stop=stop)
    assert calls == [0, 1]
    assert energy == 1.0
    assert sorted(set(rows[:, col_resSeq].astype(int))) == list(range(1, 12))


@pytest.mark.parametrize('environment_overlaps', [True, False])
def test_build_in_fragments_environment_overlaps(
        expl_rows,
        monkeypatch,
        environment_overlaps,
        ):
    """Test the overlap rule can leave out the environment."""
    first = int(expl_rows[:, col_resSeq].astype(int).min())
    last = int(expl_rows[:, col_resSeq].astype(int).max())
    pieces = [(first, first + 8), (first + 7, last)]
    environment = np.array(expl_rows, copy=True)
    xyz = environment[:, cols_coords].astype(float) + 100
    environment[:, cols_coords] = np.round(xyz, 3).astype(str)
    environment[:, col_segid] = 'ENV'
    seen = []

    def overlaps(rows, fixed, *args, **kwargs):
        seen.append(bool(np.any(fixed[:, col_segid] == 'ENV')))
        return 0

    def draw(k):
        return _moved(_chain_rows(expl_rows, *pieces[k]), k + 1), 0.0

    monkeypatch.setattr(ldrs_helper, 'attachment_overlaps', overlaps)
    build_in_fragments(
        draw,
        2,
        'C-IDR',
        environment=environment,
        environment_overlaps=environment_overlaps,
        )
    assert seen == [environment_overlaps]


def test_with_environment_keeps_environment_residues(expl_rows):
    """Test a template residue numbered like the chain's first is kept."""
    first = int(expl_rows[:, col_resSeq].astype(int).min())
    chain = _chain_rows(expl_rows, first, first + 8)
    # an ion with the chain ID and number of the chain's first residue
    ion = np.array(chain[:1], copy=True)
    ion[0, col_name] = ion[0, col_resName] = 'MG'
    ion[0, cols_coords] = ['50.0', '50.0', '50.0']
    target = attachment_anchor(chain, 'N-IDR')
    fixed = ldrs_helper._with_environment(chain, ion, 'N-IDR')
    kept, _ = ldrs_helper._without_replaced(fixed, 'N-IDR', target)
    assert 'MG' in kept[:, col_name]
    assert len(kept) == len(fixed) - np.count_nonzero(chain[:, col_resSeq] == '1')  # noqa: E501


def test_attachment_overlaps(expl_rows):
    """Test the residue next to the anchor is checked, the peptide bond not."""
    resseq = expl_rows[:, col_resSeq].astype(int)
    first = int(resseq.min())
    folded = Structure(data_folder / 'EXPL_A.pdb')
    folded.build()
    folded._data_array = expl_rows[resseq >= first + 6]
    rows = folded.data_array
    target = np.array([
        rows[(rows[:, col_resSeq].astype(int) == res) & (rows[:, col_name] == name)][0, cols_coords]  # noqa: E501
        for res, name in ((first + 6, 'C'), (first + 7, 'N'), (first + 7, 'CA'))  # noqa: E501
        ]).astype(float)
    fragment = _chain_rows(expl_rows, first, first + 7)
    aligned = align_coords(_moved(fragment, 2), target, 'N-IDR')
    assert attachment_overlaps(aligned, rows, 'N-IDR', target) == 0
    # put a side-chain atom of the IDR residue next to the folded domain,
    # which count_clashes leaves out, on a folded C-alpha
    names = aligned[:, col_name]
    numbers = aligned[:, col_resSeq].astype(int)
    pushed = np.array(aligned, copy=True)
    atom = np.flatnonzero((numbers == 6) & (names == 'CB'))[0]
    pushed[atom, cols_coords] = rows[rows[:, col_name] == 'CA'][3, cols_coords]  # noqa: E501
    assert attachment_overlaps(pushed, rows, 'N-IDR', target) > 0
    assert count_clashes(pushed, folded, 'N-IDR', 10**6)[0] == count_clashes(aligned, folded, 'N-IDR', 10**6)[0]  # noqa: E501


def test_attachment_overlaps_alternate_locations(expl_rows):
    """Test alternate locations of the bonded residue are no overlaps."""
    resseq = expl_rows[:, col_resSeq].astype(int)
    first = int(resseq.min())
    rows = expl_rows[resseq >= first + 6]
    numbers = rows[:, col_resSeq].astype(int)
    target = np.array([
        rows[(numbers == res) & (rows[:, col_name] == name)][0, cols_coords]
        for res, name in ((first + 6, 'C'), (first + 7, 'N'), (first + 7, 'CA'))  # noqa: E501
        ]).astype(float)
    aligned = align_coords(_moved(_chain_rows(expl_rows, first, first + 7), 2), target, 'N-IDR')  # noqa: E501
    # a second location of the residue bonded to the IDR, 0.7 A away,
    # after the first one
    bonded = rows[numbers == first + 7]
    alternate = np.array(bonded, copy=True)
    xyz = bonded[:, cols_coords].astype(float) + 0.4
    alternate[:, cols_coords] = np.round(xyz, 3).astype(str)
    rows = np.concatenate((rows[numbers <= first + 7], alternate, rows[numbers > first + 7]))  # noqa: E501
    assert attachment_overlaps(aligned, rows, 'N-IDR', target) == 0


@pytest.mark.parametrize('case', ['N-IDR', 'C-IDR'])
def test_build_in_fragments(expl_rows, case, monkeypatch):
    """Test fragments are attached in order and dead ends are rebuilt."""
    first = int(expl_rows[:, col_resSeq].astype(int).min())
    last = int(expl_rows[:, col_resSeq].astype(int).max())
    if case == 'C-IDR':
        pieces = [(first, first + 5), (first + 4, first + 10), (first + 9, last)]  # noqa: E501
    else:
        pieces = [(first + 9, last), (first + 4, first + 10), (first, first + 5)]  # noqa: E501
    calls = []

    def draw(k):
        calls.append(k)
        rows = _moved(_chain_rows(expl_rows, *pieces[k]), len(calls))
        if k == 1 and calls.count(1) == 1:
            rows[:, col_segid] = 'BAD'
        return rows, float(k)

    # the first conformer of fragment 1 leaves an end nothing fits onto
    real = ldrs_helper.attachment_overlaps

    def overlaps(rows, fixed, *args, **kwargs):
        return 99 if np.any(fixed[:, col_segid] == 'BAD') else real(rows, fixed, *args, **kwargs)  # noqa: E501

    monkeypatch.setattr(ldrs_helper, 'attachment_overlaps', overlaps)
    rows, energy, rebuilt = build_in_fragments(draw, 3, case, max_trials=3)  # noqa: E501
    assert calls == [0, 1, 2, 2, 2, 1, 2]
    assert rebuilt == 1
    assert energy == 3.0
    whole = _chain_rows(expl_rows, first, last)
    assert list(rows[:, col_name]) == list(whole[:, col_name])
    assert list(rows[:, col_resSeq]) == list(whole[:, col_resSeq])
    # the first fragment stays where it was built: compare internal
    # distances
    xyz = rows[:, cols_coords].astype(float)
    ref = whole[:, cols_coords].astype(float)
    d1 = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    d2 = np.linalg.norm(ref[:, None] - ref[None], axis=-1)
    assert np.allclose(d1, d2, atol=0.1)
