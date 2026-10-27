"""Test LDRS helpers."""
from itertools import islice, product

import numpy as np
import pytest

from idpconfgen.core.build_definitions import build_bend_CA_C_O, distance_C_O
from idpconfgen.core.definitions import vdW_radii_tsai_1999
from idpconfgen.core.exceptions import IDPConfGenException
from idpconfgen.ldrs_helper import (
    _random_combinations,
    atom_vdw_radii,
    backbone_anchor,
    combinations_clash_check,
    create_all_combinations,
    create_combinations,
    inter_chain_cc,
    is_hydrogen,
    next_seeker,
    packing_context,
    prepare_ccd_context,
    psurgeon,
    splice_folded_rows,
    )
from idpconfgen.libs.libstructure import (
    Structure,
    col_chainID,
    col_name,
    col_resName,
    col_resSeq,
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
            row = (numbers == residue) & (rows[:, col_name] == atom)
            return rows[row][0, cols_coords].astype(float)

        c, o = xyz(term, 'C'), xyz(term, 'O')
        ca, n = xyz(term, 'CA'), xyz(term + 1, 'N')
        assert abs(np.linalg.norm(o - c) - distance_C_O) < 2e-3
        cos = np.dot(ca - c, o - c) / np.linalg.norm(ca - c) / np.linalg.norm(o - c)  # noqa: E501
        assert abs(np.arccos(cos) - 2 * build_bend_CA_C_O) < 2e-3
        normal = np.cross(ca - c, n - c)
        assert abs(np.dot(normal / np.linalg.norm(normal), o - c)) < 5e-3


def test_atom_vdw_radii_and_hydrogens():
    """Test radii follow count_clashes and fall back to atom names."""
    radii = atom_vdw_radii(['C', '', 'XX'], ['CA', '1HB', 'ZN9'], default=1.5)
    assert radii[0] == vdW_radii_tsai_1999['C']
    assert radii[1] == vdW_radii_tsai_1999['H']
    assert radii[2] == 1.5
    assert list(is_hydrogen(['H', '', 'C', ''], ['H', 'HA', 'CB', 'N'])) == [
        True, True, False, False,
        ]


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
