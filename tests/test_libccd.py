"""Test torsion-space loop closure (libccd)."""
import math

import numpy as np
import pytest

from idpconfgen.libs import libccd


def _place(a, b, c, bond, angle, torsion):
    """Place a point from three previous ones (NeRF)."""
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([
        -bond * math.cos(angle),
        bond * math.sin(angle) * math.cos(torsion),
        bond * math.sin(angle) * math.sin(torsion),
        ])
    return c + d[0] * bc + d[1] * m + d[2] * n


def _build(phi, psi, omega=None):
    """Build N, CA, C, O of a chain from torsions in radians."""
    size = len(phi)
    omega = np.full(size, math.pi) if omega is None else omega
    X = np.zeros((size, 4, 3))
    X[0, 1] = (1.458, 0, 0)
    angle = math.radians(111.2)
    X[0, 2] = X[0, 1] + 1.525 * np.array([-math.cos(angle), math.sin(angle), 0])  # noqa: E501
    for i in range(1, size):
        X[i, 0] = _place(*X[i - 1, :3], 1.329, math.radians(116.2), psi[i - 1])  # noqa: E501
        X[i, 1] = _place(X[i - 1, 1], X[i - 1, 2], X[i, 0], 1.458, math.radians(121.7), omega[i - 1])  # noqa: E501
        X[i, 2] = _place(X[i - 1, 2], X[i, 0], X[i, 1], 1.525, math.radians(111.2), phi[i])  # noqa: E501
    for i in range(size):
        X[i, 3] = _place(*X[i, :3], 1.231, math.radians(120.8), psi[i] + math.pi)  # noqa: E501
    return X


def _loop_case(nidr, seed, noise=8.0):
    """Return a reference chain and a perturbed copy to close onto it."""
    rng = np.random.default_rng(seed)
    size = nidr + 4
    phi = np.radians(rng.choice([-65, -120, -75], size) + rng.normal(0, 10, size))  # noqa: E501
    psi = np.radians(rng.choice([-40, 130, 145], size) + rng.normal(0, 10, size))  # noqa: E501
    reference = _build(phi, psi)
    candidate = _build(
        phi + np.radians(rng.normal(0, noise, size)),
        psi + np.radians(rng.normal(0, noise, size)),
        )
    return reference, candidate


def _random_frame(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([
        [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],  # noqa: E501
        [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],  # noqa: E501
        [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d],  # noqa: E501
        ])


def test_kabsch_recovers_rigid_motion():
    """Test superposition returns the applied rotation and translation."""
    rng = np.random.default_rng(1)
    points = rng.normal(size=(6, 3))
    rot = _random_frame(rng)
    moved = points @ rot.T + np.array([1.0, -2.0, 3.0])
    rotation, translation, rmsd = libccd.kabsch(points, moved)
    assert rmsd < 1e-10
    assert np.allclose(rotation, rot)
    assert np.allclose(points @ rotation.T + translation, moved)


def test_dihedral_and_backbone_phi_psi():
    """Test torsions measured on a built chain equal the input ones."""
    rng = np.random.default_rng(2)
    phi = np.radians(rng.uniform(-170, 170, 10))
    psi = np.radians(rng.uniform(-170, 170, 10))
    chain = _build(phi, psi)
    measured = libccd.backbone_phi_psi(chain, 1, 8)
    assert np.allclose(libccd.wrap_angle(measured[:, 0] - phi[1:9]), 0, atol=1e-8)  # noqa: E501
    assert np.allclose(libccd.wrap_angle(measured[:, 1] - psi[1:9]), 0, atol=1e-8)  # noqa: E501


@pytest.mark.parametrize(
    'nres,expected',
    [(10, (0, 1, 2)), (40, (2, 1, 3)), (200, (3, 2, 4)), (5, (0, 1, 2))],
    )
def test_length_schedule(nres, expected):
    """Test pivot schedule by loop length."""
    assert libccd.length_schedule(nres) == expected


def test_rotate_torsions_changes_each_torsion_exactly():
    """Test simultaneous rotations equal sequential C-terminal-first ones."""
    rng = np.random.default_rng(3)
    phi = np.radians(rng.uniform(-170, 170, 12))
    psi = np.radians(rng.uniform(-170, 170, 12))
    chain = _build(phi, psi)
    X = chain.reshape(-1, 3).copy()
    search = libccd._LoopSearch(
        X.copy(), 12, 2, 9, X[40:43].copy(),
        libccd.CCDParameters(), None, np.zeros(12, dtype=bool),
        )
    pivots = search.pivots(2, 9, 1)
    keys, ia, ib, starts = libccd._pivot_arrays(pivots)
    anchors = X[ia].copy()
    axes = X[ib] - anchors
    axes /= np.linalg.norm(axes, axis=1)[:, None]
    step = np.radians(rng.normal(0, 20, len(pivots)))

    together = X.copy()
    libccd.rotate_torsions(together, anchors, axes, starts, step)

    one_by_one = X.copy()
    for (_key, a_i, b_i, s_i), theta in zip(pivots, step):
        a = one_by_one[a_i]
        u = one_by_one[b_i] - a
        u /= np.linalg.norm(u)
        rot = libccd.rotation_matrices(u[None, :], np.array([theta]))[0]
        one_by_one[s_i:] = (one_by_one[s_i:] - a) @ rot.T + a
    assert np.allclose(together, one_by_one, atol=1e-9)

    before = libccd.backbone_phi_psi(chain, 2, 9)
    after = libccd.backbone_phi_psi(together.reshape(12, 4, 3), 2, 9)
    expected = np.zeros_like(before)
    for (residue, kind), theta in zip(keys, step):
        expected[residue - 2, kind] = theta
    assert np.allclose(libccd.wrap_angle(after - before), expected, atol=1e-9)

    # bond lengths are untouched
    flat_before = X[:-1] - X[1:]
    flat_after = together[:-1] - together[1:]
    assert np.allclose(
        np.linalg.norm(flat_before, axis=1),
        np.linalg.norm(flat_after, axis=1),
        )


def test_torsion_displacement_is_the_linearisation():
    """Test linearised displacements against a tiny finite rotation."""
    rng = np.random.default_rng(4)
    chain = _build(
        np.radians(rng.uniform(-170, 170, 10)),
        np.radians(rng.uniform(-170, 170, 10)),
        )
    X = chain.reshape(-1, 3)
    starts = np.array([6, 11, 19, 30])
    ia = starts - 2
    ib = starts - 1
    anchors = X[ia].copy()
    axes = X[ib] - anchors
    axes /= np.linalg.norm(axes, axis=1)[:, None]
    step = np.radians(rng.normal(0, 1, len(starts))) * 1e-3
    atoms = np.arange(len(X))
    moved = X.copy()
    libccd.rotate_torsions(moved, anchors, axes, starts, step)
    linear = libccd.torsion_displacement(X, atoms, anchors, axes, starts, step)
    assert np.allclose(moved - X, linear, atol=1e-9)
    assert np.all(linear[:6] == 0)


def test_clash_environment_pairs_match_brute_force():
    """Test vectorised loop-environment pairs against all distances."""
    rng = np.random.default_rng(5)
    env_xyz = rng.uniform(0, 12, size=(300, 3))
    env_radii = rng.uniform(1.2, 1.9, 300)
    lower = np.zeros(300, dtype=bool)
    lower[:5] = True
    env = libccd.ClashEnvironment(env_xyz, env_radii, lower_mask=lower)
    loop_xyz = rng.uniform(0, 12, size=(40, 3))
    loop_radii = rng.uniform(1.4, 1.7, 40)
    residue = np.repeat(np.arange(2, 12), 4)

    k, e, cut, dist = env.pairs(loop_xyz, loop_radii, residue, 2, 11, 0.65, 0.5)  # noqa: E501
    d = np.linalg.norm(loop_xyz[:, None] - env_xyz[None], axis=2)
    c = 0.65 * (loop_radii[:, None] + env_radii[None])
    near = d < c + 0.5
    near[(residue == 2)[:, None] & lower[None, :]] = False
    assert set(zip(k, e)) == set(zip(*np.nonzero(near)))
    assert np.allclose(dist, d[k, e])
    assert np.allclose(cut, c[k, e])

    k2, e2, cut2 = env.overlaps(loop_xyz, loop_radii, residue, 2, 11, 0.65, 0.5)  # noqa: E501
    assert np.allclose(cut2, cut + 0.5)
    assert env.count_overlaps(loop_xyz, loop_radii, residue, 2, 11, 0.65) == int(np.sum(d[~((residue == 2)[:, None] & lower[None, :])] < c[~((residue == 2)[:, None] & lower[None, :])]))  # noqa: E501


def test_self_overlaps():
    """Test nonlocal backbone contacts of a chain folded onto itself."""
    params = libccd.CCDParameters()
    chain = _build(np.radians(np.full(12, -60.0)), np.radians(np.full(12, -45.0)))  # noqa: E501
    assert libccd.count_self_overlaps(chain, 0, 11, params) == 0
    squeezed = chain.copy()
    squeezed[9] = squeezed[2] + 0.3  # residue 9 almost on top of residue 2
    assert libccd.count_self_overlaps(squeezed, 0, 11, params) > 0


def test_contact_model():
    """Test amide H and C-beta atoms and bond-separation exclusions."""
    params = libccd.CCDParameters()
    sequence = 'AAGAPAAAAAAA'
    chain = _build(np.radians(np.full(12, -60.0)), np.radians(np.full(12, -45.0)))  # noqa: E501
    model = libccd.ContactModel(2, 9, sequence)
    flat = chain.reshape(-1, 3)
    coords = model.coordinates(flat)
    # 8 residues: 4 backbone atoms each, H except Pro 6, CB except Gly 4
    assert len(coords) == 8 * 4 + 7 + 7
    assert np.allclose(coords[:32], flat[8:40])
    n_h = np.linalg.norm(coords[32:39] - flat[4 * model.h_res], axis=1)
    assert np.allclose(n_h, 0.9)
    ca_cb = np.linalg.norm(coords[39:] - flat[4 * model.cb_res + 1], axis=1)
    assert np.all((ca_cb > 1.45) & (ca_cb < 1.6))
    # an ideal helix has no overlaps; only pairs > 3 bonds apart count
    assert len(model.self_pairs(coords, params)[0]) == 0
    i, j, cut, dist = model.self_pairs(coords, params, margin=3.0)
    bonds = np.abs(model.chain[i] - model.chain[j]) + model.branch[i] + model.branch[j]  # noqa: E501
    assert np.all(bonds > 3)
    assert np.all(model.follows[i] <= model.follows[j])
    # O(i) and N(i + 2) on top of each other: a local clash
    squeezed = flat.copy()
    squeezed[4 * 6] = squeezed[4 * 4 + 3] + 0.5
    i, j, cut, dist = model.self_pairs(model.coordinates(squeezed), params)
    assert len(i) > 0
    assert np.all(dist < cut)


def test_damped_step_respects_cap():
    """Test frozen torsions never exceed the cap."""
    rng = np.random.default_rng(6)
    A = rng.normal(size=(9, 20))
    b = rng.normal(size=9) * 10
    total = np.zeros(20)
    total[:4] = math.radians(59.0)
    step = libccd._damped_step(A, b, 0.09, 0.0, total, math.radians(60.0))
    assert np.all(np.abs(total + step) <= math.radians(60.0) + 1e-12)
    free = libccd._damped_step(A, b, 0.09, 0.0, np.zeros(20), 10.0)
    gram = A.T @ A + 0.09 * np.eye(20)
    assert np.allclose(free, np.linalg.solve(gram, A.T @ b))


def test_junction_failure():
    """Test junction windows."""
    params = libccd.CCDParameters()
    good = {
        f'{side}_{key}': value
        for side in ('lower', 'upper')
        for key, value in (
            ('c_n', 1.33), ('ca_c_n', 116.0), ('o_c_n', 123.0),
            ('c_n_ca', 121.0), ('omega', 179.0),
            )
        }
    assert libccd.junction_failure(good, params) is None
    bad = dict(good, upper_c_n=1.9)
    assert libccd.junction_failure(bad, params) == 'junction_upper_distance'
    # 0.07 A too long, accepted by the old 1.18-1.48 A window
    stretched = dict(good, upper_c_n=1.40)
    assert libccd.junction_failure(stretched, params) == 'junction_upper_distance'  # noqa: E501
    bent = dict(good, upper_c_n_ca=109.0)
    assert libccd.junction_failure(bent, params) == 'junction_upper_angle'
    assert libccd.junction_failure(dict(good, lower_o_c_n=float('nan')), params) is None  # noqa: E501
    cis = dict(good, lower_omega=5.0)
    assert libccd.junction_failure(cis, params) == 'junction_lower_omega'
    assert libccd.junction_failure(cis, params, lower_is_pro=True) is None
    assert libccd.junction_failure(dict(good, upper_omega=-165.0), params) is None  # noqa: E501
    twisted = dict(good, upper_omega=155.0)
    assert libccd.junction_failure(twisted, params) == 'junction_upper_omega'


@pytest.mark.parametrize(
    'omega,next_is_pro,allow,expected',
    [
        (180.0, False, False, 'trans'),
        (-179.0, False, False, 'trans'),
        (161.0, False, False, 'trans'),
        (159.0, False, False, 'twisted'),
        (-150.0, True, False, 'twisted'),
        (90.0, True, False, 'twisted'),
        (5.0, True, False, 'cis'),
        (-19.0, True, False, 'cis'),
        (5.0, False, False, 'nonpro_cis'),
        (5.0, False, True, 'cis'),
        (25.0, True, True, 'twisted'),
        (540.0, False, False, 'trans'),
        ],
    )
def test_omega_state(omega, next_is_pro, allow, expected):
    """Test the cis/trans basins of omega."""
    state = libccd.omega_state(
        omega,
        next_is_pro=next_is_pro,
        tolerance=20.0,
        allow_nonpro_cis=allow,
        )
    assert state == expected


def test_peptide_omegas_and_failures():
    """Test omega measurement and the count of disallowed peptides."""
    rng = np.random.default_rng(11)
    size = 10
    phi = np.radians(rng.uniform(-150, -50, size))
    psi = np.radians(rng.uniform(-60, 150, size))
    omega = np.radians(rng.uniform(170, 190, size))
    omega[3] = math.radians(2.0)  # cis
    omega[6] = math.radians(120.0)  # twisted
    chain = _build(phi, psi, omega)
    measured = libccd.peptide_omegas(chain, 0, size - 1)
    assert np.allclose(
        libccd.wrap_angle(np.radians(measured) - omega[:-1]), 0, atol=1e-8,
        )
    params = libccd.CCDParameters()
    next_is_pro = [False] * (size - 1)
    assert libccd.omega_failures(measured, next_is_pro, params) == {
        'nonpro_cis': 1, 'twisted': 1,
        }
    next_is_pro[3] = True
    assert libccd.omega_failures(measured, next_is_pro, params) == {
        'nonpro_cis': 0, 'twisted': 1,
        }
    assert len(libccd.peptide_omegas(chain, 4, 4)) == 0


def test_backbone_geometry_and_violations():
    """Test bond lengths and angles against the geometry windows."""
    chain = _build(np.radians(np.full(8, -65.0)), np.radians(np.full(8, 140.0)))  # noqa: E501
    geometry = libccd.backbone_geometry(chain, 0, 7)
    assert set(geometry) == set(libccd.GEOMETRY_WINDOWS)
    assert np.allclose(geometry['C_N'], 1.329)
    assert np.allclose(geometry['N_CA_C'], 111.2)
    assert np.allclose(geometry['O_C_N'], 360.0 - 116.2 - 120.8)
    assert len(geometry['C_N']) == 7 and len(geometry['N_CA']) == 8
    assert libccd.geometry_violations(chain, 0, 7) == 0
    broken = chain.copy()
    bond = chain[5, 0] - chain[4, 2]
    broken[5:] += 0.5 * bond / np.linalg.norm(bond)  # a chain break before 5
    assert libccd.geometry_violations(broken, 0, 7) >= 1
    # only the range is checked
    assert libccd.geometry_violations(broken, 5, 7) == 0


def test_idealize_backbone_keeps_torsions():
    """Test idealization changes bond geometry but no torsion."""
    rng = np.random.default_rng(12)
    size = 14
    sequence = 'AGPKLMNAAGSTPE'
    phi = np.radians(rng.uniform(-170, 170, size))
    psi = np.radians(rng.uniform(-170, 170, size))
    omega = np.radians(rng.normal(180, 6, size))
    chain = _build(phi, psi, omega)
    # a fragment junction: residues from 7 on shifted and turned
    rot = _random_frame(rng)
    pivot = chain[7, 0]
    distorted = chain.copy()
    distorted[7:] = (chain[7:] - pivot) @ (0.2 * rot + 0.8 * np.eye(3)).T + pivot  # noqa: E501
    ideal = libccd.idealize_backbone(distorted, sequence)
    assert libccd.geometry_violations(ideal, 0, size - 1) == 0
    for torsions in (
            lambda bb: libccd.backbone_phi_psi(bb, 1, size - 2),
            lambda bb: np.radians(libccd.peptide_omegas(bb, 0, size - 1)),
            ):
        before, after = torsions(distorted), torsions(ideal)
        assert np.allclose(libccd.wrap_angle(after - before), 0, atol=1e-8)
    geometry = libccd.backbone_geometry(ideal, 0, size - 1)
    for key in ('C_N', 'CA_C_N', 'O_C_N', 'C_N_CA', 'C_O', 'CA_C_O'):
        assert np.allclose(geometry[key], libccd.IDEAL_GEOMETRY[key], atol=0.2)  # noqa: E501
    # superimposed onto residue 1; an ideal chain is a fixed point
    assert np.allclose(ideal[1, :3], distorted[1, :3], atol=0.05)
    again = libccd.idealize_backbone(ideal, sequence)
    assert np.allclose(again, ideal, atol=1e-8)


def _loop_torsions(nidr, seed):
    rng = np.random.default_rng(seed)
    size = nidr + 4
    phi = np.radians(rng.choice([-65, -120, -75], size) + rng.normal(0, 10, size))  # noqa: E501
    psi = np.radians(rng.choice([-40, 130, 145], size) + rng.normal(0, 10, size))  # noqa: E501
    noise = (
        np.radians(rng.normal(0, 8.0, size)),
        np.radians(rng.normal(0, 8.0, size)),
        )
    return phi, psi, noise


def test_close_loop_checks_source_omega():
    """Test a source with a non-proline cis peptide is not closed."""
    phi, psi, (dphi, dpsi) = _loop_torsions(12, 1)
    reference = _build(phi, psi)
    omega = np.full(16, math.pi)
    omega[6] = 0.0
    candidate = _build(phi + dphi, psi + dpsi, omega)
    result = libccd.close_loop(candidate, 'A' * 16, reference[1], reference[14])  # noqa: E501
    assert not result.accepted
    assert result.reason == 'omega'
    assert result.stats['source_nonpro_cis'] == 1
    # before a proline the cis peptide is allowed
    sequence = 'A' * 7 + 'P' + 'A' * 8
    result = libccd.close_loop(candidate, sequence, reference[1], reference[14])  # noqa: E501
    assert result.reason != 'omega'
    assert result.stats['source_nonpro_cis'] == 0
    params = libccd.CCDParameters(allow_nonpro_cis=True)
    result = libccd.close_loop(candidate, 'A' * 16, reference[1], reference[14], params=params)  # noqa: E501
    assert result.reason != 'omega'


def test_close_loop_idealizes_broken_sources():
    """Test sources outside the geometry windows are rebuilt or rejected."""
    phi, psi, (dphi, dpsi) = _loop_torsions(30, 4)
    size = 34
    reference = _build(phi, psi)
    candidate = _build(phi + dphi, psi + dpsi)
    broken = candidate.copy()
    broken[15:] += np.array([0.0, 0.4, 0.3])  # stretched C(14)-N(15)
    # the rigid tail carries OXT and moves with residue 33
    oxt = broken[-1, 2] + np.array([1.2, 0.3, 0.1])
    never = libccd.CCDParameters(idealize='never')
    result = libccd.close_loop(broken, 'A' * size, reference[1], reference[32], params=never)  # noqa: E501
    assert result.reason == 'backbone_geometry'
    assert result.stats['source_geometry_violations'] > 0

    result = libccd.close_loop(
        broken, 'A' * size, reference[1], reference[32], extra=oxt,
        )
    assert result.stats['idealized'] == 1
    assert result.accepted, result.reason
    closed = result.backbone
    assert libccd.geometry_violations(closed, 1, size - 3) == 0
    # OXT kept its place relative to the (rebuilt) last residue
    assert np.isclose(
        np.linalg.norm(result.extra[0] - closed[-1, 2]),
        np.linalg.norm(oxt - broken[-1, 2]),
        atol=0.01,
        )
    # clean sources are never rebuilt with 'auto'
    clean = libccd.close_loop(candidate, 'A' * size, reference[1], reference[32])  # noqa: E501
    assert clean.stats['idealized'] == 0
    always = libccd.CCDParameters(idealize='always')
    rebuilt = libccd.close_loop(candidate, 'A' * size, reference[1], reference[32], params=always)  # noqa: E501
    assert rebuilt.stats['idealized'] == 1
    with pytest.raises(ValueError):
        libccd.close_loop(
            candidate, 'A' * size, reference[1], reference[32],
            params=libccd.CCDParameters(idealize='sometimes'),
            )


def test_close_loop_checks_the_bond_to_the_upper_anchor():
    """Test the source bond that becomes the upper junction is checked."""
    phi, psi, (dphi, dpsi) = _loop_torsions(30, 4)
    reference = _build(phi, psi)
    candidate = _build(phi + dphi, psi + dpsi)
    # stretch C(31)-N(32), the last IDR residue to the anchor residue
    broken = candidate.copy()
    bond = candidate[32, 0] - candidate[31, 2]
    broken[32:] += 0.25 * bond / np.linalg.norm(bond)
    never = libccd.CCDParameters(idealize='never')
    result = libccd.close_loop(broken, 'A' * 34, reference[1], reference[32], params=never)  # noqa: E501
    assert result.reason == 'backbone_geometry'
    result = libccd.close_loop(broken, 'A' * 34, reference[1], reference[32])  # noqa: E501
    assert result.stats['idealized'] == 1
    assert result.accepted, result.reason
    assert abs(result.stats['upper_c_n'] - libccd.IDEAL_GEOMETRY['C_N']) < 0.05  # noqa: E501
    # a cis peptide bond to the anchor residue is caught on the source
    omega = np.full(34, math.pi)
    omega[31] = 0.0
    cis = _build(phi + dphi, psi + dpsi, omega)
    result = libccd.close_loop(cis, 'A' * 34, reference[1], reference[32])
    assert result.reason == 'omega'


def test_close_loop_rotation_maps_residue_zero():
    """Test the returned superposition moves residue-0 atoms correctly."""
    phi, psi, (dphi, dpsi) = _loop_torsions(30, 8)
    reference = _build(phi, psi)
    candidate = _build(phi + dphi, psi + dpsi)
    broken = candidate.copy()
    broken[20:] += np.array([0.5, 0.0, 0.0])
    for source, atol in ((candidate, 1e-6), (broken, 0.05)):
        result = libccd.close_loop(source, 'A' * 34, reference[1], reference[32])  # noqa: E501
        moved = source[0] @ result.rotation.T + result.translation
        assert np.allclose(moved, result.backbone[0], atol=atol)


def test_ramachandran_reference():
    """Test empirical phi/psi occupancy."""
    rng = np.random.default_rng(7)
    angles = np.zeros((5000, 3))
    angles[:, 1] = np.radians(rng.normal(-65, 8, 5000))
    angles[:, 2] = np.radians(rng.normal(-40, 8, 5000))
    rama = libccd.RamachandranReference(angles, 'A' * 5000)
    assert not rama.is_outlier('A', math.radians(-65), math.radians(-40))
    assert rama.is_outlier('A', math.radians(60), math.radians(-120))
    # residue types without data use the pooled reference
    assert not rama.is_outlier('W', math.radians(-65), math.radians(-40))


@pytest.mark.parametrize('solver', libccd.SOLVERS)
@pytest.mark.parametrize('nidr,seed', [(12, 1), (12, 7), (30, 4), (30, 8)])
def test_close_loop(solver, nidr, seed):
    """Test closing a perturbed chain back onto its anchors."""
    reference, candidate = _loop_case(nidr, seed)
    size = nidr + 4
    rng = np.random.default_rng(seed)
    rot = _random_frame(rng)
    candidate = candidate @ rot.T + 5.0  # arbitrary initial frame
    params = libccd.CCDParameters(solver=solver)
    result = libccd.close_loop(
        candidate,
        'A' * size,
        reference[1],
        reference[size - 2],
        params=params,
        )
    assert result.accepted, result.reason
    assert result.stats['solver'] == solver
    assert result.stats['final_rmsd'] <= params.target_rmsd
    assert result.stats['phipsi_rms_delta'] <= params.max_phipsi_rms_delta
    closed = result.backbone
    assert np.allclose(closed[1, :3], reference[1, :3], atol=1e-6)
    assert np.allclose(closed[1, 3], reference[1, 3], atol=1e-6)

    # rigid rotations keep every bond length
    def bonds(chain):
        flat = chain.reshape(-1, 3)
        return np.linalg.norm(flat[1:] - flat[:-1], axis=1)
    assert np.allclose(bonds(closed), bonds(candidate), atol=1e-6)


def test_close_loop_dls_spreads_smaller_changes_than_ccd():
    """Test the least-squares solver changes torsions less than CCD."""
    deltas = {}
    for solver in libccd.SOLVERS:
        values = []
        for seed in (0, 4, 5, 10):
            reference, candidate = _loop_case(30, seed)
            result = libccd.close_loop(
                candidate,
                'A' * 34,
                reference[1],
                reference[32],
                params=libccd.CCDParameters(solver=solver),
                )
            assert result.accepted
            values.append(result.stats['phipsi_rms_delta'])
        deltas[solver] = np.mean(values)
    assert deltas['dls'] < deltas['ccd']


@pytest.mark.parametrize('seed', [6, 8, 10])
def test_close_loop_routes_around_environment(seed):
    """Test DLS closes a loop around atoms sitting on its natural path."""
    reference, candidate = _loop_case(30, seed, noise=15.0)
    size = 34
    # environment atoms on the middle of the reference loop
    wall = reference[12:22, 1] + 0.2
    env = libccd.ClashEnvironment(wall, np.full(len(wall), 1.7))
    result = libccd.close_loop(
        candidate,
        'A' * size,
        reference[1],
        reference[size - 2],
        environment=env,
        )
    assert result.accepted, result.reason
    coords = result.backbone[2:size - 2].reshape(-1, 3)
    residue = np.repeat(np.arange(2, size - 2), 4)
    overlaps = env.count_overlaps(
        coords,
        np.tile(libccd._BB_RADII, size - 4),
        residue,
        2,
        size - 3,
        0.65,
        )
    assert overlaps == 0
    # the reference path itself would overlap the wall
    ref_coords = reference[2:size - 2].reshape(-1, 3)
    assert env.count_overlaps(
        ref_coords,
        np.tile(libccd._BB_RADII, size - 4),
        residue,
        2,
        size - 3,
        0.65,
        ) > 0


def test_close_loop_invalid_input():
    """Test inputs that cannot be closed."""
    reference, candidate = _loop_case(12, 1)
    result = libccd.close_loop(candidate, 'A' * 3, reference[1], reference[14])  # noqa: E501
    assert not result.accepted
    assert result.reason == 'invalid_input'
    with pytest.raises(ValueError):
        libccd.close_loop(
            candidate, 'A' * 16, reference[1], reference[14],
            params=libccd.CCDParameters(solver='other'),
            )


def _lysines(chain, residues, chis):
    """Add ideal C-beta and lysine side chains to backbone residues."""
    names, numbers, coords = [], [], []
    for r, (n, ca, c, o) in enumerate(chain):
        atoms = [('N', n), ('CA', ca), ('C', c), ('O', o)]
        if r in residues:
            b, cc = ca - n, c - ca
            a = np.cross(b, cc)
            cb = -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca
            chi = chis[residues.index(r)]
            cg = _place(n, ca, cb, 1.52, math.radians(114.0), chi[0])
            cd = _place(ca, cb, cg, 1.52, math.radians(111.0), chi[1])
            ce = _place(cb, cg, cd, 1.52, math.radians(111.0), chi[2])
            nz = _place(cg, cd, ce, 1.49, math.radians(111.0), chi[3])
            atoms += [('CB', cb), ('CG', cg), ('CD', cd), ('CE', ce), ('NZ', nz)]  # noqa: E501
        for name, xyz in atoms:
            names.append(name)
            numbers.append(r + 1)
            coords.append(xyz)
    return np.array(names), np.array(numbers), np.array(coords)


def test_covalent_bonds_and_chi_pivots():
    """Test bonds stay within residues and chi rotations move their tips."""
    chain = _build(np.radians(np.full(6, -65.0)), np.radians(np.full(6, 140.0)))  # noqa: E501
    chi = [np.radians((180, 180, 180, 180))] * 2
    names, numbers, coords = _lysines(chain, [1, 4], chi)
    # an atom of residue 3 placed on top of residue 2's CA is not bonded
    coords = coords.copy()
    coords[np.flatnonzero((numbers == 4) & (names == 'O'))[0]] = (
        coords[np.flatnonzero((numbers == 3) & (names == 'CA'))[0]] + 0.3
        )
    bonds = libccd.covalent_bonds(coords, names, numbers)
    pairs = {tuple(sorted((int(i), int(j)))) for i, j in bonds}
    for i, j in pairs:
        same = numbers[i] == numbers[j]
        peptide = {names[i], names[j]} == {'C', 'N'} and abs(numbers[i] - numbers[j]) == 1  # noqa: E501
        assert same or peptide
    # 6 residues: 3 bonds each + 5 peptide bonds + 5 per lysine side chain,
    # minus the C=O of residue 4 whose O was moved away
    assert len(pairs) == 6 * 3 + 5 + 2 * 5 - 1
    pivots = libccd.chi_pivots(names, numbers, 'AKAAKA', [2, 5, 3], bonds)
    sizes = [len(p[3]) for p in pivots if p[0] == 2]
    assert sizes == [4, 3, 2, 1]
    assert {p[0] for p in pivots} == {2, 5}  # alanine 3 has no chi
    for residue, start, _end, moving in pivots:
        assert names[start] in ('CA', 'CB', 'CG', 'CD')
        assert np.all(numbers[moving] == residue)
        assert 'CB' not in names[moving]


def test_relax_side_chains_clears_overlaps():
    """Test chi rotations move side chains off fixed atoms only."""
    chain = _build(np.radians(np.full(8, -65.0)), np.radians(np.full(8, 140.0)))  # noqa: E501
    chis = [np.radians((180, 180, 180, 180)), np.radians((-60, 180, 180, 180))]  # noqa: E501
    names, numbers, coords = _lysines(chain, [2, 5], chis)
    tips = np.flatnonzero(np.isin(names, ('CE', 'NZ')))
    wall = coords[tips] + np.array([0.3, 0.0, 0.0])
    env = libccd.ClashEnvironment(wall, np.full(len(wall), 1.7))
    relaxed, stats = libccd.relax_side_chains(
        coords, names, numbers, 'AAKAAKAA', [3, 6], environment=env,
        )
    assert stats['overlaps_before'] > 0
    assert stats['overlaps_after'] == 0
    assert 0 < stats['max_chi_change'] <= 60.0 + 1e-6
    fixed = ~np.isin(names, ('CG', 'CD', 'CE', 'NZ'))
    assert np.allclose(relaxed[fixed], coords[fixed])
    # bond lengths and angles of the side chains are kept
    for residue in (3, 6):
        side = [np.flatnonzero((numbers == residue) & (names == n))[0] for n in ('CA', 'CB', 'CG', 'CD', 'CE', 'NZ')]  # noqa: E501
        for a, b in zip(side[:-1], side[1:]):
            assert np.isclose(
                np.linalg.norm(relaxed[a] - relaxed[b]),
                np.linalg.norm(coords[a] - coords[b]),
                )
    rad = 1.7
    for k in tips:
        r_k = libccd.vdW_radii_tsai_1999[names[k][0]]
        assert np.min(np.linalg.norm(wall - relaxed[k], axis=1)) >= 0.70 * (r_k + rad) - 1e-6  # noqa: E501

    # NaN rows are ignored and nothing moves without overlaps
    clean = coords.copy()
    clean[0] = np.nan
    same, stats = libccd.relax_side_chains(clean, names, numbers, 'AAKAAKAA', [3, 6])  # noqa: E501
    assert stats['overlaps_before'] == 0 and stats['iterations'] == 0
    assert np.allclose(same[1:], coords[1:]) and np.isnan(same[0]).all()
