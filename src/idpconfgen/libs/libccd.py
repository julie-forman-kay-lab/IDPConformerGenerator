"""
Close loop-IDRs between two fixed anchors in torsion space.

A loop-IDR (``Linker-IDR`` in the ``ldrs`` subclient) bridges two fixed
folded domains. Instead of searching for two independently grown halves
whose ends happen to meet (``next_seeker``), this module closes every
IDPConformerGenerator backbone individually:

1.  A backbone is built for the IDR plus two overlapping folded residues on
    each side (the same sequence ``next_seeker`` uses). Its bond lengths
    and angles must be within :data:`GEOMETRY_WINDOWS` (otherwise it is
    rebuilt with ideal geometry from its own torsions, see
    :func:`idealize_backbone`), and every peptide bond must be trans, or
    cis before a proline (:func:`omega_state`).
2.  N, CA, C of the residue preceding the IDR are superimposed onto the
    folded domain, and psi of that residue is rotated so its carbonyl O
    matches the folded one (preserves the folded peptide plane).
3.  phi/psi angles inside the IDR are changed by clash-aware damped least
    squares (Levenberg-Marquardt inverse kinematics) until N, CA, C of the
    residue following the IDR (the endpoint) land on the folded domain.
    Each iteration solves for the smallest simultaneous change of all IDR
    phi/psi angles that moves the endpoint toward its target and pushes
    apart backbone atoms coming within a margin of each other or of the
    fixed environment. Steps are bounded per torsion and by the largest
    atom displacement, and torsions that reach a cumulative cap are frozen
    while the others are re-solved. The correction is spread thinly along
    the loop and the chain is steered around the folded domains instead of
    through them.
4.  The closed loop is accepted only if the endpoint residuals, the change
    of every phi/psi relative to the source conformer, the bond geometry
    and omega of both peptide junctions, nonlocal self/environment
    contacts and the observed (phi, psi) occupancy of the torsion database
    all pass.

Rotations act on rigid sub-trees around N-CA (phi) or CA-C (psi) bonds, so
bond lengths, bond angles, omega angles and chirality of the source
conformer are kept exactly; only phi/psi inside the IDR change, and only
by bounded amounts, so the torsion sampling of IDPConformerGenerator is
preserved.

All coordinates are in Angstrom; angle parameters are in degrees.
"""
import math
from dataclasses import dataclass, field
from statistics import mean

import numpy as np
from scipy.spatial import cKDTree

from idpconfgen.core import build_definitions as _bgeo
from idpconfgen.core.definitions import vdW_radii_tsai_1999


BB_ATOMS = ('N', 'CA', 'C', 'O')
"""Backbone atom order per residue in the CCD coordinate arrays."""

N_, CA_, C_, O_ = 0, 1, 2, 3
_PHI, _PSI = 0, 1
_BB_RADII = np.array([
    vdW_radii_tsai_1999['N'],
    vdW_radii_tsai_1999['C'],
    vdW_radii_tsai_1999['C'],
    vdW_radii_tsai_1999['O'],
    ])


IDEALIZE_MODES = ('auto', 'always', 'never')
"""How :func:`close_loop` treats the bond geometry of its input."""

GEOMETRY_WINDOWS = {
    'N_CA': (1.35, 1.56),
    'CA_C': (1.40, 1.63),
    'C_O': (1.15, 1.33),
    'C_N': (1.20, 1.47),
    'N_CA_C': (98.0, 128.0),
    'CA_C_O': (102.0, 138.0),
    'CA_C_N': (98.0, 145.0),
    'O_C_N': (98.0, 145.0),
    'C_N_CA': (102.0, 145.0),
    }
"""
Accepted backbone bond lengths (Angstrom) and angles (degrees).

``C_N``, ``CA_C_N``, ``O_C_N`` and ``C_N_CA`` refer to the peptide bond
with the next residue. These are the windows of the backbone geometry
check of the AlphaFlex loop builders.
"""

_CA_C_O = 2.0 * _bgeo.build_bend_CA_C_O  # as IDPConformerGenerator builds O
IDEAL_GEOMETRY = {
    'N_CA': mean(_bgeo.distances_N_CA.values()),
    'CA_C': mean(_bgeo.distances_CA_C.values()),
    'C_O': _bgeo.distance_C_O,
    'C_N': mean(_bgeo.distances_C_Np1.values()),
    'N_CA_C': math.degrees(mean(_bgeo.bend_angles_N_CA_C.values())),
    'CA_C_O': math.degrees(_CA_C_O),
    'CA_C_N': math.degrees(mean(_bgeo.bend_angles_CA_C_Np1.values())),
    'C_N_CA': math.degrees(mean(_bgeo.bend_angles_Cm1_N_CA.values())),
    }
"""
Ideal backbone bond lengths and angles of IDPConformerGenerator.

Averages over residue types of the ``fixed`` bond geometry (Angstrom and
degrees); the peptide group is planar, so O-C-N completes 360 degrees.
"""
IDEAL_GEOMETRY['O_C_N'] = (
    360.0 - IDEAL_GEOMETRY['CA_C_O'] - IDEAL_GEOMETRY['CA_C_N']
    )


@dataclass
class CCDParameters:
    """
    Parameters controlling torsion-space loop closure.

    Groups of fields:

    * ``idealize``: bond geometry of the source (see :func:`close_loop`);
    * ``target_*`` and ``weights``: endpoint acceptance and weights of
      N, CA and C of the residue following the IDR;
    * ``max_phipsi_*``: accepted change of phi/psi from the source;
    * ``omega_tolerance`` and ``allow_nonpro_cis``: peptide bonds;
    * ``dls_*``: the damped least-squares solver;
    * ``polish_*`` and ``max_total_rotation``: least-squares polishing of
      near-closed loops;
    * ``junction_*``: geometry of both peptide junctions;
    * ``clash_scale``, ``side_chain_radius``, ``self_*`` and
      ``max_*_clashes``: contacts (see :class:`ContactModel`);
    * ``max_rama_outlier_fraction``: database (phi, psi) occupancy.

    The phi/psi limits, the contact criteria, the peptide bond (omega)
    policy and the backbone geometry windows come from the AlphaFlex loop
    builders. The endpoint targets
    are tighter, so that the peptide junctions keep ideal geometry
    without a final minimization. Angles are in degrees and distances in
    Angstrom.
    """

    # source bond geometry: 'auto' rebuilds it ideal only if it falls
    # outside GEOMETRY_WINDOWS, 'always' always, 'never' rejects it then
    idealize: str = 'auto'

    # acceptance: residuals of N, CA, C of the residue after the IDR; the
    # N and CA residuals bound the geometry of the upper junction
    target_rmsd: float = 0.10
    target_n_residual: float = 0.05
    target_max_residual: float = 0.20
    # relative weights of N, CA and C in the closure objective
    weights: tuple = (8.0, 4.0, 2.0)

    # acceptance: change of phi/psi inside the IDR relative to the source
    max_phipsi_rms_delta: float = 25.0
    max_phipsi_delta: float = 90.0
    # polishing budget: accumulated rotation allowed per torsion
    max_total_rotation: float = 90.0

    # acceptance: every peptide bond of the IDR and both junctions must be
    # within `omega_tolerance` of trans (180) or, before a proline (or for
    # any residue if `allow_nonpro_cis`), of cis (0)
    omega_tolerance: float = 20.0
    allow_nonpro_cis: bool = False

    # damped least-squares solver
    dls_iterations: int = 600
    dls_stall: int = 100
    dls_damping: float = 0.3
    dls_prior: float = 0.05
    dls_max_step: float = 4.0
    dls_max_displacement: float = 1.0
    dls_max_move: float = 2.0
    dls_torsion_cap: float = 60.0
    dls_clash_weight: float = 1.0
    dls_self_margin: float = 0.4
    dls_env_margin: float = 0.8
    # finishing pass of near-closed loops: C-terminal window, no prior
    dls_finish_window: int = 16
    dls_finish_iterations: int = 150
    dls_finish_stall: int = 40
    dls_finish_margin: float = 0.1

    # damped least-squares endpoint polishing of near-closed loops
    polish: bool = True
    polish_trigger_rmsd: float = 3.0
    polish_iterations: int = 60
    polish_max_step: float = 10.0
    polish_damping: float = 1.0
    polish_tolerance: float = 0.05
    polish_clash_weight: float = 2.0
    polish_clash_margin: float = 0.05

    # acceptance: peptide junctions with the folded domains, deviation of
    # the C-N bond and of the CA-C-N, O-C-N and C-N-CA angles from
    # IDEAL_GEOMETRY (omega follows `omega_tolerance`)
    junction_distance_tolerance: float = 0.05
    junction_angle_tolerance: float = 10.0

    # acceptance: nonlocal contacts of the IDR (see ContactModel)
    clash_scale: float = 0.65
    # radius of the C-beta pseudo-atom that keeps room for side chains
    side_chain_radius: float = 2.4
    self_ca_cutoff: float = 2.85
    self_ignore_window: int = 2
    max_self_clashes: int = 0
    max_env_clashes: int = 0

    # acceptance: fraction of IDR residues that may move into (phi, psi)
    # bins not populated in the torsion database
    max_rama_outlier_fraction: float = 0.05


@dataclass
class CCDResult:
    """Outcome of :func:`close_loop`."""

    accepted: bool
    reason: str
    backbone: np.ndarray = None
    extra: np.ndarray = None
    rotation: np.ndarray = None
    translation: np.ndarray = None
    stats: dict = field(default_factory=dict)


def kabsch(moving, target):
    """
    Superimpose `moving` onto `target` (proper rotation only).

    Parameters
    ----------
    moving, target : np.ndarray, shape (N, 3)
        Paired points, N >= 3.

    Returns
    -------
    rotation : np.ndarray, shape (3, 3)
    translation : np.ndarray, shape (3,)
        ``moving @ rotation.T + translation`` fits `target`.
    rmsd : float
        RMSD after superposition.
    """
    moving = np.asarray(moving, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    cm = moving.mean(axis=0)
    ct = target.mean(axis=0)
    u, _s, vt = np.linalg.svd((moving - cm).T @ (target - ct))
    d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    rotation = vt.T @ np.diag((1.0, 1.0, d)) @ u.T
    translation = ct - rotation @ cm
    fitted = moving @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def dihedral(p0, p1, p2, p3):
    """Calculate the dihedral p0-p1-p2-p3 in radians (IUPAC sign)."""
    b0 = np.asarray(p0, dtype=np.float64) - p1
    b1 = np.asarray(p2, dtype=np.float64) - p1
    b2 = np.asarray(p3, dtype=np.float64) - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w))


def bend_angle(p0, p1, p2):
    """Calculate the p0-p1-p2 angle in degrees."""
    v0 = np.asarray(p0, dtype=np.float64) - p1
    v1 = np.asarray(p2, dtype=np.float64) - p1
    cos = np.dot(v0, v1) / (np.linalg.norm(v0) * np.linalg.norm(v1))
    return math.degrees(math.acos(min(1.0, max(-1.0, cos))))


def backbone_phi_psi(backbone, first, last):
    """
    Calculate phi and psi angles of residues `first` to `last`.

    Parameters
    ----------
    backbone : np.ndarray, shape (L, 4, 3)
        N, CA, C, O coordinates per residue. Residues `first` - 1 and
        `last` + 1 must exist.

    first, last : int
        Zero-based indexes of the residue range (inclusive).

    Returns
    -------
    np.ndarray, shape (last - first + 1, 2)
        Phi and psi angles in radians.
    """
    out = np.empty((last - first + 1, 2), dtype=np.float64)
    for k, i in enumerate(range(first, last + 1)):
        out[k, 0] = dihedral(
            backbone[i - 1, C_],
            backbone[i, N_],
            backbone[i, CA_],
            backbone[i, C_],
            )
        out[k, 1] = dihedral(
            backbone[i, N_],
            backbone[i, CA_],
            backbone[i, C_],
            backbone[i + 1, N_],
            )
    return out


def wrap_angle(angle):
    """Wrap angle(s) in radians to the interval [-pi, pi)."""
    return (np.asarray(angle) + np.pi) % (2 * np.pi) - np.pi


def rotation_matrices(axes, angles):
    """
    Build rotation matrices about unit vectors (Rodrigues formula).

    Parameters
    ----------
    axes : np.ndarray, shape (P, 3)
        Unit rotation axes.

    angles : np.ndarray, shape (P,)
        Right-handed rotation angles in radians.

    Returns
    -------
    np.ndarray, shape (P, 3, 3)
    """
    cos, sin = np.cos(angles), np.sin(angles)
    x, y, z = axes[:, 0], axes[:, 1], axes[:, 2]
    one = 1.0 - cos
    rot = np.empty((len(angles), 3, 3))
    rot[:, 0, 0] = cos + x * x * one
    rot[:, 0, 1] = x * y * one - z * sin
    rot[:, 0, 2] = x * z * one + y * sin
    rot[:, 1, 0] = y * x * one + z * sin
    rot[:, 1, 1] = cos + y * y * one
    rot[:, 1, 2] = y * z * one - x * sin
    rot[:, 2, 0] = z * x * one - y * sin
    rot[:, 2, 1] = z * y * one + x * sin
    rot[:, 2, 2] = cos + z * z * one
    return rot


def rotate_torsions(coords, anchors, axes, starts, angles):
    """
    Change several torsions at once by rotating coordinate suffixes.

    Torsion `k` turns ``coords[starts[k]:]`` about the axis through
    ``anchors[k]`` along ``axes[k]``. Each atom moves by the composition of
    the rotations of every torsion upstream of it, so every torsion changes
    exactly by its angle, as if the rotations were applied one by one from
    the C-terminal end.

    Parameters
    ----------
    coords : np.ndarray, shape (M, 3)
        Coordinates, modified in place.

    anchors, axes : np.ndarray, shape (P, 3)
        A point on each torsion axis and its unit direction, both taken
        from `coords` before the change.

    starts : np.ndarray of int, shape (P,)
        Index of the first atom moved by each torsion.

    angles : np.ndarray, shape (P,)
        Changes in radians.
    """
    if len(starts) == 0:
        return
    order = np.argsort(starts, kind='stable')
    first_moved = starts[order]
    rot = rotation_matrices(axes[order], angles[order])
    points = anchors[order]
    lin = np.empty((len(order), 3, 3))
    shift = np.empty((len(order), 3))
    cur_lin = np.eye(3)
    cur_shift = np.zeros(3)
    for j in range(len(order)):
        cur_shift = cur_lin @ (points[j] - rot[j] @ points[j]) + cur_shift
        cur_lin = cur_lin @ rot[j]
        lin[j] = cur_lin
        shift[j] = cur_shift
    idx = np.arange(first_moved[0], len(coords))
    seg = np.searchsorted(first_moved, idx, side='right') - 1
    coords[idx] = np.einsum('aij,aj->ai', lin[seg], coords[idx]) + shift[seg]


def torsion_displacement(coords, atoms, anchors, axes, starts, angles):
    """
    Linearised displacement of selected atoms for small torsion changes.

    Parameters
    ----------
    coords : np.ndarray, shape (M, 3)

    atoms : np.ndarray of int, shape (K,)
        Indexes of the atoms of interest.

    anchors, axes, starts, angles
        As in :func:`rotate_torsions`.

    Returns
    -------
    np.ndarray, shape (K, 3)
        ``sum_k angles[k] * axes[k] x (coords[a] - anchors[k])`` over the
        torsions moving each atom `a`.
    """
    out = np.zeros((len(atoms), 3))
    if len(starts) == 0:
        return out
    order = np.argsort(starts, kind='stable')
    spin = angles[order][:, None] * axes[order]
    spin_sum = np.cumsum(spin, axis=0)
    moment_sum = np.cumsum(np.cross(spin, anchors[order]), axis=0)
    count = np.searchsorted(starts[order], atoms, side='right')
    moved = count > 0
    k = count[moved] - 1
    out[moved] = np.cross(spin_sum[k], coords[atoms[moved]]) - moment_sum[k]
    return out


class ClashEnvironment:
    """
    Heavy atoms the loop backbone must not overlap with.

    Parameters
    ----------
    coords : np.ndarray, shape (M, 3)
        Coordinates of the fixed heavy atoms (folded domains, other
        chains, membrane...).

    radii : np.ndarray, shape (M,)
        Van der Waals radii of the fixed atoms.

    lower_mask, upper_mask : np.ndarray of bool, shape (M,), optional
        Atoms of the folded residues bonded to the first and last IDR
        residue. Contacts between those residues are covalent and are
        ignored.
    """

    def __init__(self, coords, radii, lower_mask=None, upper_mask=None):
        self.coords = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
        self.radii = np.asarray(radii, dtype=np.float64).reshape(-1)
        assert self.coords.shape[0] == self.radii.shape[0]
        size = self.coords.shape[0]
        self.lower_mask = (
            np.zeros(size, dtype=bool) if lower_mask is None
            else np.asarray(lower_mask, dtype=bool)
            )
        self.upper_mask = (
            np.zeros(size, dtype=bool) if upper_mask is None
            else np.asarray(upper_mask, dtype=bool)
            )
        self.max_radius = float(self.radii.max()) if size else 0.0
        self.tree = cKDTree(self.coords) if size else None

    def overlaps(
            self,
            coords,
            radii,
            residue,
            first,
            last,
            scale,
            offset=0.0,
            ):
        """
        Find loop-environment pairs closer than `scale` times their vdW sum.

        Parameters
        ----------
        coords : np.ndarray, shape (K, 3)
            Loop atoms.

        radii : np.ndarray, shape (K,)
            Van der Waals radii of the loop atoms.

        residue : np.ndarray of int, shape (K,)
            Residue index of each loop atom.

        first, last : int
            Residue indexes of the first and last IDR residue; their
            contacts with `lower_mask` and `upper_mask` atoms are ignored.

        scale : float
            Fraction of the vdW radii sum below which a pair overlaps.

        offset : float
            Distance added to the scaled radii sum (a tolerance).

        Returns
        -------
        tuple of np.ndarray
            Loop atom indexes, environment atom indexes and the distance
            cutoff of each overlapping pair.
        """
        loop_idx, env_idx, cut, _dist = self.pairs(
            coords,
            radii,
            residue,
            first,
            last,
            scale,
            offset,
            )
        return loop_idx, env_idx, cut + offset

    def pairs(self, coords, radii, residue, first, last, scale, margin=0.0):
        """
        Find loop-environment pairs closer than their cutoff plus `margin`.

        Same arguments as :meth:`overlaps`, with `margin` in place of
        `offset`.

        Returns
        -------
        tuple of np.ndarray
            Loop atom indexes, environment atom indexes, the cutoff
            ``scale * (r_i + r_j)`` and the distance of each pair.
        """
        empty = (np.zeros(0, int), np.zeros(0, int), np.zeros(0), np.zeros(0))
        if self.tree is None or len(coords) == 0:
            return empty
        coords = np.asarray(coords, dtype=np.float64)
        radii = np.asarray(radii, dtype=np.float64)
        residue = np.asarray(residue)
        reach = scale * (float(radii.max()) + self.max_radius) + margin
        hits = self.tree.query_ball_point(coords, reach)
        sizes = np.fromiter((len(h) for h in hits), int, len(hits))
        total = int(sizes.sum())
        if total == 0:
            return empty
        loop_idx = np.repeat(np.arange(len(coords)), sizes)
        env_idx = np.fromiter((i for h in hits for i in h), int, total)
        bonded = (
            ((residue[loop_idx] == first) & self.lower_mask[env_idx])
            | ((residue[loop_idx] == last) & self.upper_mask[env_idx])
            )
        loop_idx, env_idx = loop_idx[~bonded], env_idx[~bonded]
        dist = np.linalg.norm(coords[loop_idx] - self.coords[env_idx], axis=1)
        cut = scale * (radii[loop_idx] + self.radii[env_idx])
        near = dist < cut + margin
        return loop_idx[near], env_idx[near], cut[near], dist[near]

    def count_overlaps(self, *args, **kwargs):
        """Count pairs reported by :meth:`overlaps`."""
        return len(self.overlaps(*args, **kwargs)[0])


class ContactModel:
    """
    Atoms of the IDR used for contact checks.

    Besides the backbone heavy atoms N, CA, C and O, the model holds the
    amide H (except for Pro) and an ideal C-beta (except for Gly) of every
    IDR residue, both built from the backbone. They are rigidly attached
    to the backbone, so they follow the torsion moves of the atom they
    hang on. The C-beta is a pseudo-atom with radius `side_chain_radius`
    that keeps room for the side chain packed afterwards; the amide H
    keeps hydrogen contacts physical, as the energy function sees them.

    Pairs up to three bonds apart along the chain are never compared.

    Parameters
    ----------
    first, last : int
        Zero-based IDR residue range (inclusive) in the loop arrays.

    sequence : str
        One-letter sequence of the whole loop fragment.

    side_chain_radius : float
        Radius of the C-beta pseudo-atom.
    """

    _H_DISTANCE = 0.9  # as IDPConfGen builds the amide H
    _CB = (-0.58273431, 0.56802827, -0.54067466)

    def __init__(self, first, last, sequence, side_chain_radius=2.4):
        self.first = first
        self.last = last
        residues = np.arange(first, last + 1)
        nres = len(residues)
        with_h = np.array([sequence[r] != 'P' for r in residues], dtype=bool)
        with_cb = np.array([sequence[r] != 'G' for r in residues], dtype=bool)
        self.h_res = residues[with_h]
        self.cb_res = residues[with_cb]

        # backbone atoms, flat indexes of the loop arrays
        self.bb_index = (4 * residues[:, None] + np.arange(4)).reshape(-1)
        res = np.concatenate((
            np.repeat(residues, 4),
            self.h_res,
            self.cb_res,
            ))
        kind = np.concatenate((
            np.tile(np.arange(4), nres),
            np.full(len(self.h_res), 4),
            np.full(len(self.cb_res), 5),
            ))
        # main-chain position of the atom each atom is, or hangs on
        offset = np.array([0, 1, 2, 2, 0, 1])
        self.chain = 3 * res + offset[kind]
        self.branch = np.array([0, 0, 0, 1, 1, 1])[kind]
        # flat index of the atom whose motion each atom follows: H the N,
        # CB the C (phi moves it, psi does not)
        follow = np.array([N_, CA_, C_, O_, N_, C_])
        self.follows = 4 * res + follow[kind]
        self.residue = res
        self.is_ca = kind == CA_
        self.radius = np.concatenate((
            np.tile(_BB_RADII, nres),
            np.full(len(self.h_res), vdW_radii_tsai_1999['H']),
            np.full(len(self.cb_res), float(side_chain_radius)),
            ))

    def coordinates(self, coords):
        """
        Build the model atoms from flat loop coordinates.

        Parameters
        ----------
        coords : np.ndarray, shape (M, 3)
            Flat N, CA, C, O coordinates of the loop (residue `r` at rows
            ``4 * r`` to ``4 * r + 3``).

        Returns
        -------
        np.ndarray, shape (K, 3)
        """
        h = self.h_res
        n_at = coords[4 * h + N_]
        bisector = (
            _unit(n_at - coords[4 * (h - 1) + C_])
            + _unit(n_at - coords[4 * h + CA_])
            )
        hydrogens = n_at + self._H_DISTANCE * _unit(bisector)

        cb = self.cb_res
        n_at = coords[4 * cb + N_]
        ca_at = coords[4 * cb + CA_]
        c_at = coords[4 * cb + C_]
        b = ca_at - n_at
        c = c_at - ca_at
        a = np.cross(b, c)
        betas = (
            self._CB[0] * a + self._CB[1] * b + self._CB[2] * c + ca_at
            )
        return np.concatenate((coords[self.bb_index], hydrogens, betas))

    def self_pairs(self, coords, params, margin=0.0):
        """
        Find IDR atom pairs closer than their cutoff plus `margin`.

        A pair overlaps if closer than ``params.clash_scale`` times the
        sum of the vdW radii or, for C-alpha pairs more than
        ``params.self_ignore_window`` residues apart, than
        ``params.self_ca_cutoff``.

        Parameters
        ----------
        coords : np.ndarray, shape (K, 3)
            Output of :meth:`coordinates`.

        params : :class:`CCDParameters`

        margin : float
            Distance added to the cutoffs to also report near contacts.

        Returns
        -------
        tuple of np.ndarray
            Model atom indexes ``i`` and ``j`` (``i`` follows fewer
            torsions than ``j``), the cutoff and the distance of each pair.
        """
        empty = (np.zeros(0, int), np.zeros(0, int), np.zeros(0), np.zeros(0))
        reach = max(
            params.self_ca_cutoff,
            2 * params.clash_scale * float(self.radius.max()),
            )
        found = cKDTree(coords).query_pairs(
            r=reach + margin,
            output_type='ndarray',
            )
        if len(found) == 0:
            return empty
        i, j = found[:, 0], found[:, 1]
        bonds = np.abs(self.chain[i] - self.chain[j]) + self.branch[i] + self.branch[j]  # noqa: E501
        keep = bonds > 3
        i, j = i[keep], j[keep]
        swap = self.follows[i] > self.follows[j]
        i, j = np.where(swap, j, i), np.where(swap, i, j)
        dist = np.linalg.norm(coords[i] - coords[j], axis=1)
        cut = params.clash_scale * (self.radius[i] + self.radius[j])
        far_ca = (
            self.is_ca[i]
            & self.is_ca[j]
            & (np.abs(self.residue[i] - self.residue[j]) > params.self_ignore_window)  # noqa: E501
            )
        cut = np.where(far_ca, np.maximum(cut, params.self_ca_cutoff), cut)
        near = dist < cut + margin
        return i[near], j[near], cut[near], dist[near]

    def env_pairs(self, coords, env, params, margin=0.0):
        """
        Find IDR atom - environment pairs closer than cutoff plus `margin`.

        Returns
        -------
        tuple of np.ndarray
            As :meth:`ClashEnvironment.pairs`, indexing model atoms.
        """
        if env is None:
            return (np.zeros(0, int), np.zeros(0, int), np.zeros(0), np.zeros(0))  # noqa: E501
        return env.pairs(
            coords,
            self.radius,
            self.residue,
            self.first,
            self.last,
            params.clash_scale,
            margin,
            )


def _unit(vectors):
    return vectors / np.linalg.norm(vectors, axis=-1, keepdims=True)


def junction_geometry(backbone, lower_anchor, upper_anchor, first, last):
    """
    Measure both peptide junctions between the IDR and the folded domains.

    Parameters
    ----------
    backbone : np.ndarray, shape (L, 4, 3)
        Closed loop backbone.

    lower_anchor, upper_anchor : np.ndarray, shape (4, 3)
        Folded N, CA, C, O of the residues preceding and following the
        IDR. O may be NaN.

    first, last : int
        Zero-based indexes of the first and last IDR residues.

    Returns
    -------
    dict
        C-N distances, the three junction angles (degrees) and omega
        (degrees) of each junction. Keys start with ``lower_`` and
        ``upper_``.
    """
    out = {}
    c_lo, o_lo, ca_lo = lower_anchor[C_], lower_anchor[O_], lower_anchor[CA_]
    n_s, ca_s = backbone[first, N_], backbone[first, CA_]
    out['lower_c_n'] = float(np.linalg.norm(n_s - c_lo))
    out['lower_ca_c_n'] = bend_angle(ca_lo, c_lo, n_s)
    out['lower_o_c_n'] = (
        bend_angle(o_lo, c_lo, n_s) if np.all(np.isfinite(o_lo))
        else float('nan')
        )
    out['lower_c_n_ca'] = bend_angle(c_lo, n_s, ca_s)
    out['lower_omega'] = math.degrees(dihedral(ca_lo, c_lo, n_s, ca_s))

    ca_e, c_e, o_e = backbone[last, CA_], backbone[last, C_], backbone[last, O_]
    n_u, ca_u = upper_anchor[N_], upper_anchor[CA_]
    out['upper_c_n'] = float(np.linalg.norm(n_u - c_e))
    out['upper_ca_c_n'] = bend_angle(ca_e, c_e, n_u)
    out['upper_o_c_n'] = bend_angle(o_e, c_e, n_u)
    out['upper_c_n_ca'] = bend_angle(c_e, n_u, ca_u)
    out['upper_omega'] = math.degrees(dihedral(ca_e, c_e, n_u, ca_u))
    return out


def junction_failure(geometry, params, lower_is_pro=False, upper_is_pro=False):
    """
    Return the first junction measure outside its window, or None.

    The C-N bond and the CA-C-N, O-C-N and C-N-CA angles must be within
    ``params.junction_distance_tolerance`` and
    ``params.junction_angle_tolerance`` of :data:`IDEAL_GEOMETRY`, and
    omega must be allowed by :func:`omega_state`.

    Parameters
    ----------
    geometry : dict
        Output of :func:`junction_geometry`.

    params : :class:`CCDParameters`

    lower_is_pro, upper_is_pro : bool
        Whether the residue after each junction is a proline, which also
        allows a cis peptide bond.

    Returns
    -------
    str or None
    """
    ideal_c_n = IDEAL_GEOMETRY['C_N']
    for side, is_pro in (('lower', lower_is_pro), ('upper', upper_is_pro)):
        dist = geometry[f'{side}_c_n']
        if not abs(dist - ideal_c_n) <= params.junction_distance_tolerance:
            return f'junction_{side}_distance'
        for key in ('ca_c_n', 'o_c_n', 'c_n_ca'):
            value = geometry[f'{side}_{key}']
            if math.isnan(value):
                continue
            ideal = IDEAL_GEOMETRY[key.upper()]
            if not abs(value - ideal) <= params.junction_angle_tolerance:
                return f'junction_{side}_angle'
        state = omega_state(
            geometry[f'{side}_omega'],
            next_is_pro=is_pro,
            tolerance=params.omega_tolerance,
            allow_nonpro_cis=params.allow_nonpro_cis,
            )
        if state not in ('trans', 'cis'):
            return f'junction_{side}_omega'
    return None


def omega_state(
        omega,
        next_is_pro=False,
        tolerance=20.0,
        allow_nonpro_cis=False,
        ):
    """
    Classify a peptide bond by its omega angle.

    A peptide bond is ``'trans'`` within
    `tolerance` of 180 degrees, ``'cis'`` within `tolerance` of 0 degrees
    when the next residue is a proline (or `allow_nonpro_cis`),
    ``'nonpro_cis'`` within `tolerance` of 0 degrees otherwise, and
    ``'twisted'`` anywhere else.

    Parameters
    ----------
    omega : float
        Omega angle in degrees.

    next_is_pro : bool
        Whether the residue after the peptide bond is a proline.

    tolerance : float
        Half width of the cis and trans basins, in degrees.

    allow_nonpro_cis : bool
        Accept cis peptide bonds before any residue.

    Returns
    -------
    str
    """
    deviation = abs(float(omega)) % 360.0
    deviation = min(deviation, 360.0 - deviation)  # 0 cis, 180 trans
    if 180.0 - deviation <= tolerance:
        return 'trans'
    if deviation <= tolerance:
        return 'cis' if (next_is_pro or allow_nonpro_cis) else 'nonpro_cis'
    return 'twisted'


def peptide_omegas(backbone, start, stop):
    """
    Measure omega of the peptide bonds ``k -> k + 1``, `start` <= k < `stop`.

    Parameters
    ----------
    backbone : np.ndarray, shape (L, 4, 3)

    start, stop : int
        Zero-based residue range of the first atom of each peptide bond.

    Returns
    -------
    np.ndarray of float
        Omega angles in degrees.
    """
    k = np.arange(start, stop)
    if len(k) == 0:
        return np.zeros(0)
    return np.degrees(_dihedrals(
        backbone[k, CA_],
        backbone[k, C_],
        backbone[k + 1, N_],
        backbone[k + 1, CA_],
        ))


def omega_failures(omegas, next_is_pro, params):
    """
    Count peptide bonds that are not allowed by :func:`omega_state`.

    Parameters
    ----------
    omegas : sequence of float
        Omega angles in degrees.

    next_is_pro : sequence of bool
        Whether the residue after each peptide bond is a proline.

    params : :class:`CCDParameters`

    Returns
    -------
    dict
        Counts of ``'nonpro_cis'`` and ``'twisted'`` peptide bonds.
    """
    out = {'nonpro_cis': 0, 'twisted': 0}
    for omega, is_pro in zip(omegas, next_is_pro):
        state = omega_state(
            omega,
            next_is_pro=is_pro,
            tolerance=params.omega_tolerance,
            allow_nonpro_cis=params.allow_nonpro_cis,
            )
        if state in out:
            out[state] += 1
    return out


def _dihedrals(p0, p1, p2, p3):
    """Vectorised :func:`dihedral` for arrays of points, in radians."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    x = np.sum(v * w, axis=-1)
    y = np.sum(np.cross(b1, v) * w, axis=-1)
    return np.arctan2(y, x)


def _angles(p0, p1, p2):
    """Bend angles p0-p1-p2 for arrays of points, in degrees."""
    u = p0 - p1
    v = p2 - p1
    cos = np.sum(u * v, axis=-1) / (
        np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1)
        )
    return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))


def backbone_geometry(backbone, start, stop):
    """
    Measure backbone bond lengths and angles of residues `start` to `stop`.

    Parameters
    ----------
    backbone : np.ndarray, shape (L, 4, 3)

    start, stop : int
        Zero-based residue range (inclusive). The peptide bonds measured
        are those between consecutive residues of the range.

    Returns
    -------
    dict of np.ndarray
        Keyed as :data:`GEOMETRY_WINDOWS`; distances in Angstrom and
        angles in degrees.
    """
    bb = backbone[start:stop + 1]
    n, ca, c, o = bb[:, N_], bb[:, CA_], bb[:, C_], bb[:, O_]
    nxt_n, nxt_ca = n[1:], ca[1:]
    return {
        'N_CA': np.linalg.norm(ca - n, axis=1),
        'CA_C': np.linalg.norm(c - ca, axis=1),
        'C_O': np.linalg.norm(o - c, axis=1),
        'C_N': np.linalg.norm(nxt_n - c[:-1], axis=1),
        'N_CA_C': _angles(n, ca, c),
        'CA_C_O': _angles(ca, c, o),
        'CA_C_N': _angles(ca[:-1], c[:-1], nxt_n),
        'O_C_N': _angles(o[:-1], c[:-1], nxt_n),
        'C_N_CA': _angles(c[:-1], nxt_n, nxt_ca),
        }


def geometry_violations(backbone, start, stop, windows=None):
    """
    Count backbone bond lengths and angles outside their windows.

    Parameters
    ----------
    backbone : np.ndarray, shape (L, 4, 3)

    start, stop : int
        Zero-based residue range (inclusive), see :func:`backbone_geometry`.

    windows : dict, optional
        Defaults to :data:`GEOMETRY_WINDOWS`.

    Returns
    -------
    int
    """
    windows = GEOMETRY_WINDOWS if windows is None else windows
    count = 0
    for key, values in backbone_geometry(backbone, start, stop).items():
        low, high = windows[key]
        count += int(np.count_nonzero(~((values >= low) & (values <= high))))
    return count


def _place(a, b, c, bond, angle, torsion):
    """Place the atom bonded to `c` from internal coordinates (NeRF)."""
    bc = c - b
    bc = bc / np.linalg.norm(bc)
    normal = np.cross(b - a, bc)
    normal = normal / np.linalg.norm(normal)
    in_plane = np.cross(normal, bc)
    return c + bond * (
        -math.cos(angle) * bc
        + math.sin(angle) * (
            math.cos(torsion) * in_plane + math.sin(torsion) * normal
            )
        )


def idealize_backbone(backbone, sequence):
    """
    Rebuild a backbone with ideal bond lengths and angles.

    phi, psi and omega of `backbone` are kept; bond lengths and angles
    become those of the IDPConformerGenerator ``fixed`` bond geometry
    (per residue type), and each carbonyl O lies in its peptide plane.
    This repairs stretched bonds and malformed peptide planes (for example
    at the fragment junctions of conformers built with ``idpconfgen build
    --long``) while keeping the sampled torsions. The result is
    superimposed onto N, CA and C of residue 1.

    Parameters
    ----------
    backbone : np.ndarray, shape (L, 4, 3)
        N, CA, C, O coordinates.

    sequence : str
        One-letter sequence of the L residues.

    Returns
    -------
    np.ndarray, shape (L, 4, 3)
    """
    bb = np.asarray(backbone, dtype=np.float64).reshape(-1, 4, 3)
    nres = len(bb)

    def value(table, aa):
        return table.get(aa, mean(table.values()))

    n, ca, c, o = bb[:, N_], bb[:, CA_], bb[:, C_], bb[:, O_]
    psi = np.empty(nres)
    psi[:-1] = _dihedrals(n[:-1], ca[:-1], c[:-1], n[1:])
    # the last psi is defined by its carbonyl, opposite to the next N
    psi[-1] = dihedral(n[-1], ca[-1], c[-1], o[-1]) - math.pi
    omega = _dihedrals(ca[:-1], c[:-1], n[1:], ca[1:])
    phi = _dihedrals(c[:-1], n[1:], ca[1:], c[1:])

    out = np.zeros_like(bb)
    aa = sequence[0]
    angle = value(_bgeo.bend_angles_N_CA_C, aa)
    out[0, CA_] = (value(_bgeo.distances_N_CA, aa), 0.0, 0.0)
    out[0, C_] = out[0, CA_] + value(_bgeo.distances_CA_C, aa) * np.array(
        [-math.cos(angle), math.sin(angle), 0.0],
        )
    for i in range(1, nres):
        prev, aa = sequence[i - 1], sequence[i]
        out[i, N_] = _place(
            out[i - 1, N_], out[i - 1, CA_], out[i - 1, C_],
            value(_bgeo.distances_C_Np1, prev),
            value(_bgeo.bend_angles_CA_C_Np1, prev),
            psi[i - 1],
            )
        out[i, CA_] = _place(
            out[i - 1, CA_], out[i - 1, C_], out[i, N_],
            value(_bgeo.distances_N_CA, aa),
            value(_bgeo.bend_angles_Cm1_N_CA, aa),
            omega[i - 1],
            )
        out[i, C_] = _place(
            out[i - 1, C_], out[i, N_], out[i, CA_],
            value(_bgeo.distances_CA_C, aa),
            value(_bgeo.bend_angles_N_CA_C, aa),
            phi[i - 1],
            )
    for i in range(nres):
        out[i, O_] = _place(
            out[i, N_], out[i, CA_], out[i, C_],
            _bgeo.distance_C_O,
            _CA_C_O,
            psi[i] + math.pi,
            )

    frame = min(1, nres - 1)
    rotation, translation, _rmsd = kabsch(out[frame, :3], bb[frame, :3])
    return out @ rotation.T + translation


class RamachandranReference:
    """
    Empirical phi/psi occupancy derived from an IDPConfGen torsion database.

    Used to reject closures that push torsions into regions not observed
    in the database. Occupancy is counted in `bin_size` degree bins and
    smoothed over the 3x3 neighbouring bins; a (phi, psi) pair is an
    outlier if its smoothed bin holds fewer than `min_fraction` of the
    observations.

    Parameters
    ----------
    angles : np.ndarray, shape (N, 3)
        Omega, phi, psi angles in radians, as returned by
        :func:`idpconfgen.libs.libfilter.aligndb`. NaN rows are ignored.

    sequence : str
        One-letter residue codes aligned with `angles` (the `primary`
        output of `aligndb`); other characters are ignored.

    bin_size : int
        Bin width in degrees; must divide 360.

    min_fraction : float
        Smoothed occupancy below which a bin is an outlier.

    min_residue_counts : int
        Residue types with fewer observations use the pooled reference
        (Gly and Pro are pooled separately from the other residues).
    """

    def __init__(
            self,
            angles,
            sequence,
            bin_size=10,
            min_fraction=1e-4,
            min_residue_counts=2000,
            ):
        assert 360 % bin_size == 0
        self.nbins = 360 // bin_size
        self.bin_size = bin_size
        angles = np.asarray(angles, dtype=np.float64)
        size = min(len(angles), len(sequence))
        letters = np.frombuffer(
            sequence[:size].encode('ascii', 'replace'),
            dtype='S1',
            ).astype(str)
        phi = np.degrees(angles[:size, 1])
        psi = np.degrees(angles[:size, 2])
        valid = np.isfinite(phi) & np.isfinite(psi)

        def occupancy(mask):
            hist = np.zeros((self.nbins, self.nbins), dtype=np.float64)
            if np.any(mask):
                i = self._bin(phi[mask])
                j = self._bin(psi[mask])
                np.add.at(hist, (i, j), 1.0)
            smooth = np.zeros_like(hist)
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    smooth += np.roll(np.roll(hist, di, axis=0), dj, axis=1)
            total = hist.sum()
            return smooth >= max(1.0, min_fraction * total) if total else None

        pooled = {
            'G': valid & (letters == 'G'),
            'P': valid & (letters == 'P'),
            }
        pooled['X'] = valid & (letters != 'G') & (letters != 'P')
        self._pooled = {k: occupancy(m) for k, m in pooled.items()}
        self._allowed = {}
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            mask = valid & (letters == aa)
            if mask.sum() >= min_residue_counts:
                self._allowed[aa] = occupancy(mask)
            else:
                key = aa if aa in 'GP' else 'X'
                self._allowed[aa] = self._pooled[key]

    def _bin(self, degrees):
        return (
            ((np.asarray(degrees) + 180.0) // self.bin_size).astype(int)
            % self.nbins
            )

    def is_outlier(self, residue, phi, psi):
        """
        Check whether (phi, psi) is unobserved for `residue`.

        Parameters
        ----------
        residue : str
            One-letter residue code.

        phi, psi : float
            Angles in radians.

        Returns
        -------
        bool
        """
        allowed = self._allowed.get(residue, self._pooled['X'])
        if allowed is None:
            return False
        i = self._bin(math.degrees(phi))
        j = self._bin(math.degrees(psi))
        return not bool(allowed[i, j])


def _pivot_arrays(pivots):
    """Split pivot tuples into keys and index arrays."""
    keys = [pv[0] for pv in pivots]
    ia = np.array([pv[1] for pv in pivots], dtype=int)
    ib = np.array([pv[2] for pv in pivots], dtype=int)
    starts = np.array([pv[3] for pv in pivots], dtype=int)
    return keys, ia, ib, starts


def _damped_step(A, b, mu, prior, total, cap, resolves=3):
    """
    Solve the damped least-squares torsion step.

    Minimises ``|A x - b|^2 + (mu - prior) |x|^2 + prior |total + x|^2``
    in the dual (row) space. Torsions whose accumulated change would
    exceed `cap` are frozen and the problem re-solved for the others;
    whatever still overshoots after the last solve is clipped to `cap`.

    Parameters
    ----------
    A : np.ndarray, shape (K, P)
    b : np.ndarray, shape (K,)
    mu : float
        Total diagonal regularisation, damping squared plus `prior`.
    prior : float
        Weight pulling the accumulated change `total` back to zero.
    total : np.ndarray, shape (P,)
        Accumulated change of each torsion.
    cap : float
        Bound on the accumulated change.
    resolves : int
        Maximum number of solves.

    Returns
    -------
    np.ndarray, shape (P,)
    """
    free = np.ones(A.shape[1], dtype=bool)
    eye = np.eye(A.shape[0]) * mu
    step = np.zeros(A.shape[1])
    for _ in range(max(1, resolves)):
        Af = A * free
        v = Af.T @ b - prior * total * free
        # (Af'Af + mu I)^-1 v by the push-through identity
        step = (v - Af.T @ np.linalg.solve(Af @ Af.T + eye, Af @ v)) / mu
        after = total + step
        grow = free & (np.abs(after) > cap) & (np.sign(step) == np.sign(after))
        if not np.any(grow):
            break
        free &= ~grow
    return np.clip(total + step, -cap, cap) - total


class _LoopSearch:
    """Mutable loop-closure search state for one loop candidate."""

    def __init__(
            self,
            coords,
            nres,
            first,
            last,
            target,
            params,
            env,
            is_pro,
            sequence=None,
            ):
        self.X = coords
        self.nres = nres
        self.first = first
        self.last = last
        self.target = target
        self.params = params
        self.env = env
        self.is_pro = is_pro
        if sequence is None:
            sequence = ''.join('P' if pro else 'A' for pro in is_pro)
        self.contacts = ContactModel(
            first,
            last,
            sequence,
            side_chain_radius=params.side_chain_radius,
            )
        self.weights = np.asarray(params.weights, dtype=np.float64)
        self.weights = self.weights / self.weights.sum()
        self.end = 4 * (last + 1)
        self.max_total = math.radians(params.max_total_rotation)
        self.delta = {}
        self.iterations = 0
        self.rmsd0 = self.endpoint()[0]
        self.closed = False

    # state handling
    def snapshot(self):
        return self.X.copy(), dict(self.delta)

    def restore(self, snap):
        self.X[:] = snap[0]
        self.delta = dict(snap[1])

    def backbone(self):
        return self.X[:4 * self.nres].reshape(self.nres, 4, 3)

    # endpoint
    def endpoint(self):
        points = self.X[self.end:self.end + 3]
        res = np.sqrt(np.sum((points - self.target) ** 2, axis=1))
        rmsd = math.sqrt(float(np.mean(res * res)))
        return rmsd, float(res[0]), float(res.max()), points

    def within(self, rmsd, n_res, max_res, target):
        return rmsd <= target[0] and n_res <= target[1] and max_res <= target[2]

    @property
    def final_target(self):
        p = self.params
        return (p.target_rmsd, p.target_n_residual, p.target_max_residual)

    # physics
    def violations(self):
        """
        Count overlaps of the IDR (see :class:`ContactModel`).

        Returns
        -------
        tuple of int
            Self overlaps and overlaps with the environment.
        """
        coords = self.contacts.coordinates(self.X)
        n_self = len(self.contacts.self_pairs(coords, self.params)[0])
        n_env = len(self.contacts.env_pairs(coords, self.env, self.params)[0])
        return n_self, n_env

    def capture(self, audit):
        """Mark the loop closed if within target and free of overlaps."""
        rmsd, n_res, max_res, _ = self.endpoint()
        if audit and self.within(rmsd, n_res, max_res, self.final_target):
            n_self, n_env = self.violations()
            if (
                    n_self <= self.params.max_self_clashes
                    and n_env <= self.params.max_env_clashes
                    ):
                self.closed = True
        return rmsd, n_res, max_res

    # pivots
    def pivots(self, first, last, stride):
        """Return phi/psi pivots of residues first..last, C-terminal first."""
        out = []
        for r in range(first, last + 1, max(1, int(stride))):
            if not self.is_pro[r]:
                out.append(((r, _PHI), 4 * r + N_, 4 * r + CA_, 4 * r + C_))
            out.append(((r, _PSI), 4 * r + CA_, 4 * r + C_, 4 * r + O_))
        out.reverse()
        return out

    def _clash_rows(
            self,
            anchors,
            axes,
            starts,
            self_margin,
            env_margin,
            zone=True,
            ):
        """
        Linearise IDR contacts for the least-squares solvers.

        Parameters
        ----------
        anchors, axes, starts : np.ndarray
            Pivot axes (see :func:`rotate_torsions`).

        self_margin, env_margin : float
            Clearance beyond the cutoff that self and environment pairs
            are pushed to.

        zone : bool
            If True, pairs already within their margin are included (they
            are repelled before they overlap); otherwise only overlapping
            pairs are.

        Returns
        -------
        jac : np.ndarray, shape (K, P)
            Derivative of each pair distance with respect to each pivot.

        gap : np.ndarray, shape (K,)
            Distance each pair must grow to clear cutoff plus margin.

        overlaps : int
            Number of pairs closer than their cutoff.
        """
        model = self.contacts
        coords = model.coordinates(self.X)
        follows = model.follows
        rows, gaps = [np.zeros((0, len(starts)))], [np.zeros(0)]
        overlaps = 0

        i, j, cut, dist = model.self_pairs(
            coords,
            self.params,
            self_margin if zone else 0.0,
            )
        overlaps += int(np.count_nonzero(dist < cut))
        moving = (
            (starts[None, :] > follows[i][:, None])
            & (starts[None, :] <= follows[j][:, None])
            )
        keep = moving.any(axis=1)
        if np.any(keep):
            i, j, cut, dist = i[keep], j[keep], cut[keep], dist[keep]
            unit = (coords[j] - coords[i]) / np.maximum(dist, 1e-6)[:, None]
            # d|xj - xi|/d(theta_k) = u_k . ((xj - p_k) x unit) for the
            # pivots moving atom j but not atom i
            arm = np.cross(
                coords[j][:, None, :] - anchors[None, :, :],
                unit[:, None, :],
                )
            rows.append(np.einsum('pc,kpc->kp', axes, arm) * moving[keep])
            gaps.append(cut + self_margin - dist)

        if self.env is not None:
            k, e, cut, dist = model.env_pairs(
                coords,
                self.env,
                self.params,
                env_margin if zone else 0.0,
                )
            overlaps += int(np.count_nonzero(dist < cut))
            if len(k):
                unit = (coords[k] - self.env.coords[e]) / np.maximum(dist, 1e-6)[:, None]  # noqa: E501
                arm = np.cross(
                    coords[k][:, None, :] - anchors[None, :, :],
                    unit[:, None, :],
                    )
                row = np.einsum('pc,kpc->kp', axes, arm)
                row *= starts[None, :] <= follows[k][:, None]
                rows.append(row)
                gaps.append(cut + env_margin - dist)

        return np.vstack(rows), np.concatenate(gaps), overlaps

    def _endpoint_rows(self, anchors, axes, weights):
        """
        Linearise the weighted endpoint.

        Returns
        -------
        jac : np.ndarray, shape (9, P)
            Derivative of the endpoint N, CA, C coordinates with respect to
            each pivot, rows scaled by ``sqrt(weights)``.

        error : np.ndarray, shape (3, 3)
            Target minus current endpoint.
        """
        scale = np.sqrt(np.repeat(np.asarray(weights, dtype=np.float64), 3))
        points = self.X[self.end:self.end + 3]
        jac = np.concatenate(
            [np.cross(axes, q - anchors).T for q in points],
            axis=0,
            )
        return jac * scale[:, None], self.target - points

    def _axes(self, ia, ib):
        anchors = self.X[ia].copy()
        axes = self.X[ib] - anchors
        axes /= np.linalg.norm(axes, axis=1)[:, None]
        return anchors, axes

    def polish(self, pivots, iterations, max_step):
        """
        Close the endpoint by damped least squares in torsion space.

        Levenberg-Marquardt steps solve for the minimum-norm change of
        `pivots` that removes the weighted endpoint error and pushes apart
        overlapping backbone atoms (self and environment). The correction
        is therefore spread thinly over the loop instead of concentrated on
        a few torsions. Keeps the result if within the final target;
        otherwise restores the starting state.

        Returns
        -------
        bool
            Whether the final target was reached and kept.
        """
        if not pivots:
            return False
        p = self.params
        max_step = math.radians(max_step)
        start = self.snapshot()
        lam = p.polish_damping
        keys, ia, ib, starts = _pivot_arrays(pivots)

        def evaluate():
            anchors, axes = self._axes(ia, ib)
            jac, error = self._endpoint_rows(anchors, axes, p.weights)
            res = np.sqrt(np.repeat(p.weights, 3)) * error.reshape(-1)
            cjac, gap, _overlaps = self._clash_rows(
                anchors,
                axes,
                starts,
                p.polish_clash_margin,
                p.polish_clash_margin,
                zone=False,
                )
            jac = np.vstack((jac, p.polish_clash_weight * cjac))
            res = np.concatenate((res, p.polish_clash_weight * gap))
            return jac, res, float(res @ res), len(gap), anchors, axes

        jac, res, cost, n_overlaps, anchors, axes = evaluate()
        for _ in range(max(0, int(iterations))):
            rmsd, n_res, max_res, _ = self.endpoint()
            if (
                    n_overlaps == 0
                    and rmsd <= p.polish_tolerance
                    and self.within(rmsd, n_res, max_res, self.final_target)
                    ):
                break
            gram = jac @ jac.T + (lam * lam) * np.eye(len(res))
            step = jac.T @ np.linalg.solve(gram, res)
            used = np.array([self.delta.get(k, 0.0) for k in keys])
            room = np.maximum(self.max_total - np.abs(used), 0.0)
            limit = np.minimum(max_step, room)
            step = np.clip(step, -limit, limit)
            before = self.snapshot()
            self._turn(keys, anchors, axes, starts, step)
            new = evaluate()
            if new[2] < cost:
                jac, res, cost, n_overlaps, anchors, axes = new
                lam = max(lam * 0.5, 1e-3)
            else:
                self.restore(before)
                lam *= 4.0
                if lam > 1e4:
                    break
        rmsd, n_res, max_res, _ = self.endpoint()
        ok = self.within(rmsd, n_res, max_res, self.final_target)
        if not ok:
            self.restore(start)
            return False
        self.capture(audit=True)
        return True

    def _turn(self, keys, anchors, axes, starts, step):
        """Change the pivot torsions by `step` and book the change."""
        rotate_torsions(self.X, anchors, axes, starts, step)
        for key, theta in zip(keys, step):
            if theta:
                self.delta[key] = self.delta.get(key, 0.0) + float(theta)

    def dls(self, pivots, finish=False):
        """
        Close the loop by clash-aware damped least squares.

        Each iteration linearises the weighted endpoint error (its target
        clamped to ``dls_max_move``) and every contact within
        ``dls_self_margin``/``dls_env_margin`` of its cutoff, and solves the
        damped least-squares problem for the change of all `pivots`. The
        step is scaled so that no torsion turns more than ``dls_max_step``
        and no atom moves more than ``dls_max_displacement``; torsions that
        would exceed ``dls_torsion_cap`` in total are frozen and the others
        re-solved. ``dls_prior`` pulls the accumulated change back toward
        the source conformer where the constraints allow.

        Parameters
        ----------
        pivots : list
            Output of :meth:`pivots`.

        finish : bool
            Finishing pass for a near-closed loop: no prior, margins of
            ``dls_finish_margin`` and the ``dls_finish_*`` budget.

        Returns
        -------
        bool
            Whether the loop closed without overlaps.
        """
        if not pivots:
            return False
        p = self.params
        if finish:
            prior = 0.0
            self_margin = env_margin = p.dls_finish_margin
            iterations, stall = p.dls_finish_iterations, p.dls_finish_stall
        else:
            prior = p.dls_prior
            self_margin, env_margin = p.dls_self_margin, p.dls_env_margin
            iterations, stall = p.dls_iterations, p.dls_stall
        keys, ia, ib, starts = _pivot_arrays(pivots)
        total = np.array([self.delta.get(k, 0.0) for k in keys])
        cap = math.radians(p.dls_torsion_cap)
        max_step = math.radians(p.dls_max_step)
        mu = p.dls_damping ** 2 + prior
        tracked = np.arange(4 * self.first, self.end + 3)
        scale = np.sqrt(np.repeat(self.weights, 3))
        best = float('inf')
        no_improve = 0
        for _ in range(max(0, int(iterations))):
            self.iterations += 1
            anchors, axes = self._axes(ia, ib)
            rmsd, n_res, max_res, _ = self.endpoint()
            cjac, gap, overlaps = self._clash_rows(
                anchors,
                axes,
                starts,
                self_margin,
                env_margin,
                )
            if overlaps == 0 and self.within(rmsd, n_res, max_res, self.final_target):  # noqa: E501
                break
            merit = rmsd + 0.5 * overlaps
            if merit < best - 1e-3:
                best = merit
                no_improve = 0
            else:
                no_improve += 1
                if no_improve > stall:
                    break

            jac, error = self._endpoint_rows(anchors, axes, self.weights)
            clamp = min(1.0, p.dls_max_move / max(max_res, 1e-9))
            A = np.vstack((jac, p.dls_clash_weight * cjac))
            b = np.concatenate((
                scale * (clamp * error).reshape(-1),
                p.dls_clash_weight * gap,
                ))
            step = _damped_step(A, b, mu, prior, total, cap)

            largest = float(np.abs(step).max())
            if largest > max_step:
                step *= max_step / largest
            moved = torsion_displacement(
                self.X,
                tracked,
                anchors,
                axes,
                starts,
                step,
                )
            largest = math.sqrt(float((moved * moved).sum(axis=1).max()))
            if largest > p.dls_max_displacement:
                step *= p.dls_max_displacement / largest
            self._turn(keys, anchors, axes, starts, step)
            total += step

        self.capture(audit=True)
        return self.closed


def _try_polish(search, pivots, params, stats):
    """
    Polish a near-closed state; keep it only if closed without clashes.

    Returns
    -------
    bool
        Whether the search is closed.
    """
    if search.closed:
        return True
    if not params.polish or search.endpoint()[0] > params.polish_trigger_rmsd:
        return False
    before = search.snapshot()
    stats['polish_attempts'] = stats.get('polish_attempts', 0) + 1
    reached = search.polish(
        pivots,
        params.polish_iterations,
        params.polish_max_step,
        )
    if reached and search.closed:
        stats['polished'] = 1
        return True
    search.restore(before)
    return False


def _search_dls(search, params, stats):
    """
    Close with clash-aware damped least squares.

    A loop that stalls near its target gets a finishing pass, first on a
    C-terminal window of ``dls_finish_window`` residues (local torsions fix
    the last fraction of an Angstrom without swinging the whole chain),
    then on the whole loop, and finally the least-squares polish.
    """
    pivots = search.pivots(search.first, search.last, 1)
    closed = search.dls(pivots)
    stats['dls_rmsd'] = search.endpoint()[0]
    if not closed and stats['dls_rmsd'] <= params.polish_trigger_rmsd:
        stats['dls_finish'] = 1
        start = max(search.first, search.last - params.dls_finish_window + 1)
        closed = (
            search.dls(search.pivots(start, search.last, 1), finish=True)
            or search.dls(pivots, finish=True)
            )
    if not closed:
        _try_polish(search, pivots, params, stats)
    return None


def _final_checks(search, raw_torsions, sequence, anchors, params, rama, stats):  # noqa: E501
    """
    Validate the current state of `search`.

    Returns
    -------
    str or None
        Rejection reason, or None if the loop is accepted.
    """
    first, last = search.first, search.last
    rmsd, n_res, max_res, _ = search.endpoint()
    stats.update({
        'final_rmsd': rmsd,
        'final_n_residual': n_res,
        'final_max_residual': max_res,
        'iterations': search.iterations,
        })
    if not search.within(rmsd, n_res, max_res, search.final_target):
        return 'final_endpoint'

    closed_bb = search.backbone()
    new_torsions = backbone_phi_psi(closed_bb, first, last)
    diff = np.degrees(np.abs(wrap_angle(new_torsions - raw_torsions)))
    stats['phipsi_rms_delta'] = float(np.sqrt(np.mean(diff ** 2)))
    stats['phipsi_max_delta'] = float(diff.max())
    if stats['phipsi_rms_delta'] > params.max_phipsi_rms_delta:
        return 'phipsi_rms_delta'
    if stats['phipsi_max_delta'] > params.max_phipsi_delta:
        return 'phipsi_max_delta'

    geometry = junction_geometry(closed_bb, *anchors, first, last)
    stats.update(geometry)
    failure = junction_failure(
        geometry,
        params,
        lower_is_pro=sequence[first] == 'P',
        upper_is_pro=sequence[last + 1] == 'P',
        )
    if failure:
        return failure

    n_self, n_env = search.violations()
    stats['self_clashes'] = n_self
    stats['env_clashes'] = n_env
    if n_self > params.max_self_clashes:
        return 'self_clash'
    if n_env > params.max_env_clashes:
        return 'environment_clash'

    if rama is not None:
        # only residues moved into (phi, psi) regions the source avoided
        moved = np.any(diff > 1.0, axis=1)
        outliers = sum(
            rama.is_outlier(sequence[first + k], *new_torsions[k])
            and not rama.is_outlier(sequence[first + k], *raw_torsions[k])
            for k in np.flatnonzero(moved)
            )
        stats['rama_outliers'] = int(outliers)
        if outliers > params.max_rama_outlier_fraction * (last - first + 1):
            return 'rama_outlier'
    return None


def close_loop(
        backbone,
        sequence,
        lower_anchor,
        upper_anchor,
        environment=None,
        params=None,
        rama=None,
        extra=None,
        ):
    """
    Close a loop-IDR conformer onto two fixed anchors.

    Before closing, the source conformer is checked from the residue
    preceding the IDR to the residue following it: bond lengths and angles
    must be within :data:`GEOMETRY_WINDOWS` (with ``params.idealize``
    ``'auto'`` a conformer outside them is first rebuilt with
    :func:`idealize_backbone`; ``'always'`` rebuilds every conformer and
    ``'never'`` rejects it), and every peptide bond must be allowed by
    :func:`omega_state`, since closure keeps omega. The closed loop must
    then pass the endpoint targets, the phi/psi limits, the junction
    geometry (:func:`junction_failure`), the contact criteria and, if
    `rama` is given, the database (phi, psi) occupancy.

    Parameters
    ----------
    backbone : np.ndarray, shape (L, 4, 3)
        N, CA, C, O coordinates of the conformer. Residues 0, 1 and
        L - 2, L - 1 overlap the folded domains; residues 2 to L - 3 are
        the IDR.

    sequence : str
        One-letter sequence of the L residues.

    lower_anchor, upper_anchor : np.ndarray, shape (4, 3)
        Folded N, CA, C, O of residue 1 (preceding the IDR) and of residue
        L - 2 (following the IDR). O of `lower_anchor` aligns the first
        peptide plane; if NaN, that step is skipped.

    environment : :class:`ClashEnvironment`, optional
        Fixed heavy atoms the loop backbone must not overlap.

    params : :class:`CCDParameters`, optional

    rama : :class:`RamachandranReference`, optional
        If given, modified residues must keep observed (phi, psi) pairs.

    extra : np.ndarray, shape (K, 3), optional
        Atoms rigidly attached to the last residue (e.g., OXT); moved
        together with it.

    Returns
    -------
    :class:`CCDResult`
        `backbone` holds the closed coordinates (also for rejected
        candidates, for diagnostics; after a ``'backbone_geometry'`` or
        ``'omega'`` rejection it is the source, possibly rebuilt and not
        superimposed, and `extra`, `rotation` and `translation` are
        None); `rotation` and `translation` map atoms rigidly attached to
        residue 0 of `backbone` onto the result. `reason` is
        ``'accepted'`` or why the loop was
        rejected: ``'backbone_geometry'``, ``'omega'``,
        ``'final_endpoint'``, ``'phipsi_rms_delta'``,
        ``'phipsi_max_delta'``, ``'junction_<side>_distance'``,
        ``'junction_<side>_angle'``, ``'junction_<side>_omega'``,
        ``'self_clash'``, ``'environment_clash'`` or ``'rama_outlier'``.
    """
    params = params or CCDParameters()
    if params.idealize not in IDEALIZE_MODES:
        raise ValueError(
            f'Unknown idealize mode {params.idealize!r}; '
            f'choose from {IDEALIZE_MODES}.'
            )
    bb = np.array(backbone, dtype=np.float64).reshape(-1, 4, 3)
    nres = bb.shape[0]
    first, last = 2, nres - 3
    stats = {'loop_length': max(0, last - first + 1)}
    if last - first + 1 < 1 or len(sequence) != nres:
        return CCDResult(False, 'invalid_input', stats=stats)
    lower_anchor = np.asarray(lower_anchor, dtype=np.float64).reshape(4, 3)
    upper_anchor = np.asarray(upper_anchor, dtype=np.float64).reshape(4, 3)

    extra = (
        np.zeros((0, 3)) if extra is None
        else np.asarray(extra, dtype=np.float64).reshape(-1, 3)
        )

    def rejected(reason):
        return CCDResult(False, reason, backbone=bb, stats=stats)

    # bond geometry of the source, from the residue preceding the IDR to
    # the residue following it; closure carries the bonds to both into
    # the junctions with the folded domains
    violations = geometry_violations(bb, first - 1, last + 1)
    stats['source_geometry_violations'] = violations
    stats['idealized'] = 0
    # frame change of residue 0, for atoms attached to it
    pre_rotation, pre_translation = np.eye(3), np.zeros(3)
    if params.idealize == 'always' or (violations and params.idealize == 'auto'):  # noqa: E501
        ideal = idealize_backbone(bb, sequence)
        pre_rotation, pre_translation, _ = kabsch(bb[0, :3], ideal[0, :3])
        if len(extra):
            rot, tra, _ = kabsch(bb[-1, :3], ideal[-1, :3])
            extra = extra @ rot.T + tra
        bb = ideal
        stats['idealized'] = 1
        violations = geometry_violations(bb, first - 1, last + 1)
    if violations:
        return rejected('backbone_geometry')

    # closure never changes omega: the peptide bonds from the residue
    # before the IDR to the residue after it are checked on the source
    problems = omega_failures(
        peptide_omegas(bb, first - 1, last + 1),
        [aa == 'P' for aa in sequence[first:last + 2]],
        params,
        )
    stats['source_nonpro_cis'] = problems['nonpro_cis']
    stats['source_twisted_omega'] = problems['twisted']
    if problems['nonpro_cis'] or problems['twisted']:
        return rejected('omega')

    X = np.vstack((bb.reshape(-1, 3), extra))

    # superimpose N, CA, C of the residue preceding the IDR
    rotation, translation, rmsd = kabsch(X[4:7], lower_anchor[:3])
    X = X @ rotation.T + translation
    stats['lower_fit_rmsd'] = rmsd

    # restore the folded carbonyl plane by rotating psi of residue 1
    if np.all(np.isfinite(lower_anchor[O_])):
        a, b = X[4 + CA_], X[4 + C_]
        u = (b - a) / np.linalg.norm(b - a)
        r = X[4 + O_] - a
        s = lower_anchor[O_] - a
        r = r - (r @ u) * u
        s = s - (s @ u) * u
        theta = math.atan2(float(np.cross(r, s) @ u), float(r @ s))
        rotate_torsions(
            X,
            a[None, :],
            u[None, :],
            np.array([4 + O_]),
            np.array([theta]),
            )
        stats['lower_o_residual'] = float(
            np.linalg.norm(X[4 + O_] - lower_anchor[O_])
            )

    source = X[:4 * nres].reshape(nres, 4, 3)
    raw_torsions = backbone_phi_psi(source, first, last)
    is_pro = np.array([aa == 'P' for aa in sequence])
    search = _LoopSearch(
        X,
        nres,
        first,
        last,
        upper_anchor[:3].copy(),
        params,
        environment,
        is_pro,
        sequence=sequence,
        )
    stats['initial_rmsd'] = search.rmsd0

    reason = _search_dls(search, params, stats)

    if reason is None:
        reason = _final_checks(
            search,
            raw_torsions,
            sequence,
            (lower_anchor, upper_anchor),
            params,
            rama,
            stats,
            )
    else:
        rmsd, n_res, max_res, _ = search.endpoint()
        stats.update({
            'final_rmsd': rmsd,
            'final_n_residual': n_res,
            'final_max_residual': max_res,
            'iterations': search.iterations,
            })

    return CCDResult(
        reason is None,
        reason or 'accepted',
        backbone=search.backbone().copy(),
        extra=search.X[4 * nres:].copy(),
        rotation=rotation @ pre_rotation,
        translation=rotation @ pre_translation + translation,
        stats=stats,
        )


# side chains of closed loops

CHI_BONDS = {
    'R': (('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'NE')),
    'N': (('CA', 'CB'), ('CB', 'CG')),
    'D': (('CA', 'CB'), ('CB', 'CG')),
    'C': (('CA', 'CB'),),
    'Q': (('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD')),
    'E': (('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD')),
    'H': (('CA', 'CB'), ('CB', 'CG')),
    'I': (('CA', 'CB'), ('CB', 'CG1')),
    'L': (('CA', 'CB'), ('CB', 'CG')),
    'K': (('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'CE')),
    'M': (('CA', 'CB'), ('CB', 'CG'), ('CG', 'SD')),
    'F': (('CA', 'CB'), ('CB', 'CG')),
    'S': (('CA', 'CB'),),
    'T': (('CA', 'CB'),),
    'W': (('CA', 'CB'), ('CB', 'CG')),
    'Y': (('CA', 'CB'), ('CB', 'CG')),
    'V': (('CA', 'CB'),),
    }
"""Bonds rotated by chi1, chi2, ... of each residue type (Pro has none)."""

_MAIN_CHAIN = ('N', 'CA', 'C', 'O', 'H')


def _element(name):
    return str(name).strip().lstrip('0123456789')[:1].upper()


def covalent_bonds(coords, names, residues):
    """
    Find the covalent bonds of a protein fragment.

    Bonds are taken within residues from distances (1.95 Angstrom, 2.1
    with S, 1.2 with H) and between C of a residue and N of the next one,
    so that overlapping atoms of different residues are never bonded.

    Parameters
    ----------
    coords : np.ndarray, shape (K, 3)

    names : sequence of str
        Atom names.

    residues : np.ndarray of int, shape (K,)
        Residue number of each atom.

    Returns
    -------
    np.ndarray of int, shape (B, 2)
    """
    coords = np.asarray(coords, dtype=np.float64)
    names = np.asarray(names).astype(str)
    residues = np.asarray(residues)
    pairs = cKDTree(coords).query_pairs(2.1, output_type='ndarray')
    out = np.zeros((0, 2), dtype=int)
    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        elements = np.array([_element(n) for n in names])
        dist = np.linalg.norm(coords[i] - coords[j], axis=1)
        with_s = (elements[i] == 'S') | (elements[j] == 'S')
        with_h = (elements[i] == 'H') | (elements[j] == 'H')
        limit = np.where(with_s, 2.1, np.where(with_h, 1.2, 1.95))
        keep = (residues[i] == residues[j]) & (dist < limit)
        out = pairs[keep]
    carbons = np.flatnonzero(names == 'C')
    nitrogens = {
        int(residues[k]): int(k) for k in np.flatnonzero(names == 'N')
        }
    peptide = [
        (int(k), nitrogens[int(residues[k]) + 1])
        for k in carbons
        if int(residues[k]) + 1 in nitrogens
        ]
    if peptide:
        out = np.vstack((out, np.array(peptide, dtype=int)))
    return out


def chi_pivots(names, residues, sequence, movable, bonds, offset=1):
    """
    List the chi rotations of the side chains of `movable` residues.

    Parameters
    ----------
    names : sequence of str
        Atom names.

    residues : np.ndarray of int
        Residue number of each atom.

    sequence : str
        One-letter sequence; residue number ``r`` is ``sequence[r -
        offset]``.

    movable : iterable of int
        Residue numbers whose chi angles may rotate.

    bonds : np.ndarray of int, shape (B, 2)
        Output of :func:`covalent_bonds`.

    offset : int
        Residue number of the first residue of `sequence`.

    Returns
    -------
    list of tuple
        ``(residue, axis_start, axis_end, moving_atoms)`` per chi, with
        atom indexes; the moving atoms are those beyond the rotated bond.
    """
    names = np.asarray(names).astype(str)
    residues = np.asarray(residues)
    neighbours = {}
    for i, j in bonds:
        neighbours.setdefault(int(i), set()).add(int(j))
        neighbours.setdefault(int(j), set()).add(int(i))
    out = []
    for residue in movable:
        residue = int(residue)
        aa = sequence[residue - offset].upper()
        index = {
            names[k]: int(k) for k in np.flatnonzero(residues == residue)
            }
        side = {k for name, k in index.items() if name not in _MAIN_CHAIN}
        for start_name, end_name in CHI_BONDS.get(aa, ()):
            if start_name not in index or end_name not in index:
                break
            start, end = index[start_name], index[end_name]
            moving, stack, seen = set(), [end], {start, end}
            while stack:
                atom = stack.pop()
                for other in neighbours.get(atom, ()):
                    if other not in seen and other in side:
                        seen.add(other)
                        moving.add(other)
                        stack.append(other)
            if moving:
                out.append((residue, start, end, np.array(sorted(moving))))
    return out


def relax_side_chains(
        coords,
        names,
        residues,
        sequence,
        movable,
        environment=None,
        env_ignore=None,
        offset=1,
        hard_scale=0.70,
        target_scale=0.78,
        iterations=60,
        max_step=10.0,
        max_total=60.0,
        damping=0.1,
        ):
    """
    Relieve overlaps of packed side chains by rotating their chi angles.

    Side-chain packers (FASPR, MCSCE) choose rotamers from a library and
    may leave atoms of a closed loop overlapping each other or the folded
    domains. Damped least-squares steps on the chi angles of the `movable`
    residues push every pair closer than `target_scale` times its vdW
    radii sum apart, until none is closer than `hard_scale` times it.
    Backbone, bond lengths and bond angles do not change, so closure and
    backbone geometry are untouched; it takes the place of a restrained
    minimization after packing.

    Parameters
    ----------
    coords : np.ndarray, shape (K, 3)
        Atoms of the loop fragment (heavy atoms and any hydrogens);
        rows with NaN are ignored.

    names : sequence of str
        Atom names.

    residues : np.ndarray of int, shape (K,)
        Residue number of each atom.

    sequence : str
        One-letter sequence, see :func:`chi_pivots`.

    movable : iterable of int
        Residue numbers whose chi angles may rotate.

    environment : :class:`ClashEnvironment`, optional
        Fixed atoms the side chains must not overlap.

    env_ignore : np.ndarray of bool, shape (K,), optional
        Atoms whose contacts with `environment` are ignored.

    offset : int
        Residue number of the first residue of `sequence`.

    hard_scale, target_scale : float
        Fractions of the vdW radii sum: pairs below `hard_scale` are
        overlaps, pairs below `target_scale` are pushed apart.

    iterations : int

    max_step, max_total : float
        Largest change of a chi angle per iteration and in total
        (degrees).

    damping : float
        Damping of the least-squares steps.

    Returns
    -------
    np.ndarray, shape (K, 3)
        The relaxed coordinates.

    dict
        ``overlaps_before``, ``overlaps_after``, ``iterations`` and
        ``max_chi_change`` (degrees).
    """
    X = np.array(coords, dtype=np.float64)
    names = np.asarray(names).astype(str)
    residues = np.asarray(residues)
    valid = np.all(np.isfinite(X), axis=1)
    index = np.flatnonzero(valid)
    sub = X[index]
    sub_names = names[index]
    sub_res = residues[index]
    radii = np.array([
        vdW_radii_tsai_1999.get(_element(n), 1.7) for n in sub_names
        ])
    bonds = covalent_bonds(sub, sub_names, sub_res)
    pivots = chi_pivots(sub_names, sub_res, sequence, movable, bonds, offset)
    stats = {'overlaps_before': 0, 'overlaps_after': 0, 'iterations': 0,
             'max_chi_change': 0.0}
    if not pivots:
        return X, stats

    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra
    size = len(sub)
    graph = csr_matrix(
        (np.ones(2 * len(bonds)),
         (np.r_[bonds[:, 0], bonds[:, 1]], np.r_[bonds[:, 1], bonds[:, 0]])),
        shape=(size, size),
        )
    near_bonds = dijkstra(graph, unweighted=True, limit=3.5) <= 3

    moves = np.zeros((len(pivots), size), dtype=bool)
    for p, (_res, _a, _b, moving) in enumerate(pivots):
        moves[p, moving] = True
    movable_atoms = moves.any(axis=0)
    ignore = (
        np.zeros(size, dtype=bool) if env_ignore is None
        else np.asarray(env_ignore, dtype=bool)[index]
        )
    env_ok = movable_atoms & ~ignore
    reach = target_scale * 2.0 * max(
        float(radii.max()),
        environment.max_radius if environment is not None else 0.0,
        )
    total = np.zeros(len(pivots))
    cap = math.radians(max_total)
    step_cap = math.radians(max_step)
    # leaf to root: later chis of a residue first, so axes stay valid
    order = sorted(range(len(pivots)), key=lambda p: (pivots[p][0], -p))

    def contacts():
        found = cKDTree(sub).query_pairs(reach, output_type='ndarray')
        i, j = found[:, 0], found[:, 1]
        keep = ~near_bonds[i, j] & (movable_atoms[i] | movable_atoms[j])
        i, j = i[keep], j[keep]
        dist = np.linalg.norm(sub[i] - sub[j], axis=1)
        vdw = radii[i] + radii[j]
        zone = dist < target_scale * vdw
        i, j, dist, vdw = i[zone], j[zone], dist[zone], vdw[zone]
        k = e = np.zeros(0, dtype=int)
        edist = evdw = np.zeros(0)
        if environment is not None and environment.tree is not None:
            atoms = np.flatnonzero(env_ok)
            hits = environment.tree.query_ball_point(sub[atoms], reach)
            sizes = np.fromiter((len(h) for h in hits), int, len(hits))
            k = np.repeat(atoms, sizes)
            e = np.fromiter((x for h in hits for x in h), int, int(sizes.sum()))
            edist = np.linalg.norm(sub[k] - environment.coords[e], axis=1)
            evdw = radii[k] + environment.radii[e]
            zone = edist < target_scale * evdw
            k, e, edist, evdw = k[zone], e[zone], edist[zone], evdw[zone]
        hard = int(np.count_nonzero(dist < hard_scale * vdw))
        hard += int(np.count_nonzero(edist < hard_scale * evdw))
        return (i, j, dist, vdw), (k, e, edist, evdw), hard

    selfs, envs, hard = contacts()
    stats['overlaps_before'] = hard
    for it in range(max(0, int(iterations))):
        if hard == 0:
            break
        stats['iterations'] = it + 1
        anchors = np.array([sub[b] for _r, _a, b, _m in pivots])
        axes = np.array([sub[b] - sub[a] for _r, a, b, _m in pivots])
        axes /= np.linalg.norm(axes, axis=1)[:, None]
        rows, gaps = [], []
        i, j, dist, vdw = selfs
        if len(i):
            unit = (sub[j] - sub[i]) / np.maximum(dist, 1e-6)[:, None]
            mi, mj = moves[:, i].T, moves[:, j].T
            arm_j = np.cross(axes[None], sub[j][:, None] - anchors[None])
            arm_i = np.cross(axes[None], sub[i][:, None] - anchors[None])
            rows.append(
                np.einsum('kpc,kc->kp', arm_j, unit) * (mj & ~mi)
                - np.einsum('kpc,kc->kp', arm_i, unit) * (mi & ~mj)
                )
            gaps.append(target_scale * vdw - dist)
        k, e, edist, evdw = envs
        if len(k):
            unit = (sub[k] - environment.coords[e]) / np.maximum(edist, 1e-6)[:, None]  # noqa: E501
            arm = np.cross(axes[None], sub[k][:, None] - anchors[None])
            rows.append(np.einsum('kpc,kc->kp', arm, unit) * moves[:, k].T)
            gaps.append(target_scale * evdw - edist)
        A = np.vstack(rows)
        b = np.concatenate(gaps)
        step = A.T @ np.linalg.solve(A @ A.T + damping * np.eye(len(b)), b)
        step = np.clip(step, -step_cap, step_cap)
        step = np.clip(total + step, -cap, cap) - total
        for p in order:
            theta = float(step[p])
            if abs(theta) < 1e-9:
                continue
            moving = pivots[p][3]
            rot = rotation_matrices(axes[p][None], np.array([theta]))[0]
            sub[moving] = (sub[moving] - anchors[p]) @ rot.T + anchors[p]
        total += step
        selfs, envs, hard = contacts()

    stats['overlaps_after'] = hard
    stats['max_chi_change'] = float(np.degrees(np.abs(total).max()))
    X[index] = sub
    return X, stats
