"""Mathemical calculations."""
import math
from math import cos, sin

import numpy as np
from numba import njit

from idpconfgen.core.build_definitions import (
    build_bend_CA_C_OXT,
    distance_C_OXT,
    )


# NOTE:
# numba.njit functions are defined at the end of the module

AXIS_111 = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    ])


def make_axis_vectors(A, B, C):
    """
    Make axis vectors.

    For np.array([x,y,z]) of connected atoms A->B->C.

    Example
    -------
        >>> av = np.array([0.000, 0.000, 0.000])
        >>> bv = np.array([1.458, 0.000, 0.000])
        >>> cv = np.array([2.009, 1.420, 0.000])
        >>> make_axis_vectors(av, bv, cv)
        >>> (array([ 0.,  0., -1.]),
        >>>  array([-1.,  0.,  0.]),
        >>>  array([-0., -1., -0.]))

    Parameters
    ----------
    A, B, C : np.array of shape (3,).
        3D coordinate points in space.

    Returns
    -------
    tuple
        Three vectors defining the of ABC plane.
    """
    AB_vect = np.subtract(A, B)

    # vector parallel to ABC plane
    parallel_ABC = np.cross(
        AB_vect,
        np.subtract(C, B),
        )

    # perpendicular to parallel
    perpendicular = np.cross(
        AB_vect,
        parallel_ABC,
        )

    Alen = np.linalg.norm(AB_vect)
    if Alen > 0.0:
        AB_vect /= Alen

    Nlen = np.linalg.norm(parallel_ABC)
    if Nlen > 0.0:
        parallel_ABC /= Nlen

    Clen = np.linalg.norm(perpendicular)
    if Clen > 0.0:
        perpendicular /= Clen

    return parallel_ABC, AB_vect, perpendicular


def rotation_to_plane(A, B, C):
    """
    Define rotations matricex for planar orientation.

    Where, A->B determine the X axis, Y is the normal vector to the ABC
        plane, and A is defined at 0,0,0.

    .. seealso::

        :func:`make_axis_vectors`

    Parameters
    ----------
    A, B, C: np.array of shape (3,)
        The 3D coordinates, where A is the parent coordinate,
        B is the Xaxis and C the Yaxis.

    Returns
    -------
    np.array of shape (3,3)
        Rotational matrix
    """
    parallel_ABC, AB_vect, perpendicular = make_axis_vectors(A, B, C)
    b = np.array([AB_vect, -perpendicular, parallel_ABC])

    v = AXIS_111

    return np.dot(b, v).T


def make_coord(theta, phi, distance, parent, xaxis, yaxis):
    """
    Make a new coordinate in space.

    .. seealso::

        :func:`make_coord_from_angles`, :func:`rotation_to_plane`.

    Parameters
    ----------
    theta : float
        The angle in radians between `parent` and `yaxis`.

    phi : float
        The torsion angles in radians between `parent-xaxis-yaxis`
        plane and new coordinate.

    distance : float
        The distance between `parent` and the new coordinate.

    parent : np.array of shape (3,)
        The coordinate in space of the parent atom (point). The parent
        atom is the one that preceeds the newly added coordinate.

    xaxis : np.array of shape (3,)
        The coordinate in space that defines the xaxis of the
        parent-xaxis-yaxis coordinates plane.

    yaxis : np.array of shape (3,)
        The coordinate in space that defines the yaxis of the
        parent-xaxis-yaxis coordinates plane.

    Returns
    -------
    np.array of shape (3,)
        The new coordinate in space.
    """
    new_coord = make_coord_from_angles(theta, phi, distance)
    rotation = rotation_to_plane(parent, xaxis, yaxis)
    new_coord = np.dot(new_coord, rotation.T)
    return new_coord + parent


def make_coord_from_angles(theta, phi, distance):
    """
    Make axis components from angles.

    Performs:
        np.array([
            distance * math.cos(phi),
            distance * math.sin(phi) * math.cos(theta),
            distance * math.sin(phi) * math.sin(theta),
            ])

    Returns
    -------
    np.array of shape (3,)
    """
    return np.array([
        distance * math.cos(phi),
        distance * math.sin(phi) * math.cos(theta),
        distance * math.sin(phi) * math.sin(theta),
        ])


def calc_torsion_angles(
        coords,
        ARCTAN2=np.arctan2,
        CROSS=np.cross,
        DIAGONAL=np.diagonal,
        MATMUL=np.matmul,
        NORM=np.linalg.norm,
        ):
    """
    Calculate torsion angles from sequential coordinates.

    Uses ``NumPy`` to compute angles in a vectorized fashion.
    Sign of the torsion angle is also calculated.

    Uses Prof. Azevedo implementation:
    https://azevedolab.net/resources/dihedral_angle.pdf

    Example
    -------
    Given the sequential coords that represent a dummy molecule of
    four atoms:

    >>> xyz = numpy.array([
    >>>     [0.06360, -0.79573, 1.21644],
    >>>     [-0.47370, -0.10913, 0.77737],
    >>>     [-1.75288, -0.51877, 1.33236],
    >>>     [-2.29018, 0.16783, 0.89329],
    >>>     ])

    A1---A2
           \
            \
            A3---A4

    Calculates the torsion angle in A2-A3 that would place A4 in respect
    to the plane (A1, A2, A3).

    Likewise, for a chain of N atoms A1, ..., An, calculates the torsion
    angles in (A2, A3) to (An-2, An-1). (A1, A2) and (An-1, An) do not
    have torsion angles.

    If coords represent a protein backbone consisting of N, CA, and C
    atoms and starting at the N-terminal, the torsion angles are given
    by the following slices to the resulting array:

    - phi (N-CA), [2::3]
    - psi (CA-C), [::3]
    - omega (C-N), [1::3]

    Parameters
    ----------
    coords : numpy.ndarray of shape (N>=4, 3)
        Where `N` is the number of atoms, must be equal or above 4.

    Returns
    -------
    numpy.ndarray of shape (N - 3,)
        The torsion angles in radians.
        If you want to convert those to degrees just apply
        ``np.degrees`` to the returned result.
    """
    # requires
    assert coords.shape[0] > 3
    assert coords.shape[1] == 3

    crds = coords.T

    # Yes, I always write explicit array indices! :-)
    q_vecs = crds[:, 1:] - crds[:, :-1]
    cross = CROSS(q_vecs[:, :-1], q_vecs[:, 1:], axis=0)
    unitary = cross / NORM(cross, axis=0)

    # components
    # u0 comes handy to define because it fits u1
    u0 = unitary[:, :-1]

    # u1 is the unitary cross products of the second plane
    # that is the unitary q2xq3, obviously applied to the whole chain
    u1 = unitary[:, 1:]

    # u3 is the unitary of the bonds that have a torsion representation,
    # those are all but the first and the last
    u3 = q_vecs[:, 1:-1] / NORM(q_vecs[:, 1:-1], axis=0)

    # u2
    # there is no need to further select dimensions for u2, those have
    # been already sliced in u1 and u3.
    u2 = CROSS(u3, u1, axis=0)

    # calculating cos and sin of the torsion angle
    # here we need to use the .T and np.diagonal trick to achieve
    # broadcasting along the whole coords chain
    # np.matmul is preferred to np.dot in this case
    # https://numpy.org/doc/stable/reference/generated/numpy.matmul.html
    cos_theta = DIAGONAL(MATMUL(u0.T, u1))
    sin_theta = DIAGONAL(MATMUL(u0.T, u2))

    # torsion angles
    return -ARCTAN2(sin_theta, cos_theta)



def calc_MSMV(data):
    """Calculate Mean, STD, Median, and Variance."""
    return np.mean(data), np.std(data), np.median(data), np.var(data)


# @njit
# def xxhamiltonian_multiplication_Q(a1, b1, c1, d1, a2, b2, c2, d2):
#     """Hamiltonian Multiplication."""
#     return (
#             (a1 * a2.T) - (b1 * b2.T) - (c1 * c2.T) - (d1 * d2.T),
#             (a1 * b2.T) + (b1 * a2.T) + (c1 * d2.T) - (d1 * c2.T),
#             (a1 * c2.T) - (b1 * d2.T) + (c1 * a2.T) + (d1 * b2.T),
#             (a1 * d2.T) + (b1 * c2.T) - (c1 * b2.T) + (d1 * a2.T),
#     )


@njit
def hamiltonian_multiplication_Q(a1, b1, c1, d1, a2, b2, c2, d2):
    """Hamiltonian Multiplication."""
    return (
        (a1 * a2) - (b1 * b2) - (c1 * c2) - (d1 * d2),
        (a1 * b2) + (b1 * a2) + (c1 * d2) - (d1 * c2),
        (a1 * c2) - (b1 * d2) + (c1 * a2) + (d1 * b2),
        (a1 * d2) + (b1 * c2) - (c1 * b2) + (d1 * a2),
        )


def calc_angle(
        v1,
        v2,
        ARCCOS=np.arccos,
        CLIP=np.clip,
        DOT=np.dot,
        NORM=np.linalg.norm,
        ):
    """Calculate the angle between two vectors."""
    # https://stackoverflow.com/questions/2827393/
    # assert v1.shape == v2.shape, (v1.shape, v2.shape)
    # assert v1.dtype == v2.dtype
    v1_u = v1 / NORM(v1)
    v2_u = v2 / NORM(v2)
    return ARCCOS(CLIP(DOT(v1_u, v2_u), -1.0, 1.0))


@njit
def calc_angle_njit(
        v1,
        v2,
        ARCCOS=np.arccos,
        DOT=np.dot,
        NORM=np.linalg.norm,
        ):
    """Calculate the angle between two vectors."""
    # https://stackoverflow.com/questions/2827393/
    v1_u = v1 / NORM(v1)
    v2_u = v2 / NORM(v2)

    dot_ncan = np.dot(v1_u, v2_u)

    if dot_ncan < -1.0:
        dot_ncan_clean = -1.0

    elif dot_ncan > 1.0:
        dot_ncan_clean = 1.0

    else:
        dot_ncan_clean = dot_ncan

    return ARCCOS(dot_ncan_clean)


@njit
def place_sidechain_template(
        bb_cnf,
        ss_template,
        CROSS=np.cross,
        NORM=np.linalg.norm,
        ):
    """
    Place sidechain templates on backbone.

    Sidechain residue template is expected to have CA already at 0,0,0.

    Parameters
    ----------
    bb_cnf : numpy nd.array, shape (3, 3), dtype=float64
        The backbone coords in the form of: N-CA-C
        Coordinates are not expected to be at any particular position.

    ss_template : numpy nd.array, shape (M, 3), dtype=float64
        The sidechain all-atom template. **Expected** to have the CA atom
        at the origin (0, 0, 0). This requirement could be easily
        removed but it is maintained for performance reasons and
        considering in the context where this function is meant
        to be used.

    Returns
    -------
    nd.array, shape (M, 3), dtype=float64
        The displaced side chain coords. All atoms are returned.
    """
    # places bb with CA at 0,0,0
    bbtmp = np.full(bb_cnf.shape, np.nan)
    bbtmp[:, :] = bb_cnf[:, :] - bb_cnf[1, :]

    # the sidechain residue template is expected to have CA
    # already at the the origin (0,0,0)
    N_CA = bbtmp[0, :]
    N_CA_ = ss_template[0, :]

    N_CA_N = calc_angle_njit(N_CA, N_CA_)

    # rotation vector
    rv = CROSS(N_CA_, N_CA)
    rvu = rv / NORM(rv)

    # aligns the N-CA vectors
    rot1 = rotate_coordinates_Q_njit(ss_template, rvu, N_CA_N)

    # starts the second rotation to align the CA-C vectors
    # calculates the cross vectors of the planes N-CA-C
    cross_cnf = CROSS(bbtmp[0, :], bbtmp[2, :])
    cross_ss = CROSS(rot1[0, :], rot1[2, :])

    # the angle of rotation is the angle between the plane normal
    angle = calc_angle_njit(cross_ss, cross_cnf)

    # plane rotation vector is the cross vector between the two plane normals
    rv = CROSS(cross_ss, cross_cnf)
    rvu = rv / NORM(rv)

    # aligns to the CA-C vector maintaining the N-CA in place
    rot2 = rotate_coordinates_Q_njit(rot1, rvu, angle)

    return rot2[:, :] + bb_cnf[1, :]


def rotate_coordinates_Q(
        coords,
        rot_vec,
        angle_rad,
        ARRAY=np.array,
        HMQ=hamiltonian_multiplication_Q,
        VSTACK=np.vstack,
        ):
    """
    Rotate coordinates by radians along an axis.

    Rotates using quaternion operations.

    Parameters
    ----------
    coords : nd.array (N, 3), dtype=np.float64
        The coordinates to rotate.

    rot_vec : (,3)
        A 3D space vector around which to rotate coords.
        Rotation vector **must** be a unitary vector.

    angle_rad : float
        The angle in radians to rotate the coords.

    Returns
    -------
    nd.array shape (N, 3), dtype=np.float64
        The rotated coordinates
    """
    # assert coords.shape[1] == 3

    b2, b3, b4 = sin(angle_rad / 2) * rot_vec
    b1 = cos(angle_rad / 2)

    c1, c2, c3, c4 = HMQ(
        b1, b2, b3, b4,
        0, coords[:, 0], coords[:, 1], coords[:, 2],
        )

    _, d2, d3, d4 = HMQ(
        c1, c2, c3, c4,
        b1, -b2, -b3, -b4,
        )

    rotated = VSTACK((d2, d3, d4)).T

    assert rotated.shape[1] == 3
    return rotated


@njit
def make_coord_Q(
        v1,
        v2,
        v3,
        distance,
        bend,
        torsion,
        ARRAY=np.array,
        CROSS=np.cross,
        NORM=np.linalg.norm,
        QM=hamiltonian_multiplication_Q,
        SIN=sin,
        COS=cos,
        ):
    """
    Make a new coords from 3 vectors using Quaternions.

    Angles must be given in radians. For a matter of performance,
    angles will not be converted to radians if given in degrees.

    Parameters
    ----------
    v1, v2, v3 : np.ndarray, shape (3,), dtype=np.float
        The vectors that define the three points required to place
        the new coordinate according to bend and torsion angles.

    distance : float
        The distance between `v3` and the new coordinate.

    bend : float
        The angle in radians for the bend angle between the three atoms.
        The actual `bend` value input in this function must be
        `(pi - bend) / 2`, this calculation must be computed outside
        this function for perfomance reasons.

    torsion : float
        The torsion angle (radians) around v2-v3 which will place
        the new coordinate correctly.
        Contrarily to the `bend` angle, do not compute any additional
        calculations and just provide the torsion angle value as is.

    Returns
    -------
    np.ndarray of shape (3,), dtype=np.float32
    """
    # transfer vectors to origin
    o1 = v1 - v2
    o2 = v3 - v2

    ocross = CROSS(o2, o1)  # changed
    u_ocross = ocross / NORM(ocross)

    # creates quaterion to rotate on the bend angle
    b2, b3, b4 = SIN(bend) * u_ocross
    b1 = COS(bend)

    # rotates the unitary of o2 according to bend angle
    uo2 = o2 / NORM(o2)
    p2, p3, p4 = uo2  # p1 is zero according to Quaternion theory
    n1, n2, n3, n4 = QM(
        *QM(b1, b2, b3, b4, 0, p2, p3, p4),
        b1, -b2, -b3, -b4,
        )

    # rotates the previous result according to torsion angle
    torsion_angle = torsion / 2
    t2, t3, t4 = SIN(torsion_angle) * uo2
    t1 = COS(torsion_angle)

    f1, f2, f3, f4 = QM(
        *QM(t1, t2, t3, t4, 0, n2, n3, n4),
        t1, -t2, -t3, -t4,
        )

    # the new rotated vector is unitary
    # extends to the correct lenth
    fov_size = ARRAY((f2 * distance, f3 * distance, f4 * distance))
    # the above implementation takes 857 ns, and is faster than
    # np.array([f2, f3, f4]) * distance which takes 1.8 us
    # given by %%timeit jupyter notebook

    # transfers the new coord created in origin
    # to the correct space position
    # performing this sum at the moment of fov_size is not faster
    # than doing it separately here
    return fov_size + v3


@njit
def make_coord_Q_planar(
        vector1,
        center_point,
        vector2,
        distance,
        bend,
        ARRAY=np.array,
        CROSS=np.cross,
        NORM=np.linalg.norm,
        QM=hamiltonian_multiplication_Q,
        SIN=sin,
        COS=cos,
        ):
    """
    Create carbonyl (C=O) coordinate for protein backbone.

    Uses rotation by quaternions logic.

    Parameters
    ----------
    vector1, center_point, vector2 : np.ndarray, shpe(3,), dtype=np.float
        The coordinates of the three point that form a plane.
        `center_point` is considered the origin of vector1 and vector2.

    distance : float
        The distance between center_point and the new point.

    bend : float
        The angle in radians for vector1 - center_point - new system.
        `bend` should be half of the desired value.

    Returns
    -------
    np.ndarray of shape (3,), dtype=np.float32
    """
    o1 = vector1 - center_point
    o2 = vector2 - center_point

    ocross = CROSS(o2, o1)
    u_ocross = ocross / NORM(ocross)

    # creates quaterion to rotate on the bend angle
    b2, b3, b4 = SIN(bend) * u_ocross
    b1 = COS(bend)

    # rotates a the unitary of o2 according to bend angle
    uo1 = o1 / NORM(o1)
    p2, p3, p4 = uo1  # p1 is zero according to Quaternion theory
    n1, n2, n3, n4 = QM(
        *QM(b1, b2, b3, b4, 0, p2, p3, p4),
        b1, -b2, -b3, -b4,
        )

    return ARRAY((n2 * distance, n3 * distance, n4 * distance)) + center_point


@njit
def make_coord_Q_COO(
        CA_term,
        C_term,
        distance=distance_C_OXT,
        bend=build_bend_CA_C_OXT,
        ARRAY=np.array,
        CROSS=np.cross,
        NORM=np.linalg.norm,
        QM=hamiltonian_multiplication_Q,
        SIN=sin,
        COS=cos,
        ):
    """
    Create the C-terminal carboxyl coordinates.

    Uses a strategy based on Quaternion rotations.

    Parameters
    ----------
    CA_term : np.ndarray, shape (3,), dtype=np.float
        The coordinates of the terminal CA atom.

    C_term : np.ndarray, shape (3,), dtype=np.float
        The coordinates of the terminal C atom.

    distance : float, optional
        The distance of the C-O and C-OXT bond lengths.
        The distance is considered the same for both atom pairs.
        Defaults to 1.27.

    bend : float, optional
        The angle in radians for the CA-C-O and CA-C-OXT bonds.
        The angle between both cases is considered the same.
        Defaults to 2 * pi / 3 (120º).
        If `bend` is given, consider the actual `bend` value must be
        `(pi - bend) / 2`, this calculation must be computed outside
        this function for perfomance reasons.

    Returns
    -------
    tuple of length 2
        np.ndarray of shape (3,), dtype=np.float32
    """
    o1 = C_term - CA_term
    # creates an inmaginary coordinate perpendicular to o1 and origin
    o2 = o1[::-1]

    # creates the cross vector that will be used to rotate the new coordinates
    ocross2 = CROSS(o1, o2)
    u_ocross2 = ocross2 / NORM(ocross2)

    # creates quaterions to rotate O and OXT on the bend angle
    # O and OXT have the same rotation angle along the same axis but
    # with opposite angle signs
    minus_bend = -bend
    b2, b3, b4 = SIN(minus_bend) * u_ocross2
    b1 = COS(minus_bend)

    c2, c3, c4 = SIN(bend) * u_ocross2
    c1 = COS(bend)

    # creates a unitary CA-C vector in the origin
    # the unitary vector is created to allow proper distance placement
    # in the final step
    # p1 is not need because it is 0 according to Quaternion math
    p2, p3, p4 = o1 / NORM(o1)

    # rotates the above copy for O and OXT
    # the rotation creates generates new coordinates
    m1, m2, m3, m4 = QM(
        *QM(b1, b2, b3, b4, 0, p2, p3, p4),
        b1, -b2, -b3, -b4,
        )

    n1, n2, n3, n4 = QM(
        *QM(c1, c2, c3, c4, 0, p2, p3, p4),
        c1, -c2, -c3, -c4,
        )

    # creates the final coordinate arrays
    O_coords = ARRAY((m2 * distance, m3 * distance, m4 * distance)) + C_term
    OXT_coords = ARRAY((n2 * distance, n3 * distance, n4 * distance)) + C_term

    return O_coords, OXT_coords


@njit
def calc_all_vs_all_dists_square(coords):
    """
    Calculate the upper half of all vs. all distances squared.

    Reproduces the operations of scipy.spatial.distance.pdist
    but without applying the sqrt().

    Parameters
    ----------
    coords : np.ndarray, shape (N, 3), dtype=np.float64

    Returns
    -------
    np.ndarray, shape ((N * N - N) // 2,), dytpe=np.float64
    """
    len_ = coords.shape[0]
    shape = ((len_ * len_ - len_) // 2,)
    results = np.empty(shape, dtype=np.float64)

    c = 1
    i = 0
    for a in coords:
        for b in coords[c:]:
            x = b[0] - a[0]
            y = b[1] - a[1]
            z = b[2] - a[2]
            results[i] = x * x + y * y + z * z
            i += 1
        c += 1

    return results


# njit available
def calc_all_vs_all_dists(coords):
    """
    Calculate the upper half of all vs. all distances.

    Reproduces the operations of scipy.spatial.distance.pdist.

    Parameters
    ----------
    coords : np.ndarray, shape (N, 3), dtype=np.float64

    Returns
    -------
    np.ndarray, shape ((N * N - N) // 2,), dytpe=np.float64
    """
    len_ = coords.shape[0]
    shape = ((len_ * len_ - len_) // 2,)
    results = np.empty(shape, dtype=np.float64)

    c = 1
    i = 0
    for a in coords:
        for b in coords[c:]:
            x = b[0] - a[0]
            y = b[1] - a[1]
            z = b[2] - a[2]
            results[i] = (x * x + y * y + z * z) ** 0.5
            i += 1
        c += 1

    return results


# def calc_vdW_AB(sigma_i, sigma_j, eps_i, eps_j, alpha=0.8):
#     """
#     non vectorized
#     """
#     rminA = (2 * sigma_i)**(1/6)
#     rminB = (2 * sigma_j)**(1/6)
#     rmin_pair = alpha * (rminA + rminB)
#     eps_pair = (eps_i * atomB_j)**0.5
#     A = eps_pair * rmin_pair**12
#     B = 2 * eps_pair * rmin_pair**6
#     return A, B


# def calc_Coulomb_ij(qi, qj, rij):
#     return qi * qj / rij


# def calc_all_Coulom(partial_charges, r_pairs, ep=4):
#     """
#     al implementar rever que no se computa contra si mismo
#     i < j
#
#     Where ep is the dieletric constant
#     """
#     coulombs = partial_charges[1:] * partial_charges[:-1]
#     coef = coulombs / r_pairs
#     energy_pair = coef / ep
#     return sum(energy_pair)


# def calc_FGB(
#         const=-0.5 * (1 / 4 - 1 / 80),
#         ):
#     """."""
#
#     coulombs = partial_charges[1:] * partial_charges[:-1]
#     return const * sum(qi*qj/fGB(rij) for i in range(1))  #complete
#
#
# def calc_fGB():
#     """."""
#     rij2 = rij**2
#     (rij2+Ri*Rj*math.exp(-rij2/4*Ri*Rj))**0.5
#     return


@njit
def calc_residue_num_from_index(i, step=3):
    """Calculate residue number from index.

    Parameters
    ----------
    i : int
        The index.

    step : int, optional
        How many residue numbers define a residue.
        Defaults to 3.

    Returns
    -------
    int
        The corresponding residue number for the index.
    """
    return i // step


# njit available
def round_radian_to_degree_bin_10(x0):
    """
    Round RADIAN to the nearest 10 degree bin.

    23 degrees round to 20 degrees.
    27 degrees round to 30 degrees.

    Parameters
    ----------
    x0 : float
        The angle in radians.

    Returns
    -------
    int
        The nearest bind of 10 degrees according to rounding rules.
    """
    x = int(round((x0 * 180 / 3.141592653589793), 0))
    mod_ = x % 10

    if mod_ < 5:
        return x - mod_

    elif mod_ > 5:
        return x + (10 - mod_)

    elif mod_ == 5:

        x10 = x // 10

        if x10 % 2:
            return x + 5
        else:
            return x - 5
    else:
        raise Exception('Code should not reach this point.')


def unit_vector(vector):
    """Calculate the unitary vector."""
    return vector / np.linalgn.norm(vector)


def init_lennard_jones_calculator(acoeff, bcoeff):
    """
    Calculate Lennard-Jones full pontential.

    The LJ potential is calculated fully and no approximations to
    proximity of infinite distance are considered.

    Parameters
    ----------
    acoeff, bcoeff : np.ndarray, shape (N, 3), dtype=np.float
        The LJ coefficients prepared already for the ij-pairs upon which
        the resulting function is expected to operate.
        IMPORTANT: it is up to the user to define the coefficients such
        that resulting energy is np.nan for non-relevant ij-pairs, for
        example, covalently bonded pairs, or pairs 2 bonds apart.

    Returns
    -------
    numba.njitted func
        Function closure with registered `acoeff`s and `bcoeff`s that
        expects an np.ndarray of distances with same shape as `acoeff`
        and `bcoeff`: (N,).
        `func` returns an integer.
    """
    @njit
    def calc_lennard_jones(distances_ij, NANSUM=np.nansum):
        ar = acoeff / (distances_ij ** 12)
        br = bcoeff / (distances_ij ** 6)
        energy_ij = ar - br
        return NANSUM(energy_ij)
    return calc_lennard_jones


def init_coulomb_calculator(charges_ij):
    """
    Calculate Coulomb portential.

    Parameters
    ----------
    charges_ij : np.ndarray, shape (N, 3), dtype=np.float
        The `charges_ij` prepared already for the ij-pairs upon which
        the resulting function is expected to operate.
        IMPORTANT: it is up to the user to define the charge such
        that resulting energy is np.nan for non-relevant ij-pairs, for
        example, covalently bonded pairs, or pairs 2 bonds apart.

    Returns
    -------
    numba.njitted func
        Function closure with registered `charges_ij` that expects an
        np.ndarray of distances with same shape as `acoeff` and `bcoeff`:
        (N,).
        `func` returns an integer.
    """
    @njit
    def calculate(distances_ij, NANSUM=np.nansum):
        return NANSUM(distances_ij / charges_ij)
    return calculate


def energycalculator_ij(distf, efuncs):
    """
    Calculate the sum of energy terms.

    This function works as a closure.

    Accepts only energy terms that compute for non-redundant ij-pairs.

    Energy terms must have distances ij as unique positional parameter,
    and should return an integer.

    Example
    -------
    >>> ecalc = energycalculator_ij(calc_ij_pair_distances, [...])
    >>> total_energy = ecalc(coords)

    Where `[...]` is a list containing energy term functions.

    See Also
    --------
    init_lennard_jones_calculator
    init_coulomb_calculator

    Parameters
    ----------
    distf : func
        The function that will be used to calculate ij-pair distances
        on each call. If performance is a must, this function should be
        fast. `distf` function should receive `coords` as unique
        argument where `coords` is a np.ndarray of shape (N, 3), where N
        is the number of atoms, and 3 represents the XYZ coordinates.
        This function should return a np.ndarray of shape
        (N * (N - 1)) / 2,), dtype=np.float.

    efuncs : list
        A list containing the energy terms functions. Energy term
        functions are prepared closures that accept the output of
        `distf` function.

    Returns
    -------
    func
        A function that accepts coords in the form of (N, 3). The
        coordinates sent to the resulting function MUST be aligned with
        the labels used to prepare the `efuncs` closures.
    """
    def calculate(coords):
        dist_ij = distf(coords)
        energy = 0
        for func in efuncs:
            energy += func(dist_ij)
        return energy
    return calculate


# njit
def sum_upper_diagonal_raw(data, result):
    """
    Calculate outer sum for upper diagonal with for loops.

    The use of for-loop based calculation avoids the creation of very
    large arrays using numpy outer derivates. This function is thought
    to be jut compiled.

    Does not create new data structure. It requires the output structure
    to be provided. Hence, modifies in place. This was decided so
    because this function is thought to be jit compiled and errors with
    the creation of very large arrays were rising. By passing the output
    array as a function argument, errors related to memory freeing are
    avoided.

    Parameters
    ----------
    data : an interable of Numbers, of length N

    result : a mutable sequence, either list of np.ndarray,
             of length N*(N-1)//2
    """
    c = 0
    len_ = len(data)
    for i in range(len_ - 1):
        for j in range(i + 1, len_):
            result[c] = data[i] + data[j]
            c += 1

    # assert result.size == (data.size * data.size - data.size) // 2
    # assert abs(result[0] - (data[0] + data[1])) < 0.0000001
    # assert abs(result[-1] - (data[-2] + data[-1])) < 0.0000001
    return


# njit available
def multiply_upper_diagonal_raw(data, result):
    """
    Calculate the upper diagonal multiplication with for loops.

    The use of for-loop based calculation avoids the creation of very
    large arrays using numpy outer derivatives. This function is thought
    to be njit compiled.

    Does not create new data structure. It requires the output structure
    to be provided. Hence, modifies in place. This was decided so
    because this function is thought to be jit compiled and errors with
    the creation of very large arrays were rising. By passing the output
    array as a function argument, errors related to memory freeing are
    avoided.

    Parameters
    ----------
    data : an interable of Numbers, of length N

    result : a mutable sequence, either list of np.ndarray,
             of length N*(N-1)//2
    """
    c = 0
    len_ = len(data)
    for i in range(len_ - 1):
        for j in range(i + 1, len_):
            result[c] = data[i] * data[j]
            c += 1

    # assert result.size == (data.size * data.size - data.size) // 2
    # assert abs(result[0] - data[0] * data[1]) < 0.0000001
    # assert abs(result[-1] - data[-2] * data[-1]) < 0.0000001
    return


def make_seq_probabilities(seq, reverse=True):
    """Make probabilites from a sequence of numbers."""
    sum_ = sum(seq)
    probs = np.array(seq) / sum_
    if reverse:
        return probs
    else:
        return probs[::-1]


calc_all_vs_all_dists_njit = njit(calc_all_vs_all_dists)
multiply_upper_diagonal_raw_njit = njit(multiply_upper_diagonal_raw)
rotate_coordinates_Q_njit = njit(rotate_coordinates_Q)
rrd10_njit = njit(round_radian_to_degree_bin_10)
sum_upper_diagonal_raw_njit = njit(sum_upper_diagonal_raw)
unit_vec_njit = njit(unit_vector)
