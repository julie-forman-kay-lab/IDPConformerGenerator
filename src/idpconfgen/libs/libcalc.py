"""Mathemical calculations."""
import math
from math import sin, cos, pi
from numba import njit

import numpy as np

from idpconfgen.core.definitions import (
    build_bend_CA_C_O,
    build_bend_CA_C_OXT,
    distance_C_O,
    distance_C_OXT,
    )


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


def get_separate_torsions(torsions_array):
    """
    Separate torsion angles accorindg to the protein backbone concept.

    Considers torsions angles for bonds in between atom pairs:
        - CA - C
        - C - N
        - N - CA

    Backbone obeys the order: N-CA-C-N-CA-C(...)

    And the firt value corresponds to a CA-C pair, because the
    first N-CA pair of the protein backbone has no torsion angle.
    """
    assert torsions_array.ndim == 1
    assert torsions_array.size % 3 == 0

    CA_C = torsions_array[::3].tolist()
    C_N = torsions_array[1::3].tolist()
    N_CA = torsions_array[2::3].tolist()

    assert len(CA_C) == len(C_N) == len(N_CA)
    return CA_C, C_N, N_CA


def validate_backbone_labels_for_torsion(labels, minimum=2):
    """
    Validate labels for torsion angle calculation.

    Assumes labels are aligned with their corresponding coordinates.
    Yet, coordinates have no scope in this function.

    Excepts only the mininal backbone labels, these are: N, CA, and C.

    Parameters
    ----------
    labels : np.array of shape (N,) or alike
        Where N % 3 equals 0.

    minimum : int
        The minimum number of residues to consider valid.
    """
    if len(labels) / 3 < minimum:
        return 'Too small segment'

    if labels[0] != 'N':
        return 'The first atom is not N, it should be!'

    if len(labels) % 3:
        return 'Number of backbone atoms is not module of 3.'

    if set(labels) != {'N', 'C', 'CA'}:
        return 'There are atoms other than N, C and CA.'

    return ''


def calc_MSMV(data):
    """Calculate Mean, STD, Median, and Variance."""
    return np.mean(data), np.std(data), np.median(data), np.var(data)


@njit
def hamiltonian_multiplication_Q(a1, b1, c1, d1, a2, b2, c2, d2):
    """Hamiltonian Multiplication."""
    return (
            a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    )


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
        calcualtions and just provide the torsion angle value as is.

    Returns
    -------
    np.ndarray of shape (3,), dtype=np.float32
    """
    # transfer vectors to origin
    o1 = v1 - v2
    o2 = v3 - v2

    ocross = CROSS(o2, o1)  #changed
    u_ocross = ocross / NORM(ocross)

    # creates quaterion to rotate on the bend angle
    b2, b3, b4 = SIN(bend) * u_ocross
    b1 = COS(bend)

    # rotates the unitary of o2 according to bend angle
    uo2 = o2  / NORM(o2)
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
def make_coord_Q_CO(
        CA_coords,
        C_coords,
        N_coords,
        distance=distance_C_O,
        bend=build_bend_CA_C_O,
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
    CA_coords, C_coods, N_coords : np.ndarray, shape (3,), dtype=np.float
        The XYZ coordinates for the CA, C and N atoms surrounding the
        C=O bond, respectively.

    distance : float, optional
        The distance of the C-O bond pair.
        Defaults to ~1.234.

    bend : float, optional
        The angle in radians for the CA-C-O.
        Defaults to ~2.1428 radians (~122.777 degress).
        If `bend` is given, consider the actual `bend` value must be
        `(pi - bend) / 2`, this calculation must be computed outside
        this function for perfomance reasons.

    Returns
    -------
    np.ndarray of shape (3,), dtype=np.float32
    """
    o1 = CA_coords - C_coords
    o2 = N_coords - C_coords

    ocross = CROSS(o1, o2)
    u_ocross = ocross / NORM(ocross)

    # creates quaterion to rotate on the bend angle
    b2, b3, b4 = SIN(bend) * u_ocross
    b1 = COS(bend)

    # rotates a the unitary of o2 according to bend angle
    uo2 = o2 / NORM(o2)
    p2, p3, p4 = uo2  # p1 is zero according to Quaternion theory
    n1, n2, n3, n4 = QM(
        *QM(b1, b2, b3, b4, 0, p2, p3, p4),
        b1, -b2, -b3, -b4,
        )

    return ARRAY((n2 * distance, n3 * distance, n4 * distance)) + C_coords


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
