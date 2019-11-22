"""Mathemical calculations."""
import numpy as np


AXIS_111 = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        ])


def make_axis_vectors(A, B, C):
    """
    Make axis vectors.
    
    For np.array([x,y,z]) of connected atoms A->B->C.
    
    Credits to @RobertVernon.

    Example:
        >>> av = np.array([0.000, 0.000, 0.000])
        >>> bv = np.array([1.458, 0.000, 0.000])
        >>> cv = np.array([2.009, 1.420, 0.000])
        >>> make_axis_vectors(av, bv, cv)
        >>> (array([ 0.,  0., -1.]), array([-1.,  0.,  0.]), array([-0., -1., -0.]))

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



def RT_to_plane(A, B, C):# PARENT, XAXIS, PLANE
    """
    Define rotations and translation matrices for planar orientation.

    Where, A->B determine the X axis, Y is the normal vector to the ABC
        plane, and A is defined at 0,0,0.

    Uses make_axis_vectors()

    Parameters
    ----------
    A, B, C: np.array of shape (3,)
        The 3D coordinates.

    Returns
    -------
    tuple
        Rotational matrix and A
    """
    parallel_ABC, AB_vect, perpendicular = make_axis_vectors(A, B, C)
    b = np.array([AB_vect, -perpendicular, parallel_ABC])
    
    v = AXIS_111
    
    return np.dot(b, v), A


def make_coord(theta, phi, rad, parent, xaxis, yaxis):
    """
    Makes a new coordinate in space.

    Parameters
    ----------
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
    new_coord = np.array([
        rad * math.cos(phi),
        rad * math.sin(phi) * math.cos(theta),
        rad * math.sin(phi) * math.sin(theta),
        ])

    rotation, translation = RT_to_plane(parent, xaxis, yaxis)
    new_coord = np.dot(new_coord, rotation.T)
    return new_coord + translation
