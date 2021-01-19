"""
Create the vertices for MEEP to use as a material
"""

from typing import List, Optional

import meep as mp
import numpy as np


def create_sinus_grating(
    ampl: float,
    periodicity: float,
    thickness: float,
    resolution: float,
    sizex: float,
    y: bool = False,
    matrix: Optional[np.ndarray] = None,
) -> List[mp.Vector3]:
    """
    Create vertices for sinusoidal grating with given parameters
    If y is given, rotate it 90 degrees
    If matrix is not None, additionally transform geometry
    """

    freq = 1 / periodicity
    vertices = []
    aux_vert = []
    dx = 2 / resolution

    def f(x):
        return ampl / 2 * np.cos(2 * np.pi * x * freq)

    # Calculate vertices
    for i in np.arange(-sizex / 2 - dx, sizex / 2 + dx + dx, dx):
        x_coord = i
        y_coord = f(i)

        # Optionally transform the coordinates.
        if matrix is not None:
            transformed = np.dot(matrix, [x_coord, y_coord])
            x_coord = transformed[0]
            y_coord = transformed[1]

        # Append everything to the list of vertices
        if y:
            vertices.append(mp.Vector3(y_coord, x_coord))
        else:
            vertices.append(mp.Vector3(x_coord, y_coord))

    # Calculate vertices with thickness added
    for i in np.arange(sizex / 2 + dx, -sizex / 2 - dx - dx, -dx):
        x_coord = i
        y_coord = f(i) + thickness

        if matrix is not None:
            transformed = np.dot(matrix, [x_coord, y_coord])
            x_coord = transformed[0]
            y_coord = transformed[1]

        if y:
            aux_vert.append(mp.Vector3(y_coord, x_coord))
        else:
            aux_vert.append(mp.Vector3(x_coord, y_coord))

    return vertices + aux_vert


def two_nps(
    radius: float,
    separation: float,
    center: mp.Vector3,
    material: mp.Medium,
    y: bool = False,
) -> List[mp.Sphere]:
    """
    Create two spherical nanoparticles used in a simulation.
    """

    # Generate the slab coordinates
    if y:
        sphere1_pos = mp.Vector3(center.x - radius - separation / 2, center.y, 0)
        sphere2_pos = mp.Vector3(center.x + radius + separation / 2, center.y, 0)

    else:
        sphere1_pos = mp.Vector3(center.x, center.y - radius - separation / 2, 0)
        sphere2_pos = mp.Vector3(center.x, center.y + radius + separation / 2, 0)

    # block = mp.Block(block_size, center=block_center, material=material)
    sphere1 = mp.Sphere(center=sphere1_pos, radius=radius, material=material)
    sphere2 = mp.Sphere(center=sphere2_pos, radius=radius, material=material)

    return [sphere1, sphere2]
