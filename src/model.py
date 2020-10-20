"""
Create the vertices for MEEP to use as a material
"""

import meep as mp
import numpy as np


def create_sinus_grating(
    ampl, periodicity, thickness, resolution, sizex, y=False, matrix=None
):
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

        # TODO Understand whatever that does.
        if matrix is not None:
            trfrmd = np.dot(matrix, [x_coord, y_coord])
            x_coord = trfrmd[0]
            y_coord = trfrmd[1]

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
            trfrmd = np.dot(matrix, [x_coord, y_coord])
            x_coord = trfrmd[0]
            y_coord = trfrmd[1]

        if y:
            aux_vert.append(mp.Vector3(y_coord, x_coord))
        else:
            aux_vert.append(mp.Vector3(x_coord, y_coord))

    # print(list(zip(vertices, aux_vert[::-1])))

    return vertices + aux_vert
