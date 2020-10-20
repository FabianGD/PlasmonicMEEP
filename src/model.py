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

    # print(list(zip(vertices, aux_vert[::-1])))

    return vertices + aux_vert


def create_np_on_mirror(radius, separation, material, size_x, size_y, y=False):

    # TODO

    # Generate the slab coordinates
    if y:
        block_width = size_x / 2
        block_size = mp.Vector3(block_width, mp.inf, mp.inf)
        block_center = mp.Vector3(size_x / 2, 0, 0)
        sphere1_pos = mp.Vector3(
            block_center.x - block_width / 2 - radius - separation, 0, 0
        )
        sphere2_pos = mp.Vector3(
            block_center.x - block_width / 2 + radius, 0, 0
        )
    else:
        block_width = size_y / 2
        block_size = mp.Vector3(mp.inf, block_width, mp.inf)
        block_center = mp.Vector3(0, 1 / block_width, 0)
        sphere_pos = mp.Vector3(
            0, block_center.x - block_width / 2 - radius - separation, 0
        )

    # block = mp.Block(block_size, center=block_center, material=material)
    sphere1 = mp.Sphere(center=sphere1_pos, radius=radius, material=material)
    sphere2 = mp.Sphere(center=sphere2_pos, radius=radius, material=material)

    return [sphere1, sphere2]
