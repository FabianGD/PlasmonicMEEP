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


def bow_y(
    height: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Vector3]:
    """Create bowtie-shaped dimer oriented along y axis.

    Args:
        height (float): The height of the individual prisms making up the dimer. 
            Required argument (no default).
        separation (float): The separation distance between the prisms. 
            Required argument (no default).
        center (mp.Vector3): The coordinates of the center of the dimer. 
            Required argument (no default).
        material (mp.Medium): The material of the dimer. Required argument (no default).

    Returns:
        List[mp.Vector3]: List of geometric structures making up the dimer being simulated.
    """

    # Calculate vertices
    side_length = np.sqrt(height ** 2 / 3)

    vertices1 = [
        mp.Vector3(center.x, center.y - separation / 2),
        mp.Vector3(center.x - side_length, center.y - separation / 2 - height),
        mp.Vector3(center.x + side_length, center.y - separation / 2 - height),
    ]
    vertices2 = [
        mp.Vector3(center.x, center.y + separation / 2),
        mp.Vector3(center.x - side_length, center.y + separation / 2 + height),
        mp.Vector3(center.x + side_length, center.y + separation / 2 + height),
    ]

    triangle1 = mp.Prism(vertices1, height=mp.inf, material=material)
    triangle2 = mp.Prism(vertices2, height=mp.inf, material=material)

    return [triangle1, triangle2]


def bow_x(
    height: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Vector3]:
    """Create bowtie-shaped dimer oriented along x axis.

    Args:
        height (float): The height of the individual prisms making up the dimer. 
            Required argument (no default).
        separation (float): The separation distance between the prisms. 
            Required argument (no default).
        center (mp.Vector3): The coordinates of the center of the dimer. 
            Required argument (no default).
        material (mp.Medium): The material of the dimer. Required argument (no default).

    Returns:
        List[mp.Vector3]: List of geometric structures making up the dimer being simulated.
    """

    # Calculate vertices
    side_length = np.sqrt(height ** 2 / 3)

    vertices1 = [
        mp.Vector3(center.x - separation / 2, center.y),
        mp.Vector3(center.x - separation / 2 - height, center.y - side_length),
        mp.Vector3(center.x - separation / 2 - height, center.y + side_length),
    ]
    vertices2 = [
        mp.Vector3(center.x + separation / 2, center.y),
        mp.Vector3(center.x + separation / 2 + height, center.y - side_length),
        mp.Vector3(center.x + separation / 2 + height, center.y + side_length),
    ]

    triangle1 = mp.Prism(vertices1, height=mp.inf, material=material)
    triangle2 = mp.Prism(vertices2, height=mp.inf, material=material)

    return [triangle1, triangle2]


def spheres_y(
    radius: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Sphere]:
    """Create dimer of spherical nanoparticles oriented along y axis.

    Args:
        radius (float): The radius of the individual spheres making up the dimer. 
            Required argument (no default).
        separation (float): The separation distance between the spheres. 
            Required argument (no default).
        center (mp.Vector3): The coordinates of the center of the dimer. 
            Required argument (no default).
        material (mp.Medium): The material of the dimer. Required argument (no default).

    Returns:
        List[mp.Vector3]: List of geometric structures making up the dimer being simulated.
    """

    # Calculate spherical coordinates
    sphere1_pos = mp.Vector3(center.x, center.y - radius - separation / 2)
    sphere2_pos = mp.Vector3(center.x, center.y + radius + separation / 2)

    sphere1 = mp.Sphere(center=sphere1_pos, radius=radius, material=material)
    sphere2 = mp.Sphere(center=sphere2_pos, radius=radius, material=material)

    return [sphere1, sphere2]


def spheres_x(
    radius: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Sphere]:
    """Create dimer of spherical nanoparticles oriented along x axis.

    Args:
        radius (float): The radius of the individual spheres making up the dimer.  
            Required argument (no default).
        separation (float): The separation distance between the spheres. 
            Required argument (no default).
        center (mp.Vector3): The coordinates of the center of the dimer. 
            Required argument (no default).
        material (mp.Medium): The material of the dimer. Required argument (no default).

    Returns:
        List[mp.Vector3]: List of geometric structures making up the dimer being simulated.
    """

    # Calculate spherical coordinates
    sphere1_pos = mp.Vector3(center.x - radius - separation / 2, center.y)
    sphere2_pos = mp.Vector3(center.x + radius + separation / 2, center.y)

    sphere1 = mp.Sphere(center=sphere1_pos, radius=radius, material=material)
    sphere2 = mp.Sphere(center=sphere2_pos, radius=radius, material=material)

    return [sphere1, sphere2]


def invertedtr_y(
    height: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Vector3]:
    """Create dimer of inverted triangular nanoparticles oriented along y axis.

    Args:
        height (float): The height of the individual prisms making up the dimer. 
            Required argument (no default).
        separation (float): The separation distance between the prisms. 
            Required argument (no default).
        center (mp.Vector3): The coordinates of the center of the dimer. 
            Required argument (no default).
        material (mp.Medium): The material of the dimer. Required argument (no default).

    Returns:
        List[mp.Vector3]: List of geometric structures making up the dimer being simulated.
    """

    # Calculate vertices
    side_length = np.sqrt(height ** 2 / 3)

    vertices1 = [
        mp.Vector3(center.x, center.y - separation / 2 - height),
        mp.Vector3(center.x - side_length, center.y - separation / 2),
        mp.Vector3(center.x + side_length, center.y - separation / 2),
    ]
    vertices2 = [
        mp.Vector3(center.x, center.y + separation / 2 + height),
        mp.Vector3(center.x - side_length, center.y + separation / 2),
        mp.Vector3(center.x + side_length, center.y + separation / 2),
    ]

    triangle1 = mp.Prism(vertices1, height=mp.inf, material=material)
    triangle2 = mp.Prism(vertices2, height=mp.inf, material=material)

    return [triangle1, triangle2]


def invertedtr_x(
    height: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Vector3]:
    """Create dimer of inverted triangular nanoparticles oriented along x axis.

    Args:
        height (float): The height of the individual prisms making up the dimer. 
            Required argument (no default).
        separation (float): The separation distance between the prisms. 
            Required argument (no default).
        center (mp.Vector3): The coordinates of the center of the dimer. 
            Required argument (no default).
        material (mp.Medium): The material of the dimer. Required argument (no default).

    Returns:
        List[mp.Vector3]: List of geometric structures making up the dimer being simulated.
    """

    # Create vertices
    side_length = np.sqrt(height ** 2 / 3)

    vertices1 = [
        mp.Vector3(center.x - separation / 2 - height, center.y),
        mp.Vector3(center.x - separation / 2, center.y - side_length),
        mp.Vector3(center.x - separation / 2, center.y + side_length),
    ]
    vertices2 = [
        mp.Vector3(center.x + separation / 2 + height, center.y),
        mp.Vector3(center.x + separation / 2, center.y - side_length),
        mp.Vector3(center.x + separation / 2, center.y + side_length),
    ]

    triangle1 = mp.Prism(vertices1, height=mp.inf, material=material)
    triangle2 = mp.Prism(vertices2, height=mp.inf, material=material)

    return [triangle1, triangle2]


MODEL_MAPPING = {
    "bowtie-y": bow_y,
    "bowtie-x": bow_x,
    "spheres-y": spheres_y,
    "spheres-x": spheres_x,
    "inv-bow-y": invertedtr_y,
    "inv-bow-x": invertedtr_x,
}

if __name__ == "__main__":
    bow_y
