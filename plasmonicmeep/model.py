"""

Create the vertices for MEEP to use as a material.

S.D.G.

"""

from functools import partial
from typing import List, Optional

import meep as mp
import numpy as np


def bow_y(
    height: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Prism]:
    """Create bowtie-shaped dimer oriented along y axis.

    Args:
        height (float): The height of the individual prisms making up the dimer.
        separation (float): The separation distance between the prisms.
        center (mp.Vector3): The coordinates of the center of the dimer.
        material (mp.Medium): The material of the dimer.

    Returns:
        List[mp.Prism]: List of geometric structures making up the dimer being simulated.
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
) -> List[mp.Prism]:
    """Create bowtie-shaped dimer oriented along x axis.

    Args:
        height (float): The height of the individual prisms making up the dimer.
        separation (float): The separation distance between the prisms.
        center (mp.Vector3): The coordinates of the center of the dimer.
        material (mp.Medium): The material of the dimer.

    Returns:
        List[mp.Prism]: List of geometric structures making up the dimer being simulated.
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


def single_facing_triangle(
        height: float, separation: float, center: mp.Vector3, material: mp.Medium, direction: str = "y"
) -> List[mp.Sphere]:
    """Create dimer of spherical nanoparticles oriented along y axis, but skipping the second nanoparticle.

    Args:
        radius (float): The height of the triangle.
        separation (float): Twice the distance between the triangle and the center.
        center (mp.Vector3): The coordinates of the center of the structure.
        material (mp.Medium): The material of the triangle.


    Returns:
        List[mp.Sphere]: Single-item list of the nanotriangle.
    """

    # TODO: Fix this!

    edge_length = np.sqrt(height ** 2 / 3)

    if direction == "y":
        vertices = [
            mp.Vector3(center.x, center.y + separation / 2),
            mp.Vector3(center.x + edge_length / 2, center.y + height + separation / 2),
            mp.Vector3(center.x - edge_length / 2, center.y + height + separation / 2),
        ]
    else:
        vertices = [
            mp.Vector3(center.x + separation / 2, center.y),
            mp.Vector3(center.x + height + separation / 2, center.y + edge_length / 2),
            mp.Vector3(center.x + height + separation / 2, center.y - edge_length / 2),
        ]

    return [mp.Prism(vertices=vertices, height=mp.inf, material=material)]



def single_sphere(
        radius: float, separation: float, center: mp.Vector3, material: mp.Medium, direction: str = "y",
) -> List[mp.Sphere]:
    """Create dimer of spherical nanoparticles oriented along y axis, but skipping the second nanoparticle.

    Args:
        radius (float): The radius of the sphere.
        separation (float): Twice the distance between the sphere and the center.
        center (mp.Vector3): The coordinates of the center of the structure.
        material (mp.Medium): The material of the sphere.

    Returns:
        List[mp.Sphere]: Single-item list of the nanosphere.
    """

    if direction == "y":
        xy = center.x, center.y + radius + separation / 2
    else:
        xy = center.x + radius + separation / 2, center.y

    return [mp.Sphere(center=mp.Vector3(*xy), radius=radius, material=material)]


def spheres_y(
    radius: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Sphere]:
    """Create dimer of spherical nanoparticles oriented along y axis.

    Args:
        radius (float): The radius of the individual spheres making up the dimer.
        separation (float): The separation distance between the spheres.
        center (mp.Vector3): The coordinates of the center of the dimer.
        material (mp.Medium): The material of the dimer.

    Returns:
        List[mp.Sphere]: List of geometric structures making up the dimer being simulated.
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
        separation (float): The separation distance between the spheres.
        center (mp.Vector3): The coordinates of the center of the dimer.
        material (mp.Medium): The material of the dimer.

    Returns:
        List[mp.Sphere]: List of geometric structures making up the dimer being simulated.
    """

    # Calculate spherical coordinates
    sphere1_pos = mp.Vector3(center.x - radius - separation / 2, center.y)
    sphere2_pos = mp.Vector3(center.x + radius + separation / 2, center.y)

    sphere1 = mp.Sphere(center=sphere1_pos, radius=radius, material=material)
    sphere2 = mp.Sphere(center=sphere2_pos, radius=radius, material=material)

    return [sphere1, sphere2]


def invertedtr_y(
    height: float, separation: float, center: mp.Vector3, material: mp.Medium,
) -> List[mp.Prism]:
    """Create dimer of inverted triangular nanoparticles oriented along y axis.

    Args:
        height (float): The height of the individual prisms making up the dimer.
        separation (float): The separation distance between the prisms.
        center (mp.Vector3): The coordinates of the center of the dimer.
        material (mp.Medium): The material of the dimer.

    Returns:
        List[mp.Prism]: List of geometric structures making up the dimer being simulated.
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
) -> List[mp.Prism]:
    """Create dimer of inverted triangular nanoparticles oriented along x axis.

    Args:
        height (float): The height of the individual prisms making up the dimer.
        separation (float): The separation distance between the prisms.
        center (mp.Vector3): The coordinates of the center of the dimer.
        material (mp.Medium): The material of the dimer.

    Returns:
        List[mp.Prism]: List of geometric structures making up the dimer being simulated.
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
    "triangle-tip-y": partial(single_facing_triangle, direction="y"),
    "triangle-tip-x": partial(single_facing_triangle, direction="x"),
    "spheres-y": spheres_y,
    "spheres-x": spheres_x,
    "inv-bow-y": invertedtr_y,
    "inv-bow-x": invertedtr_x,
    "sphere-y": partial(single_sphere, direction="y"),
    "sphere-x": partial(single_sphere, direction="x"),
    }
