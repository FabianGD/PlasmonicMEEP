"""

S.D.G.

"""

from pathlib import Path
from typing import Union

import yaml
import attr
import meep as mp


def _unit_converter(value):
  ...


@attr.s()
class Point2D:
    """
    Point implementation as attrs class.
    """
    x: float = attr.ib(converter=_unit_converter, default=0.0)
    y: float = attr.ib(converter=_unit_converter, default=0.0)

    def as_vector3(self):
        return mp.Vector3(self.x, self.y, 0.0)


@attr.s()
class Geometry:
    """

    """
    anchor: str = attr.ib()
    material: Union[str, mp.Medium] = attr.ib()


@attr.s()
class Sphere(Geometry):
    """

    """
    size: Point2D = attr.ib(default=(Point2D()))
    center: str = attr.ib()


@attr.s()
class Cell:
    """

    """
    size: Point2D = attr.ib()
    origin: str = attr.ib(validator=attr.validators.in_(["center", "lower left", "upper left"]))


def read_geom(file: Path):
    """

    """

    with open(file, "r") as f:
        data = yaml.load(f.read(), loader=yaml.FullLoader)


if __name__ == "__main__":

    data = """
    resolution: 1000
    cell:
      size:
        x: 1.0
        y: 1.0
      origin: center
    geometry:
      - sphere:
        position:
          x: 0.0
          y: 0.0
        anchor:
        material:
    """