"""
Argument parsers for the plasmonicmeep scripts

S.D.G
"""

import argparse
from functools import partial
from typing import List, Any, Optional, Sequence, Type
from meep import Vector3

from .model import MODEL_MAPPING


class VectorAction(argparse.Action):
    """
    Custom argparse.Action subclass to generate a meep.Vector3 directly
    from input.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):

        vec = Vector3(*values)
        setattr(namespace, self.dest, vec)


def positive_type(value: str, rtype: Type, other_allowed: List[Any] = None) -> Any:
    """Type conversion for positive numericals (floats/ints). Allows specified other values.

    Args:
        value (str): The value to convert and check
        rtype (Type): The type conversion function to use.
        other_allowed (List[Any], optional): Other allowed values for
            special cases, such as "auto" or "-1" etc. Defaults to [].

    Raises:
        argparse.ArgumentTypeError: If the value provided is not in `other_allowed` and
            negative.
        ValueError: If the value provided is not in `other_allowed` and could not be
            converted to the requested value.

    Returns:
        Any: Either the value as is (if in `other_allowed`) or the value after
            conversion. It will also return the converted value if the converted
            value is in `other_allowed`.
    """

    if other_allowed is None:
        other_allowed = []

    if value in other_allowed:
        return value

    try:
        conv_val = rtype(value)
    except ValueError as e:
        raise ValueError(
            "Could not convert value '{}' to type '{}'".format(value, rtype.__name__)
        ) from e

    if conv_val in other_allowed:
        return conv_val
    elif conv_val <= 0:
        raise argparse.ArgumentTypeError(
            "{} has to be a positive {} value.".format(value, rtype.__name__)
        )

    return conv_val


def fdtd_argparsing(args: Optional[Sequence[str]] = None):
    """
    Argument parser for the FDTD module.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Main FDTD module. Computes spectra and densities of a given model. "
            "May be run in parallel using MPI."
        )
    )

    parser.add_argument(
        "structure",
        choices=tuple(MODEL_MAPPING.keys()),
        type=str,
        help=(
            "Structure to use for the FDTD simulation. Can be easily extended by editing plasmonicmeep/model.py."
        )
        # default="spheres_y"
    )

    parser.add_argument(
        "-R",
        "--radius",
        default=0.05,
        type=float,
        help=(
            "Radius of the structure in µm. For the triangular structures, the radius "
            "of the circumscribing circle. Defaults to 0.05."
        )
    )
    parser.add_argument(
        "-S",
        "--separation",
        default=0.005,
        type=float,
        help=(
            "Separation of the (two) structure in µm. Defaults to 0.005."
        )
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        help=(
            "Get more verbose output. Maximum two times."
        ),
        default=0
    )

    parser.add_argument(
        "-n",
        "--calc-name",
        type=str,
        default="PlasmonicMEEP",
        help=(
            "Name or identifier of the calculation. Used as filename_prefix in the MEEP simulation."
        ),
    )

    parser.add_argument(
        "-x",
        "--sizex",
        type=partial(positive_type, rtype=float),
        default=0.5,
        help="The size of the box in µm in the x direction. Defaults to 0.5",
    )

    parser.add_argument(
        "-y",
        "--sizey",
        type=partial(positive_type, rtype=float),
        default=0.5,
        help="The size of the box in µm in the y direction. Defaults to 0.5",
    )

    parser.add_argument(
        "-r",
        "--resolution",
        type=partial(positive_type, rtype=int),
        default=200,
        help="The resolution of the box (nr. of pixels per µm). Defaults to 200.",
    )

    parser.add_argument(
        "-f",
        "--frequency",
        type=partial(positive_type, rtype=float),
        default=1.5,
        help=(
            "Change the central frequency of the incident laser field. "
            "Frequency is given in units of 1/µm. Defaults to 1.5 (667 nm)."
        ),
    )

    parser.add_argument(
        "-w",
        "--freq-width",
        type=partial(positive_type, rtype=float),
        default=1.5,
        help=(
            "Change the frequency width of the incident laser field. "
            "Pulse with is given in units of µm. Defaults to 1.5."
        ),
    )

    parser.add_argument(
        "-u",
        "--run-until",
        default=20.0,
        type=partial(positive_type, rtype=float),
        help=(
            "Set the finish time of the simulations in MEEP time units. Defaults to 20.0"
        ),
    )

    parser.add_argument(
        "-p",
        "--point",
        type=float,
        default=[0, 0],
        nargs=2,
        metavar=("X", "Y"),
        help=(
            "Specify the location of the single point on which the "
            "field is output from. Defaults to the center of the 2D box, "
            "which is at [0, 0]."
        ),
    )

    parser.add_argument(
        "-t",
        "--cfield-time-step",
        default="auto",
        type=partial(positive_type, rtype=float, other_allowed=["auto"]),
        help=(
            "Specify a custom time step to write the field in the entire cell to file. "
            "By default or by specifying 'auto', the time step is calculated "
            "by $dt = 0.5 / (cfreq + fwidth * 0.5)$. With the defaults, which is around "
            "0.22 for the default values."
        ),
    )

    parser.add_argument("-o", "--output", default="./data/", help="Output folder.")

    # Boolean args
    parser.add_argument(
        "-g",
        "--show-geometry",
        action="store_true",
        help="Plot only geometry and exit.",
    )

    parser.add_argument(
        "-s",
        "--spectrum",
        action="store_true",
        help=(
            "At the end of the simulation, plot and save transmission/reflectance/loss spectra. "
            "This might fail for MPI runs or on clusters. "
            "Adds complexity and runtime to the calculation."
        ),
    )

    parser.add_argument(
        "-c",
        "--enable-complex",
        action="store_true",
        help=(
            "Enable calculation of complex fields, useful when trying "
            "to retrieve the phase around the nanostructure"
        ),
        default=False,
    )

    parser.add_argument(
        "--disable-cell-field",
        action="store_true",
        help=(
            "Disable writing the entire field to file. Saves some computational effort."
        ),
    )
    parser.add_argument(
        "--disable-single-point",
        action="store_true",
        help=(
            "Disable writing the single point field to file. Saves some computational effort."
        ),
    )

    return parser.parse_args(args=args)
