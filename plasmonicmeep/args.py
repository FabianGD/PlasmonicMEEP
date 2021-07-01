"""
S.D.G

Argument parsers for the plasmonicmeep scripts

"""

import argparse


def fdtd_argparsing():
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(
        description=(
            "Main FDTD module. Computes spectra and densities of a given model. "
            "May be run in parallel using MPI."
        )
    )

    parser.add_argument(
        "-n",
        "--calc-name",
        default="plasmonic",
        type=str,
        help=(
            "Name or identifier of the calculation. Used as filename_prefix in the MEEP simulation."
        ),
    )

    parser.add_argument(
        "-x",
        "--sizex",
        type=float,
        default=0.5,
        help="The size of the box in µm in the x direction.",
    )

    parser.add_argument(
        "-y",
        "--sizey",
        type=float,
        default=0.5,
        help="The size of the box in µm in the y direction.",
    )

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=200,
        help="The resolution of the box (nr. of pixels per µm?)",
    )

    parser.add_argument(
        "-f",
        "--frequency",
        type=float,
        default=1.5,
        help=(
            "Change the central frequency of the incident laser field. "
            "Frequency is given in units of 1/µm."
        ),
    )

    parser.add_argument(
        "-w",
        "--freq-width",
        type=float,
        default=1.5,
        help=(
            "Change the frequency width of the incident laser field. "
            "Pulse with is given in units of µm."
        ),
    )

    parser.add_argument(
        "-u",
        "--run-until",
        default=20.0,
        type=float,
        help=(
            "Set the finish time of the simulations in MEEP time units. Defaults to 20.0"
        ),
    )

    parser.add_argument(
        "-t",
        "--cfield-time-step",
        default="auto",
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

    return parser.parse_args()
