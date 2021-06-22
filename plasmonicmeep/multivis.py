"""
Visualise the FDTD results.
"""

import argparse
import sys
from pathlib import Path

import yaml
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pint

try:
    from .multiviewer import multi_slice_viewer
except ImportError:
    from multiviewer import multi_slice_viewer


def argparsing():
    """
    Module-level argparsing
    """
    parser = argparse.ArgumentParser(
        description="Visualize calculated field enhancement"
    )
    parser.add_argument("file", type=str, help="YAML input file to use")

    # parser.add_argument("-k", "--xyskip", type=int, help="Number of pixels to skip")

    # parser.add_argument(
    #     "-p", "--mplstyle", type=str, help="Give a matplotlib style to use. Optional."
    # )

    # parser.add_argument(
    #     "-d",
    #     "--direction",
    #     help="Whether to take a cut along x or y",
    #     type=str,
    #     default="x",
    #     choices=["x", "y"],
    # )

    return parser.parse_args()

def open_yaml(file: Path) -> dict:
    """

    """

    with open(file, "r") as f:
        data = yaml.load(f.read())

    return data


def main():
    """
    Main function for the visualiser
    """

    args = argparsing()
    inputpath = Path(args.file)

    inputs = open_yaml(inputpath)
    direction = str(inputs.get("direction"))

    if direction not in ["x", "y"]:
        raise ValueError("Direction wrong.")

    mplstyle = inputs.get("mplstyle", None)
    if mplstyle is not None:
        plt.style.use(mplstyle)

    # Set up a plot
    fig, ax = plt.subplots(nrows=1, ncols=1)

    for input_dict in inputs["files"]:

        file = Path(input_dict["file"])
        label = str(input_dict["label"])
        xyskip = int(input_dict.get("skip", 30))
        approx_freq = pint.Quantity(input_dict.get("freq", None)).to("nm", "sp").m
        print(approx_freq)
        # Prepare frequency h5 file
        h5f = h5py.File(file, "r")

        keys = list(h5f.keys())
        dset = [h5f[key] for key in keys]
        # There should not be any more datasets
        dset = dset[0]


        # Get the frequency info
        freqs = dset.attrs["freqs"]

        # Skip N pixels from the sides to eliminate enhancement on boundaries
        skip_x = xyskip
        skip_y = xyskip

        # Load data into memory
        dset = np.asarray(dset).transpose(1, 0, -1)

        # Get the nearest frequency
        freq_idx = np.argmin(approx_freq - (1000 / freqs))
        midx_idx = (dset.shape[0] - 2 * skip_x) // 2
        midy_idx = (dset.shape[1] - 2 * skip_y) // 2

        if direction == "y":
            ax.plot(dset[midx_idx, skip_y:-skip_y, freq_idx], label=label)

        else:
            ax.plot(dset[skip_x:-skip_x, midy_idx, freq_idx], label=label)

        ax.set_xlabel(f"Coordinate {direction} / px")
        ax.set_ylabel("$|\\vec{E}$/$\\vec{E_0}|^2$")

    # visualize data
    fig.tight_layout()

    fig.savefig(inputpath.parent / "MultiPlot.svg")
    fig.savefig(inputpath.parent / "MultiPlot.png", dpi=300)

    plt.show()


if __name__ == "__main__":
    main()
