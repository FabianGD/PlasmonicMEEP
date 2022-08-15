"""
Visualise the FDTD results.
"""

import argparse
import logging
import sys
from pathlib import Path

import matplotlib as mpl

import matplotlib.pyplot as plt
import numpy as np
import tqdm

from .plotting import (
    multiplot_enhancement,
    open_h5file,
    plot_phasemap,
    read_h5ds_direct,
)
from .multiviewer import multi_slice_viewer


def argparsing(*args):
    """
    Module-level argparsing
    """
    parser = argparse.ArgumentParser(
        description="Visualize calculated field enhancement"
    )
    parser.add_argument("file", type=str, help="File to use")
    parser.add_argument(
        "-k", "--xyskip", type=int, help="Number of pixels to skip", default=30
    )
    parser.add_argument(
        "-p", "--mplstyle", type=str, help="Give a matplotlib style to use. Optional."
    )
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        default=False,
        help="Get an interactive plot of the maps to scroll through.",
    )
    parser.add_argument(
        "-f",
        "--min-freq",
        type=float,
        default=0.5,
        help=(
            "Minimum frequency (in 1/µm) to plot. "
            "Default is 0.5 which is equivalent to 2000 nm wavelength."
        ),
    )
    parser.add_argument(
        "--every-nth",
        type=int,
        default=1,
        dest="every",
        help="Number of grid points to skip.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help=(
            "Optionally, specify an output directory. Defaults to the folder in which the "
            "input file is located."
        ),
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        help=(
            "Optional. Resolution of the simulation. This allows rendering the "
            "maps with nanometer scales instead of pixels."
        ),
    )
    parser.add_argument(
        "-c",
        "--colormap",
        type=str,
        default="viridis",
        help="Matplotlib colormap to use for plotting the maps.",
    )
    parser.add_argument(
        "--disable-phase",
        action="store_true",
        help=(
            "Disable phase map generation. Phase maps are based on Hilbert transform of "
            "the real-valued signal. This relies on the *-norm.h5 and *-ref.h5 data "
            "to be present in the same folder as the enhancement file."
        ),
    )
    parser.add_argument(
        "--disable-maps",
        action="store_true",
        dest="disable_maps",
        help="Whether to disable intensity map generation altogether.",
    )
    parser.add_argument(
        "--disable-spectrum",
        action="store_true",
        dest="disable_spectrum",
        help="Whether to disable enhancement spectrum generation altogether.",
    )

    args = parser.parse_args(*args)

    return args


def enhancement_spectrum(
    ax,
    dset,
    freqs: np.ndarray,
    slice_xy: slice,
    plot_slice: slice = slice(-1),
    percent: float = 99.9,
):

    enhancement = np.zeros_like(freqs)
    for freq_idx in tqdm.trange(len(freqs), desc="Enhancement"):

        tmp_data = read_h5ds_direct(
            dset, slice_selector=np.s_[slice_xy, slice_xy, freq_idx]
        )

        # taking 99.9 percentile instead of maximum to eliminate
        # hotspot enhancement
        enhancement[freq_idx] = np.percentile(tmp_data, percent)
        del tmp_data

    (line, ) = ax.plot(1000 / freqs[plot_slice], enhancement[plot_slice])

    return line


def main(*args):
    """
    Main function for the visualiser
    """

    args = argparsing(*args)
    inputfile = Path(args.file).expanduser().resolve()

    if not args.output_dir:
        output_dir = inputfile.parent
    else:
        output_dir: Path = args.output_dir.expanduser().resolve()
        output_dir.mkdir(exist_ok=True, parents=True)

    if not args.interactive:
        mpl.use("agg")

    if args.mplstyle is not None:
        mplstyle = Path(args.mplstyle)
        plt.style.use(mplstyle)

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # try:
    #     inputfile.parent.glob(".yml")

    # Skip N pixels from the sides to eliminate enhancement on boundaries
    slice_xy = slice(args.xyskip, -args.xyskip, args.every)

    map_figsize = (6, 5)

    # Prepare frequency h5 file
    h5file, dset = open_h5file(inputfile)

    # Get the frequency info
    freqs = dset.attrs["freqs"]

    # skip the nonphysical low-frequncy part of spectra
    skip_freq = np.argmin(np.abs(args.min_freq - freqs))

    s_freqs = freqs[skip_freq:]

    # Calculate enhancement at each frequency

    if not args.disable_spectrum:
        logger.info("Calculating enhancement over frequency.")
        specfig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

        enhancement_spectrum(
            ax, dset, freqs, slice_xy=slice_xy, plot_slice=np.s_[skip_freq:]
        )

        ax.set_xlabel("Wavelength / nm")
        ax.set_ylabel("$|\\vec{E}$/$\\vec{E_0}|^2$")

        logger.info("Done calculating enhancement over frequency.")

        # Print spectrum file
        specfile = output_dir / "Spectrum_Enhancement_over_Wavelength.png"
        specfig.savefig(specfile)
        specfig.savefig(specfile.with_suffix(".svg"))

    if not args.disable_maps:
        # Plot the field enhancement maps
        data = read_h5ds_direct(
            dset, slice_selector=np.s_[slice_xy, slice_xy, skip_freq:]
        ).transpose(1, 0, -1)

        multiplot_enhancement(
            data,
            freqs=s_freqs,
            folder=output_dir,
            subfolder="maps",
            extensions=[".png", ".svg"],
            cmap=args.colormap,
            figsize=map_figsize,
            resolution=args.resolution,
        )

    if not args.disable_phase:
        for polarization in ["ex", "ey"]:
            phasefig = plot_phasemap(
                inputfile,
                slice_xy=slice_xy,
                polarization=polarization,
                figsize=map_figsize,
            )

            phase_fname = (
                output_dir / f"Map_Phase_Hilbert_{polarization.capitalize()}.png"
            )
            phasefig.savefig(phase_fname)
            phasefig.savefig(phase_fname.with_suffix(".svg"))

    if args.interactive:
        # Visualize data with the multiviewer
        data = read_h5ds_direct(
            dset, slice_selector=np.s_[slice_xy, slice_xy, skip_freq:]
        ).transpose(1, 0, -1)

        multi_slice_viewer(
            volume=data,
            index_function=lambda x: 1000 / freqs[skip_freq:][x],
        )
        plt.show()

    h5file.close()


if __name__ == "__main__":
    main()
