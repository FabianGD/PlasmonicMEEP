"""
Visualise the FDTD results.
"""

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .multiviewer import multi_slice_viewer


def argparsing():
    """
    Module-level argparsing
    """
    parser = argparse.ArgumentParser(
        description="Visualize calculated field enhancement"
    )
    parser.add_argument("file", type=str, help="File to use")
    parser.add_argument(
        "-s", "--save", action="store_true", help="Save calculated spectra to file"
    )
    parser.add_argument(
        "-x", "--excel", action="store_true", help="Save to excel insted of CSV"
    )
    parser.add_argument("-k", "--xyskip", type=int, help="Number of pixels to skip")
    parser.add_argument(
        "-p", "--mplstyle", type=str, help="Give a matplotlib style to use. Optional."
    )
    parser.add_argument(
        "-r", "--resolution", type=int, help="Resolution of the simulation, in px / µm."
    )

    args = parser.parse_args()

    return args


def main():
    """
    Main function for the visualiser
    """
    args = argparsing()

    if args.mplstyle is not None:
        mplstyle = Path(args.mplstyle)
        plt.style.use(mplstyle)

    print("Opening {} as data file...".format(args.file))

    # Prepare frequency h5 file
    file = h5py.File(args.file, "r")
    keys = list(file.keys())
    dset = [file[key] for key in keys]
    # There should not be any more datasets
    dset = dset[0]

    # Get the frequency info
    freqs = dset.attrs["freqs"]
    freqlen = freqs.shape[0]

    enhancement = np.zeros_like(freqs)
    data = dset[:]

    # Skip N pixels from the sides to eliminate enhancement on boundaries
    if args.xyskip:
        skip_x = args.xyskip
        skip_y = args.xyskip
    else:
        skip_x = 10
        skip_y = 0

    for freq_idx in range(freqlen):
        tmp_data = data[skip_x:-skip_x, skip_y:-skip_y, freq_idx]
        # taking 99.9 percentile instead of maximum to eliminate
        # hotspot enhancement
        enhancement[freq_idx] = np.percentile(tmp_data, 99.9)

    # skip the nonphysical low-frequncy part of spectra
    skip_freq = np.argmin(np.abs(0.5 - freqs))

    # Save the spectrum to file
    if args.save:
        df = pd.DataFrame()
        df["lambda"] = 1000 / freqs[skip_freq:]
        df["enhancement"] = enhancement[skip_freq:]

        fpath = Path(args.file).parent
        fname = Path(args.file).stem

        if args.excel:
            df.to_excel(fpath / "".join(["spectra_", str(fname), ".xlsx"]), index=False)
        else:
            df.to_csv(fpath / "".join(["spectra_", str(fname), ".csv"]), index=False)

        sys.exit()

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(1000 / freqs[skip_freq:], enhancement[skip_freq:])
    ax.set_xlabel("Wavelength / nm")
    ax.set_ylabel("$|\\vec{E}$/$\\vec{E_0}|^2$")

    fig.tight_layout()
    plt.show()

    print(args.resolution)

    data = data.transpose(1, 0, -1)

    # visualize data
    multi_slice_viewer(
        data[skip_x:-skip_x, skip_y:-skip_y, :],
        index_function=lambda x: 1000 / freqs[x],
        resolution=args.resolution
    )

    plt.show()


if __name__ == "__main__":
    main()
