"""
Calculate the power spectrum in chunks.
"""

import argparse
import logging

import h5py
import numpy as np
from joblib import Parallel, delayed


def argparsing():
    """
    Build an argument parser and parse arguments.
    """
    parser = argparse.ArgumentParser(
        description="Perform FFT on time-domain EM data to get EM enhancement"
    )
    parser.add_argument("file", type=str, help="Field file to use")
    parser.add_argument("reffile", type=str, help="Reference file to use")
    parser.add_argument("savefile", type=str, help="File to save results")
    parser.add_argument("-s", "--size", type=int, help="Batch size")
    parser.add_argument("-f", "--freq", type=float, help="Central gaussian frequency")
    parser.add_argument("-w", "--width", type=float, help="Width of gaussian frequency")
    parser.add_argument(
        "--loglevel",
        type=str,
        choices=["DEBUG", "INFO", "WARNING"],
        help="Set the logging level. One of ['DEBUG', 'INFO', 'WARNING']",
        default="INFO",
    )
    args = parser.parse_args()

    return args


def get_power(data, interp=4, nproc=4, parallel=True):
    """
    Convert a chunk of data from time domain to squared frequency domain (power spectrum)
    """

    timelen = data.shape[2]

    # length of data in frequency domain
    freqlen = (timelen * interp) // 2

    if parallel:
        data = np.array_split(data, nproc, axis=1)
        fdata = Parallel(n_jobs=nproc, prefer="threads")(
            delayed(np.fft.fft)(d, n=freqlen * 2, axis=-1) for d in data
        )
        fdata = np.concatenate(fdata, axis=1) / timelen

    else:
        fdata = np.fft.fft(data, n=freqlen * 2, axis=-1) / timelen

    fdata = 2 * np.abs(fdata[:, :, 0:freqlen])
    power = np.square(fdata)

    del data

    return power


def main():
    """Main function. To be enhanced."""

    # Some argparsing.
    args = argparsing()
    cfreq = args.freq
    freq_width = args.width

    logging.basicConfig(level=args.loglevel)

    # Read all the data from the data and reference files.
    logging.info("Opening {} as data file...".format(args.file))
    dfile = h5py.File(args.file, "r", rdcc_nbytes=500 * 1024 ** 2)
    # read all datasets in the hf5 file
    keys = list(dfile.keys())
    data = [dfile[key] for key in keys]

    # read all datasets in the hf5 file
    rfile = h5py.File(args.reffile, "r", rdcc_nbytes=500 * 1024 ** 2)
    reference = [rfile[key] for key in keys]

    print(reference)

    logging.info("Reference shape: {}".format(reference[0].shape))
    logging.info("Components: {}".format(keys))
    logging.info("Data shape: {}".format(data[0].shape))

    # The lengths for later use
    x_length, y_length, t_length = data[0].shape

    # Parameter calculation and defaults
    # optional zero padding of data to increase resolution, which is equivalent to interpolation
    interp = 4
    freq_length = (t_length * interp) // 2  # length of data in frequency
    period = 0.5 / (cfreq + freq_width * 0.5)  # period of maximal frequncy component
    freq_max = 0.5 / period  # maxumal frequnce in spectrum (Nyquist-Shannon)
    freqs = np.linspace(0, 1, freq_length) * freq_max

    # domain is 2 times smaller due to FFT symmetry
    logging.info(
        "Temporal length of data: {}, "
        "length of freq domain: {}, "
        "Maximum frequency: {}".format(t_length, freq_length, freq_max)
    )

    # Open new output h5 file
    filename = args.savefile
    outfile = h5py.File(filename, "w", rdcc_nbytes=500 * 1024 ** 2)
    dset = outfile.create_dataset(
        "enhancement", (x_length, y_length, freq_length), dtype="f"
    )
    dset.attrs["freqs"] = freqs

    # Batch size in the x direction
    batch_size = min([getattr(args, "size", 50), x_length])

    # iterate over chunks, calculating power spectra and saving to file
    for i in range(0, x_length, batch_size):

        logging.info("Iteration {} from {}".format(i, x_length))
        bshape = data[0][i : i + batch_size, ...].shape

        enh_data = np.zeros((bshape[0], bshape[1], freq_length))
        ref_data = np.zeros((bshape[0], bshape[1], freq_length))

        for component_e, component_r in zip(data, reference):
            logging.info("Processing component {} of data.".format(component_e.name))
            logging.info("Processing component {} of ref.".format(component_r.name))

            enh_data += get_power(component_e[i : i + batch_size, ...])
            ref_data += get_power(component_r[i : i + batch_size, ...])

        dset[i : i + batch_size, ...] = enh_data[...] / ref_data[...]

    # Close all the files.
    dfile.close()
    rfile.close()
    outfile.close()


if __name__ == "__main__":
    main()
