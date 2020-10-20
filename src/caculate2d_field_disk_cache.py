"""
Calculate the power spectrum in chunks.
"""

import argparse

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

    args = parser.parse_args()

    return args


def get_power(data, interp=4, nproc=4, parallel=True):
    """
    Convert a chunk of data from time domain to squared frequency domain (power spectrum)
    """

    timelen = data.shape[2]

    # length of data in frequency domain
    flen = (timelen * interp) // 2

    if parallel:
        data = np.array_split(data, nproc, axis=1)
        fdata = Parallel(n_jobs=nproc, prefer="threads")(
            delayed(np.fft.fft)(d, flen * 2, axis=-1) for d in data
        )
        fdata = np.concatenate(fdata, axis=1) / timelen

    else:
        fdata = np.fft.fft(data, flen * 2, axis=-1) / timelen

    fdata = 2 * np.abs(fdata[:, :, 0:flen])
    power = np.square(fdata)

    del data

    return power


def main():
    """Main function. To be enhanced."""
    args = argparsing()

    cfreq = args.freq
    freq_width = args.width

    data = h5py.File(args.file, "r", rdcc_nbytes=500 * 1024 ** 2)
    print("Opening {} as data file...".format(args.file))

    keys = list(data.keys())

    # read all datasets in the hf5 file
    data = [data[key] for key in keys]
    print("Components: {}".format(keys))
    print("Data shape: {}".format(data[0].shape))

    rfile = h5py.File(args.reffile, "r", rdcc_nbytes=500 * 1024 ** 2)
    reference = [rfile[key] for key in keys]
    print("Reference shape: {}".format(reference[0].shape))

    timelen = data[0].shape[2]  # length of temporal dimension
    # optional zero padding of data to increase resilution, which is equivalent to interpolation
    interp = 4
    freqlen = (timelen * interp) // 2  # length of data in frequency

    # domain is 2 times smaller due to FFT symmetry
    print("Temporal length of data {}, length of freq domain {}".format(timelen, freqlen))
    period = 0.5 / (cfreq + freq_width * 0.5)  # period of maximal frequncy component
    fmax = 0.5 / period  # maxumal frequnce in spectrum (Nyquist-Shannon)
    print("Freq max {}".format(fmax))

    freqs = np.linspace(0, 1, freqlen) * fmax  # change linspace to arange 0...f_nyq

    filename = args.savefile

    # Length in x & y & t direction
    xytlen = data[0].shape


    # Batch size in the x direction
    batch_size = min([getattr(args, "size", 50), xytlen[0]])

    # Open h5 file
    file = h5py.File(filename, "w", rdcc_nbytes=500 * 1024 ** 2)
    dset = file.create_dataset("enhancement", (*xytlen[:-1], freqlen), dtype="f")
    dset.attrs["freqs"] = freqs

    # iterate over chunks, calculating power spectra and saving to file
    for i in range(0, xytlen[0], batch_size):

        print("Iteration {} from {}".format(i, xytlen[0]))
        bshape = data[0][i : i + batch_size, ...].shape

        enh_data = np.zeros((bshape[0], bshape[1], freqlen))
        ref_data = np.zeros((bshape[0], bshape[1], freqlen))

        for component_e, component_r in zip(data, reference):
            print(component_e.name)
            print(component_r.name)

            enh_data += get_power(component_e[i : i + batch_size, ...])
            ref_data += get_power(component_r[i : i + batch_size, ...])

        dset[i : i + batch_size, ...] = enh_data[...] / ref_data[...]

    file.close()


if __name__ == "__main__":
    main()
