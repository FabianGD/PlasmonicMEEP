"""
Plotting utilities
"""

from dataclasses import dataclass
from functools import partial
import logging
from pathlib import Path
from typing import Optional, Tuple, Union, Iterable

import matplotlib as mpl

import h5py
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scipy.signal as ssi

from matplotlib.colors import LogNorm

from .utils import PlasmonicMEEPInputError


def single_plot(
    data: Union[npt.ArrayLike, h5py.Dataset],
    freq: float,
    extent: Tuple[float, float, float, float],
    cmap: str,
    vmax: float,
    extensions: Iterable[str],
    specfolder: Path,
    pattern: str,
    **fig_kwargs,
):
    """Plot and write single maps from the array-like data.

    Args:
        data (Union[npt.ArrayLike, h5py.Dataset]): Data to plot.
            2D array-like structure.

        freq (float): Frequency to use. Will be rendered into the suptitle.

        extent (Tuple[float, float, float, float]): 4-Tuple describing the extent
            of the imshow.

        cmap (str): Matplotlib colormap for the `imshow` and the colorbar.

        vmax (float): Value that defines the maximum of the colormap.

        extensions (Iterable[str]): List of file extensions.
            Needs to include the dot (".").

        specfolder (Path): Folder in which to plot will be saved.

        pattern (str): File name pattern.
            Parameter `wvl` will be rendered with the wavelength in nanometer.
    """

    # Calculate the wavelength
    wvl = 1000 / freq

    # Set up the figure and the axis
    fig, ax = plt.subplots(constrained_layout=True, **fig_kwargs)

    # Do the plotting
    im = ax.imshow(data, extent=extent, cmap=cmap, norm=LogNorm(vmin=1e0, vmax=vmax))
    fig.suptitle("$\\lambda = {wvl:.0f}\\,$nm".format(wvl=wvl))
    fig.colorbar(im, ax=ax)

    for extension in extensions:
        filename = (specfolder / pattern.format(wvl=wvl)).with_suffix(extension)
        fig.savefig(filename)

    plt.close(fig)


def get_extent(
    shape: Iterable[int], resolution: Optional[int] = None
) -> Optional[Tuple[float, float, float, float]]:

    """# TODO

    Returns:
        [type]: [description]
    """

    # Calculate extent, if required
    if resolution:
        max_x = shape[0] / (resolution * 1e-3)
        max_y = shape[1] / (resolution * 1e-3)
        extent = (-max_x, max_x, -max_y, max_y)
    else:
        extent = None

    return extent


def label_ax_coordinates(
    ax: mpl.Axes, extent: Optional[Tuple[float, float, float, float]] = None
) -> None:
    """
    # TODO

    Args:
        ax (mpl.Axes): [description]
        extent (Optional[Tuple[float, float, float, float]], optional): [description]. Defaults to None.
    """

    if extent:
        label = "Coordinate {coord} / nm"
    else:
        label = "Coordinate {coord} / px"

    ax.set_xlabel(label.format(coord="x"))
    ax.set_ylabel(label.format(coord="y"))


def multiplot_enhancement(
    map_data: Union[npt.ArrayLike, h5py.Dataset],
    freqs: Iterable[float],
    folder: Path,
    resolution: Optional[int] = None,
    pattern: str = "map_{wvl:.0f}nm",
    subfolder: Optional[str] = None,
    extensions: Iterable[str] = [".png"],
    cmap: str = "viridis",
    vmax: Optional[Optional[float]] = None,
    **fig_kwargs,
) -> None:
    """Plot the enhancement data for the map_data/freqs iterators.

    Args:
        map_data (Union[npt.ArrayLike, h5py.Dataset]): 3D arraylike structure.
            Iterated data has to be on the last dimension with length `n`.

        freqs (Iterable[float]): Frequencies, usually given in `µm`.
            Length should be `n` (length of the last dim. of map_data).

        folder (Path): Data folder used for writing the plots.

        resolution (Optional[int], optional): Resolution of the simulation.
            Defaults to None.

        pattern (str, optional): Filename pattern.
            `wvl` will be substituted for the wavelength in nm units.
            Defaults to "map_{wvl:.0f}nm".

        subfolder (Optional[str], optional): Whether the plots shall be put
            in a separate folder. Defaults to None.

        extensions (Iterable[str], optional): File name extensions.
            Has to include the dot. Defaults to [".png"].

        cmap (str, optional): Matplotlib colormap. Defaults to "viridis".

        vmax (Optional[Optional[float]], optional): Maximum value of the colormap.
            Defaults to None (each plot scales individually).
    """

    # Create the subfolder if wanted
    if subfolder:
        specfolder = folder / subfolder
        specfolder.mkdir(exist_ok=True)
    else:
        specfolder = folder

    print(f"Plotting {map_data.shape[-1]} maps.")

    func = partial(
        single_plot,
        extent=get_extent(map_data.shape, resolution),
        cmap=cmap,
        vmax=vmax,
        extensions=extensions,
        specfolder=specfolder,
        pattern=pattern,
        **fig_kwargs,
    )

    # Plot the frames. Parallel does not work satisfactorily for large files.
    for args in zip(map_data.transpose(-1, 0, 1), freqs):
        func(*args)


def find_data(folder: Path) -> Tuple[Path, Path]:

    """# TODO

    Raises:
        PlasmonicMEEPInputError: [description]

    Returns:
        [type]: [description]
    """

    norm = list(folder.glob("*-norm.h5"))
    ref = list(folder.glob("*-ref.h5"))

    if len(norm) != 1 or len(ref) != 1:
        raise PlasmonicMEEPInputError(
            "Either no or more than one 'norm' or 'ref' file(s) was not (were) "
            "found in the input folder. Exiting."
        )

    return norm[0], ref[0]


def inner(a: npt.ArrayLike, b: npt.ArrayLike) -> npt.ArrayLike:
    """# TODO"""
    return np.einsum("...a,...a->...", a, b)


def get_phasedata(
    file: Path, polarization: str = "ex", xy_skip: int = 30
) -> npt.ArrayLike:

    """# TODO

    Returns:
        [type]: [description]
    """

    with h5py.File(file, "r") as f:

        print("Opened the file.")

        keys = list(f.keys())
        attrs = dict(f.attrs)

        print("File attributes: {!r}".format(attrs))

        dshape = f[keys[0]].shape
        print("Dataset shape: ", dshape)

        # Assert square dataset
        shape = dshape[0]
        assert shape == dshape[1]

        # Read the real part of the dataset
        try:
            real_ds = f[polarization + ".r"]
        except KeyError:
            real_ds = f[polarization]

        # Getting the real array and transposing it
        # old data is stored 90 deg. rotated
        slice_real = np.array(real_ds[xy_skip:-xy_skip, xy_skip:-xy_skip, :]).transpose(
            (1, 0, 2)
        )

        # Calculate the analytical signal using Hilbert transform
        analytical_signal = ssi.hilbert(slice_real, axis=-1)

        del slice_real

    return analytical_signal


@dataclass
class PhaseData:
    """# TODO

    Returns:
        [type]: [description]
    """

    norm: npt.ArrayLike
    ref: npt.ArrayLike

    _phase: Optional[npt.ArrayLike] = None

    def _calc_phase(self, arctan: bool = True) -> npt.ArrayLike:
        phaseH = inner(self.ref, np.conj(self.norm)) / np.sqrt(
            inner(self.ref, np.conj(self.ref)) * inner(self.norm, np.conj(self.norm))
        )
        if arctan:
            phase = np.arctan(phaseH.imag / phaseH.real)
        else:
            phase = np.angle(phaseH)

        return phase

    @property
    def phase(self) -> npt.ArrayLike:
        if not self._phase:
            self._phase = self._calc_phase()

        return self._phase


def plot_phasemap(
    inputfile: Path,
    skip_xy: int,
    polarization: str = "ey",
    resolution: Optional[int] = None,
    **fig_kwargs,
) -> mpl.Figure:
    """
    Plot a phasemap by finding the field data files from the specified file,
    calculating and plotting the relative phase using Hilbert transform and
    returning the matplotlib figure.

    Args:
        inputfile (Path): Path to a file in the same folder as the field files.

        skip_x (int): Number of pixels to skip on either side of both dimensions.

        polarization (str, optional): Polarization component to plot the phase from.
            Defaults to "ey".

        resolution (int, optional): Resolution used in the FDTD simulation.
            Defaults to None.

        **fig_kwargs: Optional arguments given to `plt.subplots`.

    Returns:
        mpl.Figure: Figure in which the phase map is plotted.
    """

    normfile, reffile = find_data(inputfile.parent)
    kwargs = dict(xy_skip=skip_xy, polarization=polarization)
    phase_data = PhaseData(
        get_phasedata(normfile, **kwargs), get_phasedata(reffile, **kwargs)
    )

    extent = get_extent(phase_data.phase.shape, resolution=resolution)

    phasefig, phaseax = plt.subplots(**fig_kwargs)
    im = phaseax.imshow(
        phase_data.phase,
        cmap="twilight_shifted",
        extent=extent,
        vmin=-np.pi / 2,
        vmax=np.pi / 2,
        interpolation="nearest",
    )
    phasefig.colorbar(im, ax=phaseax)
    label_ax_coordinates(phaseax, extent)

    return phasefig


def open_h5file(h5file: Path) -> h5py.Dataset:
    """Open the first found hdf5 dataset found in the file. Useful for single-dataset files like the enhancement data files.

    Args:
        h5file (Path): Path to the HDF5 data file.

    Returns:
        h5py.Dataset: The first dataset found in the file.
    """

    logging.info("Opening '{!s}' as data file...".format(h5file))

    file = h5py.File(h5file, "r")
    keys = list(file.keys())
    dset = [file[key] for key in keys]
    # There should not be any more datasets
    dset = dset[0]

    return dset
