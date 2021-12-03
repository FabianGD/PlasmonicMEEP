"""
Plotting utilities
"""

from functools import partial
from pathlib import Path
from typing import Optional, Tuple, Union, Iterable

import h5py
import matplotlib.pyplot as plt
import numpy.typing as npt
from matplotlib.colors import LogNorm


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

    if extent:
        label = "Coordinate {coord} / nm"
    else:
        label = "Coordinate {coord} / px"

    ax.set_xlabel(label.format(coord="x"))
    ax.set_ylabel(label.format(coord="y"))

    for extension in extensions:
        filename = (specfolder / pattern.format(wvl=wvl)).with_suffix(extension)
        fig.savefig(filename)

    plt.close(fig)


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
):
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

    # Calculate extent, if required
    if resolution:
        max_x = map_data.shape[0] / (resolution * 1e-3)
        max_y = map_data.shape[1] / (resolution * 1e-3)
        extent = (-max_x, max_x, -max_y, max_y)
    else:
        extent = None

    func = partial(
        single_plot,
        extent=extent,
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
