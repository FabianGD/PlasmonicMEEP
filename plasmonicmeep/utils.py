"""
Utils for Plasmonic meep calculations

S.D.G
"""

import logging
from pathlib import Path
from typing import Union

import h5py
import numpy as np


def meep_path(folder: Union[Path, str], prefix: str, dset: str):
    """Get the Path to a MEEP dataset file

    Args:
        folder (Union[Path, str]): Output folder. Set via sim.use_output_directory()
        prefix (str): Prefix. Get it via sim.get_filename_prefix()
        dset (str): Data set descriptor used in mp.to_appended step function.

    Returns:
        pathlib.Path: Path to the h5 file.
    """

    if prefix != "":
        pth = Path(folder) / f"{prefix}-{dset}.h5"
    else:
        pth = Path(folder) / f"{dset}.h5"

    return pth


def vec3_to_nparray(vec):
    """
    Utility to convert a Vector3 to a (3)-shaped np.ndarray
    """
    return np.asarray([vec.x, vec.y, vec.z])


def append_attrs(folder: Union[Path, str], prefix: str, dset: str, **kwargs) -> None:
    """
    Custom function to append attributes to the h5py data files.
    This function has to be called from the MPI master!
    """

    try:
        with h5py.File(meep_path(folder, prefix, dset), "a") as h5f:
            for key, val in kwargs.items():
                h5f.attrs[key] = val
                logging.info(
                    "Added custom attribute {} with value {} to {}".format(
                        key, val, repr(h5f)
                    )
                )

    except Exception as e:
        raise Exception(f"Errored with msg: {e}") from e


def append_timegrid(
    folder: Union[Path, str], prefix: str, dset: str, n_timesteps: float, timestep
) -> None:
    """
    Custom function to append an approximated time grid to the h5py data files.
    This function has to be called from the MPI master!
    """

    shapes = []

    with h5py.File(meep_path(folder, prefix, dset), "a") as h5f:

        for dataset in h5f.values():
            shapes.append(dataset.shape[-1])

        logging.info("Fields have time domain lengths: {}".format(repr(shapes)))

        if shapes.count(shapes[0]) == len(shapes):
            t_arr = np.arange(0, n_timesteps + 1) * timestep
            h5f.create_dataset("tgrid", data=t_arr)

            logging.info(
                "Calculated time grid: {} ... {} with {} steps.".format(
                    t_arr[0], t_arr[-1], len(t_arr)
                )
            )

        else:
            logging.warning("Time domain step numbers in the dataset are not the same.")

    return None
