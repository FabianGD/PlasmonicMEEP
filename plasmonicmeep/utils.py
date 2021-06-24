"""
Utils for Plasmonic meep calculations

S.D.G
"""

from pathlib import Path

import h5py
import numpy as np


def vec3_to_nparray(vec):
    """
    Utility to convert a Vector3 to a (3)-shaped np.ndarray
    """
    return np.asarray([vec.x, vec.y, vec.z])


def append_attrs(output, prefix, dset, **kwargs):
    """
    Custom function to append attributes to the h5py data files.
    This function has to be called from the MPI master!
    """

    file = Path(output) / f"{prefix}-{dset}.h5"

    try:
        with h5py.File(file, "a") as h5f:
            for key, val in kwargs.items():
                h5f.attrs[key] = val
                print(f"Added custom attribute {key} with value {val} to {h5f!r}")

    except Exception as e:
        raise Exception(f"Errored with msg: {e}") from e
