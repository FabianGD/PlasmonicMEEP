"""Unit tests for the

Returns:
    [type]: [description]
"""

import shutil
import tempfile
from collections import namedtuple
from pathlib import Path

import h5py
import meep as mp
import numpy as np
import pytest

from plasmonicmeep.utils import (  # pylint: disable=import-error
    append_attrs,
    append_timegrid,
    meep_path,
    vec3_to_nparray,
)


_H5File = namedtuple("h5file", ["path", "prefix", "dset", "timestep", "n_timesteps"])


@pytest.mark.parametrize(
    "folder,prefix,dset,expected",
    [
        ("/tmp", "test", "test-center", Path("/tmp/test-test-center.h5")),
        ("/tmp", "", "test-center", Path("/tmp/test-center.h5")),
    ],
)
def test_meep_path(folder, prefix, dset, expected):

    got = meep_path(folder, prefix, dset)
    assert got == expected


@pytest.fixture(scope="function", name="h5f")
def h5file() -> _H5File:
    """
    Fixture that copies a workable meep h5 file to a temporary folder and returns a H5File NT
    """
    orig_file = Path(__file__).parent / "data/plas-meep-norm-center_no-tgrid.h5"
    tmpdir = Path(tempfile.mkdtemp())
    newpath = tmpdir / orig_file.name

    shutil.copy(orig_file, newpath)

    # print(f"h5 file copied to {newpath!r}")

    f = _H5File(
        path=newpath,
        prefix="plas-meep",
        dset="norm-center_no-tgrid",
        timestep=0.025,
        n_timesteps=8000,
    )

    return f


def test_append_tgrid(h5f: _H5File) -> None:  #
    """
    Test appending a time grid to a h5f file.

    Args:
        h5f (_H5File) h5file namedtuple, getting that from the fixture.
    """
    pref = h5f.prefix
    dset = h5f.dset
    folder = h5f.path.parent

    with h5py.File(h5f.path, "r") as h5:
        dss = list(h5.keys())

    append_timegrid(
        folder=folder,
        prefix=pref,
        dset=dset,
        n_timesteps=h5f.n_timesteps,
        timestep=h5f.timestep,
    )

    with h5py.File(h5f.path, "r") as h5:
        new_dss = list(h5.keys())
        new_n_dts = [h5[i].shape[-1] for i in new_dss]

    assert len(new_dss) == len(dss) + 1
    assert new_n_dts.count(new_n_dts[0]) == len(new_dss)  # All have the same length
    assert new_n_dts[0] == h5f.n_timesteps + 1


def test_append_attrs(h5f: _H5File) -> None:
    """
    Test appending attrs to a h5f file.

    Args:
        h5f (_H5File) h5file namedtuple, getting that from the fixture.
    """

    pref = h5f.prefix
    dset = h5f.dset
    folder = h5f.path.parent

    with h5py.File(h5f.path, "r") as h5:
        attrs = dict(h5.attrs)

    to_append = {
        "test_int": 9999,
        "test_str": "Hello",
        "test_float": 99.995,
        "test_array": np.array([1]),
    }

    append_attrs(folder=folder, prefix=pref, dset=dset, **to_append)

    with h5py.File(h5f.path, "r") as h5:
        new_attrs = dict(h5.attrs)

    assert {**attrs, **to_append} == new_attrs

    if attrs == dict():
        assert to_append == new_attrs


@pytest.mark.parametrize(
    "vector, expected",
    [
        (mp.Vector3(), np.array([0, 0, 0])),
        (mp.Vector3(1), np.array([1, 0, 0])),
        (mp.Vector3(1, 1), np.array([1, 1, 0])),
        (mp.Vector3(1, 1, 1), np.array([1, 1, 1])),
        (mp.Vector3(1, -10, 1), np.array([1, -10, 1])),
    ],
)
def test_vec3_to_nparray(vector: mp.Vector3, expected: np.ndarray):
    """Unit test the vec3_to_ndarray function.

    Args:
        vector (mp.Vector3): Input meep vector.
        expected (np.ndarray): Expected result.
    """

    assert np.all(vec3_to_nparray(vector) == expected)
