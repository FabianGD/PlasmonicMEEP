"""
Testing functions for the argparsing and functional tests for PlasmonicMEEP.
"""

import pytest

import matplotlib as mpl
from plasmonicmeep.args import positive_type  # pylint: disable=import-error
from plasmonicmeep.fdtd import main  # pylint: disable=import-error
from plasmonicmeep.utils import PlasmonicMEEPInputError  # pylint: disable=import-error


@pytest.mark.parametrize(
    "value,rtype,other_allowed,expected",
    [
        pytest.param("-10", float, None, None, marks=pytest.mark.xfail),
        pytest.param("-10", int, None, None, marks=pytest.mark.xfail),
        pytest.param("0", int, None, None, marks=pytest.mark.xfail),
        pytest.param("auto", int, None, None, marks=pytest.mark.xfail),
        pytest.param("1.3", int, None, None, marks=pytest.mark.xfail),
        pytest.param("-1", int, [-1], -1),
        pytest.param("-1", int, ["-1"], "-1"),
        pytest.param("10", int, None, 10),
        pytest.param("1.3", float, None, 1.3),
        pytest.param("1E9", float, None, 1e9),
        pytest.param("1000000", int, None, 1_000_000),
        pytest.param("auto", int, ["auto"], "auto"),
    ],
)
def test_positive_type(value, rtype, other_allowed, expected):
    """
    Test cases for the utils.positive_type() function
    """
    rval = positive_type(value, rtype, other_allowed)

    assert rval == expected

    if other_allowed is None:
        other_allowed = []

    if not value in other_allowed and not rval in other_allowed:
        assert isinstance(rval, rtype)


@pytest.mark.parametrize(
    "args",
    [
        pytest.param("-vg".split()),
        pytest.param("-vvvg".split()),
        pytest.param("--disable-single-point".split()),
        pytest.param("-csg".split()),
        pytest.param("-x 1.0 -y 1.0".split()),
        pytest.param("--point 0.1 0.1 --run-until 10 -w 1.1 -f 1.2 -r 200".split()),
        pytest.param(
            "-x -1.0".split(),
            marks=pytest.mark.xfail,
        ),
        pytest.param(
            "--disable-cell-field --disable-single-point".split(),
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_main_c_args(args):
    """Functional tests for the FDTD module."""
    mpl.use("agg")

    try:
        main(argv=args)
    except SystemExit as e:
        if e.code is None or e.code == 0:  # This indicates normal termination
            pass
        elif e.code == 2:  # Exit code 2 is from argument parsing errors.
            raise PlasmonicMEEPInputError() from e
        else:
            raise e
