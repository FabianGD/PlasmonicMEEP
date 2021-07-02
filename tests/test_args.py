import pytest

from plasmonicmeep.args import positive_type  # pylint: disable=import-error


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
