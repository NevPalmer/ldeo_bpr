"""Test for arg_parser.py."""

from pathlib import Path

import numpy as np
import pytest

from ldeo_bpr import parse_arguments


@pytest.mark.parametrize(
    "arg_name, arg_in, arg_out",
    (
        # Tuples of CLI argument name with valid input arg value and their
        # resulting arg value.
        ("apgini", "./ParosAPG.ini", Path("./ParosAPG.ini")),
        ("loggerini", "./APGlogger.ini", Path("./APGlogger.ini")),
        ("tempsmth", "5001", 5001),
        ("decimate", "5", 5),
        ("clkstart", "2018-10-12_07:11:30", np.datetime64("2018-10-12T07:11:30")),
        ("beginwndw", "2019-04-29_03:00:00.0", np.datetime64("2019-04-29T03:00:00.0")),
        ("endwndw", "2019-04-29_04:00:00.0", np.datetime64("2019-04-29T04:00:00.0")),
        ("bininterval", "24H", np.timedelta64(86400000, "ms")),
        ("gpssynctime", "2019-307_19:49:30", np.datetime64("2019-11-03T19:49:30")),
        ("synctickcount", "0x03CBB30B84", 16302410628),
        ("plot", "rd", {"format": "r", "output": "d"}),
        ("fmttime", "d", "d"),
        ("outfile", "./out/GNS22-PH.csv", Path("./out/GNS22-PH.csv")),
        ("network", "XX", "XX"),
        ("network", "Ab", "AB"),
        ("network", "d3", "D3"),
        ("station", "sssss", "SSSSS"),
    ),
)
def test_valid_cli_args(arg_name, arg_in, arg_out):
    """Test when valid CLI args provided."""
    args = (
        f"--infile . --snapg 999999 --version CSAC2013 " f"--{arg_name} {arg_in}"
    ).split()
    result_args = parse_arguments(args)
    print(result_args)
    assert getattr(result_args, arg_name) == arg_out


def test_valid_miniseed_cli_args():
    """Test when valid miniseed CLI args provided."""
    args = (
        "--infile . --snapg 999999 --version CSAC2013 "
        "--mseedpath ./mseed/ --network XX --station SSSSS"
    ).split()
    result_args = parse_arguments(args)
    print(result_args)
    assert result_args.mseedpath == Path("./mseed/")


@pytest.mark.parametrize(
    "arg_name, arg_in",
    (
        # Tuples of CLI argument name and invalid input arg value.
        ("bininterval", "24F"),
        ("bininterval", "XH"),
        ("network", "xxx"),
        ("network", "ABC"),
        ("network", "A$"),
        ("network", "G"),
        ("network", ""),
        ("station", "ssss"),
        ("station", "AAAAAA"),
        ("station", "XXX$#"),
        ("mseedpath", "./mseed/"),
    ),
)
def test_invalid_cli_args(arg_name, arg_in):
    """Test when invalid CLI args provided."""
    args = (
        f"--infile . --snapg 999999 --version CSAC2013 --{arg_name} {arg_in}"
    ).split()
    with pytest.raises(SystemExit):
        (parse_arguments(args))


@pytest.mark.parametrize(
    "args_str",
    (
        "--infile . --snapg 999999",
        "--infile . --version CSAC2013",
        "--snapg 999999 --version CSAC2013",
        "",
    ),
)
def test_reqd_cli_args_missing(args_str):
    """Test when required CLI args are missing."""
    args = args_str.split()
    with pytest.raises(SystemExit):
        (parse_arguments(args))


@pytest.mark.parametrize(
    "arg_name, arg_out",
    (
        # Tuples of CLI argument name and default output arg value.
        ("apgini", Path("./ParosAPG.ini")),
        ("loggerini", Path("./APGlogger.ini")),
        ("tempsmth", 1),
        ("decimate", 0),
        ("clkstart", np.datetime64("2000-01-01T00:00:00")),
        ("beginwndw", None),
        ("endwndw", None),
        ("bininterval", np.timedelta64(0, "ms")),
        ("gpssynctime", None),
        ("synctickcount", None),
        ("plot", {"format": "n", "output": ""}),
        ("fmttime", "s"),
        ("outfile", None),
        ("network", None),
        ("station", None),
        ("mseedpath", None),
    ),
)
def test_default_cli_args(arg_name, arg_out):
    """Test defaults when only required CLI args are provided."""
    args = ("--infile . --snapg 999999 --version CSAC2013 ").split()
    result_args = parse_arguments(args)
    print(result_args)
    assert getattr(result_args, arg_name) == arg_out


if __name__ == "__main__":
    pytest.main()
