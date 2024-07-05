"""Module to parse arguments from command line."""

import datetime as dt
import re
from argparse import ArgumentParser, Namespace
from pathlib import Path

import numpy as np

import ldeo_bpr as bpr


def parse_arguments(args=None) -> Namespace:
    """Parse command line arguments."""
    # Default values and choices for reading params from command line.
    apg_ini = Path("./ParosAPG.ini")
    package_path: Path = Path(__file__).parents[1]
    if not apg_ini.is_file():
        apg_ini: Path = package_path / "ParosAPG.ini"
    logger_ini = Path("./APGlogger.ini")
    if not logger_ini.is_file():
        logger_ini: Path = package_path / "APGlogger.ini"
    logger_versions: list[str] = ["CSAC2013", "Seascan2018", "TEST"]
    clk_start: str = "2000-01-01_00:00:00"  # 'YYYY-MM-DD_hh:mm:ss'
    out_filename: Path | None = None
    mseed_path: Path | None = None
    tmptr_smth_fctr: int = 1
    decmt_intvl: int = 0

    # Read in parameters from command line
    helpdesc = (
        "Reads a raw APG data file and outputs decimated pressure data."
        "Two .ini files are required which must contain configuration values "
        "for the specific Paroscientific pressure transducer and the correct "
        "version of APG logger board used."
    )
    parser = ArgumentParser(description=helpdesc)
    parser.add_argument(
        "-i",
        "--infile",
        help="Full path and filename of raw APG input file.",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-a",
        "--apgini",
        help="Full path and filename for Paros APG "
        f'configuration settings. Default: "{apg_ini}"',
        type=Path,
        default=apg_ini,
    )
    parser.add_argument(
        "-s",
        "--snapg",
        help="Serial number of the Paroscientific APG used. "
        "This must correspond to the serial number of an "
        "entry in the apgini file.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-l",
        "--loggerini",
        help=f"Full path and filename for the APG logger "
        f"board configuration settings. "
        f'Default: "{logger_ini}"',
        type=Path,
        default=logger_ini,
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Specify the version/firmware of the APG logger board used.",
        choices=logger_versions,
        type=str,
        required=True,
    )
    parser.add_argument(
        "-d",
        "--decimate",
        help=f"Required sample interval in seconds for "
        f"pressure decimation. Zero for no decimation. "
        f"Value must equal a single digit integer of seconds "
        f"or minutes or a multiple of 5 or 10."
        f'Default: "{decmt_intvl}"',
        type=int,
        default=decmt_intvl,
    )
    parser.add_argument(
        "-t",
        "--tempsmth",
        help=f"Temperature smoothing factor (must be an odd "
        f"integer). 5001 gives sensible smoothing. 50001 "
        f"gives better smoothing for Seascan logger but is "
        f'slow. Default: "{tmptr_smth_fctr}"',
        type=int,
        default=tmptr_smth_fctr,
    )
    parser.add_argument(
        "-c",
        "--clkstart",
        help=f"Precise date and time when the logger clock "
        f'was started. Format: "YYYY-MM-DDThh:mm:ss" '
        f'Default: "{clk_start}"',
        type=dtstr_to_dt64,
        default=clk_start,
    )
    parser.add_argument(
        "-b",
        "--beginwndw",
        help="Date and time to begin data extraction. "
        "Assumes beginning of file if omitted. "
        'Format: "YYYY-MM-DDThh:mm:ss.s"',
        type=dtstr_to_dt64,
    )
    parser.add_argument(
        "-e",
        "--endwndw",
        help="Date and time to end data extraction. Assumes "
        "end of file if omitted. "
        'Format: "YYYY-MM-DDThh:mm:ss.s"',
        type=dtstr_to_dt64,
    )
    parser.add_argument(
        "-B",
        "--bininterval",
        help="The Bin Interval defines the period of data "
        "that will be processed at each iteration. Each bin "
        "period processed will be appended to the output CSV "
        "file if specified. If specified, multiple MiniSEED "
        "files will be created, one for each Bin extracted."
        'Format: "##[DHM]" where ## is an integer and '
        "character D, H or M indicates Days, Hours or "
        "Minutes.",
        type=intvlstr_to_dt64,
        default=np.timedelta64(0, "ms"),
    )
    parser.add_argument(
        "-g",
        "--gpssynctime",
        help="Precise date and time from GPS clock for "
        "syncronising end time. No clock drift adjustment is "
        'made if omitted. Format: "YYYY-DDD_hh:mm:ss"',
        type=dtstr_to_dt64,
    )
    parser.add_argument(
        "-y",
        "--synctickcount",
        help="The hexidecimal tick count that corresponds to "
        "GPSSYNCTIME. If GPSSYNCTIME is specified and "
        "SYNCTICKCOUNT is omitted, then it is assumed that an "
        "artificial frequency was inserted precisely "
        "at GPSSYNCTIME. This parameter is ignored if "
        'GPSSYNCTIME is omitted. Format: "0xHHHHHHHHHH"',
        type=lambda x: int(x, 0),
    )
    parser.add_argument(
        "-w",
        "--network",
        help=(
            f"Network name to be used in MiniSEED file header. Max "
            f"{bpr.NWK_NAME_LEN} characters."
        ),
        type=nwk_name,
        default=None,
    )
    parser.add_argument(
        "-n",
        "--station",
        help=(
            f"Station name to be used in MiniSEED file header. Max "
            f"{bpr.STN_NAME_LEN} characters."
        ),
        type=stn_name,
        default=None,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Full path and filename for output file. No file "
        "will be generated if not specified.",
        type=Path,
        default=out_filename,
    )
    parser.add_argument(
        "-m",
        "--mseedpath",
        help="Full path of location to save MiniSEED file(s). No file(s) will "
        "be generated if not specified. If specifed, then both"
        "-n/--station and -w/--network are required.",
        type=Path,
        default=mseed_path,
    )
    parser.add_argument(
        "-f",
        "--fmttime",
        help="Specify the format to be used for presenting "
        "time in outputs and plots, to be displayed as either "
        "(s)econds or (d)ate-time.",
        choices=["s", "d"],
        default="s",
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="Generate and display a timeline plot. Either "
        "display only the (f)inal smoothed/decimated result "
        "or additionally dispaly the (r)aw data in the "
        "background or (n)ot generate any plot at all. "
        "Also specify whether to (s)ave as a file, (d)isplay to "
        "the screen or output as (b)oth.",
        choices=["n", "fs", "fd", "fb", "rs", "rd", "rb"],
        default="n",
    )
    parser.add_argument(
        "--noisefilt",
        help="Apply a binned median to pressure & temperature data and remove"
        "any data points that are greater than a predefined range from the"
        "median values.",
        action="store_true",
        default=False,
    )
    args = parser.parse_args(args)

    args.plot = {"format": args.plot[:1], "output": args.plot[1:]}

    if args.mseedpath and not (args.network and args.station):
        parser.error(
            "Arguments -n/--station and -w/--network are both required because "
            "-m/--mseedpath was specified."
        )
    return args


def dtstr_to_dt64(datetime_str: str) -> np.datetime64:
    """Convert a formatted string to a Numpy datetime64."""
    # Standardise the delimiters
    standardised_str = re.sub("[-: _/tT]", "_", datetime_str)
    match1 = re.match(r"^\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2}$", standardised_str)
    match2 = re.match(r"^\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2}\.\d+$", standardised_str)
    match3 = re.match(r"^\d{4}_\d{3}_\d{2}_\d{2}_\d{2}$", standardised_str)
    if match1:
        pydt = dt.datetime.strptime(standardised_str, "%Y_%m_%d_%H_%M_%S")
    elif match2:
        pydt = dt.datetime.strptime(standardised_str, "%Y_%m_%d_%H_%M_%S.%f")
    elif match3:
        pydt = dt.datetime.strptime(standardised_str, "%Y_%j_%H_%M_%S")
    else:
        raise ValueError(f"'{datetime_str}' is not a valid datetime string.")
    return np.datetime64(pydt).astype("datetime64[ms]")


def intvlstr_to_dt64(intvl_str: str) -> np.timedelta64:
    """Convert a time interval string to a Numpy timedelta64.

    Format of string is "##[DHM]". Where ## is an integer and character
    D, H or M indicates Days, Hours or Minutes.
    """
    match = re.match(r"^(\d+)([DHMdhm])$", intvl_str)
    if match:
        value = int(match.group(1))
        unit = match.group(2).lower()
        if unit == "d":
            unit = "D"
        return np.timedelta64(value, unit).astype("timedelta64[ms]")
    else:
        raise ValueError(f"'{intvl_str}' is not a valid time interval string.")


def nwk_name(nwk_str: str) -> str:
    """Validation function for --network arg."""
    match = re.match(r"^[A-Za-z0-9_]+$", nwk_str)
    # match = re.match(r"^\w+$", nwk_str)
    if not match or len(match.group(0)) != bpr.NWK_NAME_LEN:
        # if not match:
        raise ValueError(f"Network name must be maximum {bpr.NWK_NAME_LEN} charcters.")
    return nwk_str.upper()


def stn_name(stn_str: str) -> str:
    """Validation function for --station arg."""
    match = re.match(r"^\w+$", stn_str)
    if not match or len(match.group(0)) != bpr.STN_NAME_LEN:
        raise ValueError(f"Station name must be maximum {bpr.STN_NAME_LEN} charcters.")
    return stn_str.upper()


def plot_flags(flags: str) -> dict[str:str]:
    """Validation functionfor --plot arg."""
    return {"format": flags[:1], "output": flags[1:]}
