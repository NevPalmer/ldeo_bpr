"""Module to parse arguments from command line."""

import re
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Iterable

import numpy as np

from . import constants as const
from . import dt64_utils
from .version import __version__


def parse_args_apg_read(args: Iterable[str] = None) -> Namespace:
    """Parse CLI arguments for apg_read.py."""
    helpdesc = (
        "Reads a raw APG data file and outputs decimated pressure data."
        "Two .ini files are required which must contain configuration values "
        "for the specific Paroscientific pressure transducer and the correct "
        "version of APG logger board used."
    )
    parser = ArgumentParser(
        parents=[
            base_parser(),
            infile_parser(),
            ini_cfg_parser(),
            processing_param_parser(),
            outfile_parser(),
            plot_param_parser(),
        ],
        description=helpdesc,
    )
    args = parser.parse_args(args)

    args.plot = {"format": args.plot[:1], "output": args.plot[1:]}

    if args.mseedpath and not (args.network and args.station):
        parser.error(
            "Arguments -n/--station and -w/--network are both required because "
            "-m/--mseedpath was specified."
        )

    return args


def parse_args_apg_read_serial(args: Iterable[str] = None) -> Namespace:
    """Parse CLI arguments for apg_read_serial.py."""
    helpdesc = (
        "Reads a raw APG data file and outputs decimated pressure data."
        "Two .ini files are required which must contain configuration values "
        "for the specific Paroscientific pressure transducer and the correct "
        "version of APG logger board used."
    )
    parser = ArgumentParser(
        parents=[
            base_parser(),
            ini_cfg_parser(),
            serial_param_parser(),
        ],
        description=helpdesc,
    )
    args = parser.parse_args(args)

    return args


def base_parser() -> None:
    """Return standard base parser common to all."""
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def infile_parser() -> None:
    """Return parser to specify input file."""
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "-i",
        "--infile",
        help="Full path and filename of raw APG input file.",
        type=Path,
        required=True,
    )

    return parser


def ini_cfg_parser() -> None:
    """Return parser for .ini configuration files."""
    parser = ArgumentParser(add_help=False)

    package_path: Path = Path(__file__).parents[1]
    apg_ini = Path(".") / const.DFLT_APG_INI
    if not apg_ini.is_file():
        apg_ini: Path = package_path / const.DFLT_APG_INI
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
        "This must correspond to an entry in the file specified by --apgini.",
        type=str,
        required=True,
    )

    logger_ini = Path(".") / const.DFLT_LOGGER_INI
    if not logger_ini.is_file():
        logger_ini: Path = package_path / const.DFLT_LOGGER_INI
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
        "--loggerversion",
        help="Specify the version/firmware of the APG logger board used."
        "This must correspond to an entry in the file specified by "
        "--loggerversion.",
        type=str,
        required=True,
    )

    return parser


def processing_param_parser() -> None:
    """Return parser to retrieve processing parameters."""
    parser = ArgumentParser(add_help=False)

    parser.add_argument(
        "-d",
        "--decimate",
        help=f"Required sample interval in seconds for "
        f"pressure decimation. Zero for no decimation. "
        f"Value must equal a single digit integer of seconds "
        f"or minutes or a multiple of 5 or 10."
        f'Default: "{const.DFLT_DECIMATE_INTVL}"',
        type=int,
        default=const.DFLT_DECIMATE_INTVL,
    )

    parser.add_argument(
        "-t",
        "--tempsmth",
        help=f"Temperature smoothing factor (must be an odd "
        f"integer). 5001 gives sensible smoothing. 50001 "
        f"gives better smoothing for Seascan logger but is "
        f'slow. Default: "{const.DFLT_TMPTR_SMTH_FCTR}"',
        type=int,
        default=const.DFLT_TMPTR_SMTH_FCTR,
    )

    parser.add_argument(
        "-c",
        "--clkstart",
        help=f"Precise date and time when the logger clock "
        f'was started. Format: "YYYY-MM-DDThh:mm:ss" or "YYYY:DDD:hh:mm:ss" '
        f"(optionally include decimal seconds)"
        f'Default: "{const.DFLT_CLK_START}"',
        type=dt64_utils.dtstr_to_dt64,
        default=const.DFLT_CLK_START,
    )

    parser.add_argument(
        "-b",
        "--beginwndw",
        help="Date and time to begin data extraction. "
        "Assumes beginning of file if omitted. "
        'Format: "YYYY-MM-DDThh:mm:ss" or "YYYY:DDD:hh:mm:ss" '
        "(optionally include decimal seconds)",
        type=dt64_utils.dtstr_to_dt64,
    )

    parser.add_argument(
        "-e",
        "--endwndw",
        help="Date and time to end data extraction. Assumes "
        "end of file if omitted. "
        'Format: "YYYY-MM-DDThh:mm:ss" or "YYYY:DDD:hh:mm:ss" '
        "(optionally include decimal seconds)",
        type=dt64_utils.dtstr_to_dt64,
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
        type=dt64_utils.intvlstr_to_delta64,
        default=np.timedelta64(0, "ms"),
    )

    parser.add_argument(
        "-g",
        "--gpssynctime",
        help="Precise date and time from GPS clock for "
        "syncronising end time. No clock drift adjustment is "
        'made if omitted. Format: "YYYY-MM-DDThh:mm:ss" or "YYYY:DDD:hh:mm:ss" '
        "(optionally include decimal seconds)",
        type=dt64_utils.dtstr_to_dt64,
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
        "--noisefilt",
        help="Apply a binned median to pressure & temperature data and remove"
        "any data points that are greater than a predefined range from the"
        "median values.",
        action="store_true",
        default=False,
    )

    return parser


def outfile_parser() -> None:
    """Return parser for output file parameters."""
    parser = ArgumentParser(add_help=False)

    parser.add_argument(
        "-w",
        "--network",
        help=(
            f"Network name to be used in MiniSEED file header. Max "
            f"{const.NWK_NAME_LEN} characters."
        ),
        type=nwk_name,
        default=None,
    )

    parser.add_argument(
        "-n",
        "--station",
        help=(
            f"Station name to be used in MiniSEED file header. Max "
            f"{const.STN_NAME_LEN} characters."
        ),
        type=stn_name,
        default=None,
    )

    out_filename: Path | None = None
    parser.add_argument(
        "-o",
        "--outfile",
        help="Full path and filename for output file. No file "
        "will be generated if not specified.",
        type=Path,
        default=out_filename,
    )

    mseed_path: Path | None = None
    parser.add_argument(
        "-m",
        "--mseedpath",
        help="Full path of location to save MiniSEED file(s). No file(s) will "
        "be generated if not specified. If specifed, then both"
        "-n/--station and -w/--network are required.",
        type=Path,
        default=mseed_path,
    )

    return parser


def plot_param_parser() -> None:
    """Return parser to retrieve plotting parameters."""
    parser = ArgumentParser(add_help=False)

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

    return parser


def serial_param_parser() -> None:
    """Return parser to retrieve plotting parameters."""
    parser = ArgumentParser(add_help=False)

    parser.add_argument(
        "--serport",
        help=f'Serial port name. Default: "{const.SER_PORT}"',
        type=str,
        default=const.SER_PORT,
    )
    parser.add_argument(
        "--serbaud",
        help=f"Serial baud rate. Default: {const.SER_BAUD}",
        type=int,
        default=const.SER_BAUD,
    )
    parser.add_argument(
        "--serparity",
        help=f'Serial parity. Default: "{const.SER_PARITY}"',
        type=str,
        default=const.SER_PARITY,
    )
    parser.add_argument(
        "--serstop",
        help=f"Serial stop bit. Default: {const.SER_STOP}",
        type=int,
        default=const.SER_STOP,
    )
    parser.add_argument(
        "--serbytesize",
        help=f"Serial byte size. Default: {const.SER_BYTESIZE}",
        type=int,
        default=const.SER_BYTESIZE,
    )

    return parser


def nwk_name(nwk_str: str) -> str:
    """Validation function for --network arg."""
    match = re.match(r"^[A-Za-z0-9_]+$", nwk_str)
    # match = re.match(r"^\w+$", nwk_str)
    if not match or len(match.group(0)) != const.NWK_NAME_LEN:
        # if not match:
        raise ValueError(
            f"Network name must be maximum {const.NWK_NAME_LEN} charcters."
        )
    return nwk_str.upper()


def stn_name(stn_str: str) -> str:
    """Validation function for --station arg."""
    match = re.match(r"^\w+$", stn_str)
    if not match or len(match.group(0)) != const.STN_NAME_LEN:
        raise ValueError(
            f"Station name must be maximum {const.STN_NAME_LEN} charcters."
        )
    return stn_str.upper()


### Choices are checked after type conversion performed so this fails.
# def plot_flags(flags: str) -> dict[str, str]:
#     """Validation functionfor --plot arg."""
#     return {"format": flags[:1], "output": flags[1:]}
