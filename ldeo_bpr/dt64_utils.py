"""Utilities for working with and converting to/from numpy datetime64."""

import datetime as dt
import re

import numpy as np


###############################################################################
# The following are conversions to and from numpy datetime64 & timedelta64
def dt64_to_pydt(dt64):
    """Convert a numpy datetime64 to a py datetime."""
    seconds_since_epoch = dt64_to_ms(dt64) / 1000
    return dt.datetime.fromtimestamp(seconds_since_epoch, tz=dt.timezone.utc)


def pydt_to_dt64(pydt):
    """Convert a py datetime to a numpy datetime64."""
    return np.datetime64(pydt).astype("int64")


def delta64_to_ms(delta64):
    """Convert a numpy timedelta64 to milliseconds."""
    return delta64.astype("timedelta64[ms]").astype("int64")


def ms_to_delta64(msec):
    """Convert numpy array of milliseconds to timedelta64."""
    return msec.astype("timedelta64[ms]")


def dt64_to_ms(dt64):
    """Convert a numpy datetime64 to milliseconds."""
    return dt64.astype("datetime64[ms]").astype("int64")


def ms_to_dt64(msec):
    """Convert milliseconds to a numpy datetime64."""
    return np.datetime64(int(msec), "ms")


###############################################################################
# The following are conversions from strings to numpy datetime64 & timedelta64
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
