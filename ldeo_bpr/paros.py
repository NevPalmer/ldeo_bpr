"""Module for paroscientific APG pressure gauges."""

import sys
from configparser import ConfigParser
from dataclasses import dataclass
from pathlib import Path


@dataclass()
class Paros:
    """Class for Paroscientific APG pressure gauge."""

    ser_no: str = None
    u: float = None
    y: tuple[float] = None
    c: tuple[float] = None
    d: tuple[float] = None
    t: tuple[float] = None

    @classmethod
    def from_file(cls, filename: str, paros_sn: str):
        """Create a Paros object from an .ini configuration file."""
        paros = cls.__new__(cls)
        paros_coefs = ("u", "y", "c", "d", "t")
        paros.ser_no = paros_sn
        paros_cfgs = ConfigParser()
        paros_cfgs.read(Path(filename))
        try:
            for coef in paros_coefs:
                coef_values = paros_cfgs[paros.ser_no][coef].split(",")
                setattr(paros, coef, tuple(float(x) for x in coef_values[::-1]))
        except KeyError as key:
            sys.exit(
                f"The file '{filename}' does not contain an entry for "
                f"APG sensor with serial number {key}."
            )
        paros.y = paros.y + (0.0,)  # Constant term for 'Y' is zero.
        return paros
