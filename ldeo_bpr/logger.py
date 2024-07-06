"""Module for LDEO BPR data loggers."""

import sys
from configparser import ConfigParser
from dataclasses import dataclass
from pathlib import Path


@dataclass()
class Logger:
    """Class for LDEO BPR data loggers.."""

    version: str = None
    head_len: int = None
    rec_len: int = None
    smpls_per_rec: int = None
    sample_epoch: int = None
    record_epoch: int = None
    clock_freq: int = None
    tp_fctr: int = None
    tp_cnst: int = None
    pp_fctr: int = None
    pp_cnst: int = None
    timing: str = None
    rec_fmt: tuple[int] = None
    fmt_field: dict[str, int] = None
    tic_bit_len: int = None

    @classmethod
    def from_file(cls, filename: str, logger_version: str):
        """Create a Logger object from an .ini configuration file."""
        logger = cls.__new__(cls)
        logger.version = logger_version
        all_logger_cfgs = ConfigParser()
        all_logger_cfgs.read(Path(filename))
        logger_cfg = dict(all_logger_cfgs[logger_version])
        logger.head_len = int(logger_cfg["head_len"])
        logger.rec_len = int(logger_cfg["rec_len"])
        logger.smpls_per_rec = int(logger_cfg["smpls_per_rec"])
        logger.sample_epoch = int(logger_cfg["epoch"])
        logger.record_epoch = logger.sample_epoch * logger.smpls_per_rec
        logger.clock_freq = int(logger_cfg["clock_freq"])
        logger.tp_fctr = _eval_exponent_str(logger_cfg["tp_fctr"])
        logger.tp_cnst = int(logger_cfg["tp_cnst"])
        logger.pp_fctr = _eval_exponent_str(logger_cfg["pp_fctr"])
        logger.pp_cnst = int(logger_cfg["pp_cnst"])
        logger.timing = logger_cfg["timing"]
        logger.rec_fmt = tuple([int(x) for x in logger_cfg["rec_fmt"].split(",")])

        fmt_field = {}
        fmt_field["tic"] = int(logger_cfg["tic_field"])
        fmt_field["tptr"] = int(logger_cfg["temperature_field"])
        fmt_field["pcore"] = int(logger_cfg["pcore_field"])
        fmt_field["pn"] = int(logger_cfg["pn_field"])
        if abs(fmt_field["pcore"] - fmt_field["pn"] + 1) != logger.smpls_per_rec:
            sys.exit(
                f"The number of samples per record "
                f"({logger.smpls_per_rec}) for logger {logger.version}, "
                f"does not match the number of Pcore and Pn fields in the "
                f"record format (Pcore through Pn inclusive = "
                f"{fmt_field['pcore']-fmt_field['pn']+1}), "
                f"as provided in file {filename}."
            )
        rec_fmt_bits = int(sum(map(abs, logger.rec_fmt)))
        if rec_fmt_bits != (logger.rec_len * 8):
            sys.exit(
                f"The total number of bits ({rec_fmt_bits}) given by the "
                f"record format {logger.rec_fmt} "
                f"does not equal the record length ({logger.rec_len} "
                f"bytes x 8), as provided in the file {filename}."
            )

        logger.fmt_field = fmt_field
        logger.tic_bit_len = logger.rec_fmt[fmt_field["tic"]]

        return logger


###############################################################################
# Miscellaneous utility functions.
def _eval_exponent_str(exp_str: str) -> int:
    """Convert a string in the form "X**Y" to an int.

    Will also accept a pure integer string ("X") as input.
    """
    try:
        return int(exp_str)
    except ValueError:
        base, exponent = exp_str.split("**", 1)
        try:
            return int(base) ** int(exponent)
        except ValueError as err:
            raise ValueError(
                f"'{exp_str}' is not a valid exponent string to be evaluated."
            ) from err
