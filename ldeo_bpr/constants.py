"""Module conatining package global constants."""

# Dictionary of flags for turning on/off trouble shooting outputs.
# Assign True to enable and  False or '' to disable.
TROUBLE_SHOOT: dict[str, bool] = {
    # Save raw data as a text file of binary 1s & 0s. (Very slow)
    "binary_out": False,
    # Save raw data as a text file of hexadecimal values. (Very slow)
    "hex_out": False,
    # Save raw records to file as integers without tick rollover removed.
    "raw_rlovr": False,
    # Save raw records to file as integers with tick rollover removed.
    "raw_no_rlovr": False,
    # Save time sync records to file as integers with tick rollover remved.
    "raw_sync": False,
    # Write a summary file showing syncronised tic counts and exit.
    "tic_sync": False,
}

# miniSEED constants
NWK_NAME_LEN = 2
STN_NAME_LEN = 5
P_CHNL_CODE = "HDH"
T_CHNL_CODE = "BKO"

# Pressure conversion factor from PSIA to Pascal.
PSIA_TO_PASCAL = 6.894757293168e3

# Default values and choices for reading params from command line.
DFLT_APG_INI: str = "ParosAPG.ini"
DFLT_LOGGER_INI: str = "APGlogger.ini"
DFLT_CLK_START: str = "2000-01-01T00:00:00"  # 'YYYY-MM-DDThh:mm:ss'
DFLT_TMPTR_SMTH_FCTR: int = 1
DFLT_DECIMATE_INTVL: int = 0
SER_PORT: str = "COM1"
SER_BAUD: int = 9600
SER_PARITY: str = "N"
SER_STOP: int = 1
SER_BYTESIZE: int = 8
