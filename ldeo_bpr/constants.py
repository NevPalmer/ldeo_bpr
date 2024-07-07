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
PRESS_CONV_FCTR = 6.894757293168e3
