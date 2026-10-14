import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters

from ldeo_bpr import constants as const

register_matplotlib_converters()

# filename = "GNS24-PM.csv"


def main():  # noqa: D103
    helpdesc = (
        "Generates a plot from a CSV file containing 'Time', 'Pressure', "
        "'Temperature' fields. The accepted input format is the same as the CSV "
        "file output from 'apg_read.py'."
    )
    parser = argparse.ArgumentParser(
        description=helpdesc,
    )

    parser.add_argument(
        "infile",
        help=(
            "Full path and filename of CSV input file containing columns named "
            "'Time', 'Pressure', 'Temperature'."
        ),
        type=Path,
    )

    parser.add_argument(
        "--skiprows",
        help=(
            "Number of rows to skip at the beginning of the CSV file. "
            "Default is None."
        ),
        type=int,
        default=None,
    )

    parser.add_argument(
        "--numrows",
        help=(
            "Number of rows to include from the CSV file. "
            "Default is to read to end of the file."
        ),
        type=int,
        default=None,
    )

    args = parser.parse_args()
    filename = args.infile
    skip = args.skiprows
    numrows = args.numrows

    df = pd.read_csv(
        filename,
        header=None,
        on_bad_lines="warn",
        skiprows=skip,  # one day = 86,400
        nrows=numrows,
        names=("Time", "Pressure", "Temperature"),
        parse_dates=["Time"],
        dtype={
            "Pressure": "float",
            "Temperature": "float",
        },
    )
    time = df["Time"]
    pressure_kpa = df["Pressure"] / 1000
    temperature_degc = df["Temperature"]

    # pres_median = pressure_kpa.median()
    # # pres_median = 10216.258
    # pres_spread = 25
    # pres_min = pres_median - (pres_spread / 2)
    # pres_max = pres_median + (pres_spread / 2)

    # tptr_median = temperature_degc.median()
    # # tptr_median = 5.7
    # tptr_spread = 2
    # tptr_min = tptr_median - (tptr_spread / 2)
    # tptr_max = tptr_median + (tptr_spread / 2)

    ## Set min-max values for plot Y-axes
    if const.PRESSURE_SPREAD == 0:
        pres_min = np.min(pressure_kpa)
        pres_max = np.max(pressure_kpa)
        p_range = pres_max - pres_min

        if p_range == 0:
            intvl = 10
        else:
            int_log = int(np.log10(p_range))
            if 10**int_log / p_range >= 0.5:
                intvl = 10**int_log / 10
            elif 10**int_log / p_range >= 0.2:
                intvl = 10**int_log / 5
            elif 10**int_log / p_range >= 0.1:
                intvl = 10**int_log / 2
            else:
                intvl = 10**int_log
        pres_min = pres_min - pres_min % intvl - intvl
        pres_max = pres_max - pres_max % intvl + 2 * intvl
    else:
        if const.PRESSURE_MID == 0:
            p_mid = np.median(pressure_kpa)
        else:
            p_mid = const.PRESSURE_MID
        pres_min = p_mid - (const.PRESSURE_SPREAD / 2)
        pres_max = p_mid + (const.PRESSURE_SPREAD / 2)


    if const.TEMPERATURE_SPREAD == 0:
        tptr_min = np.min(temperature_degc)
        tptr_max = np.max(temperature_degc)
        t_range = tptr_max - tptr_min
        if t_range == 0:
            intvl = 0.1
        else:
            int_log = int(np.log10(t_range))
            if 10**int_log / t_range >= 0.5:
                intvl = 10**int_log / 10
            elif 10**int_log / t_range >= 0.2:
                intvl = 10**int_log / 5
            elif 10**int_log / t_range >= 0.1:
                intvl = 10**int_log / 2
            else:
                intvl = 10**int_log
        tptr_min = tptr_min - tptr_min % intvl - intvl
        tptr_max = tptr_max - tptr_max % intvl + 2 * intvl
    else:
        if const.TEMPERATURE_MID == 0:
            t_mid = np.median(temperature_degc)
        else:
            t_mid = const.TEMPERATURE_MID
        tptr_min = t_mid - (const.TEMPERATURE_SPREAD / 2)
        tptr_max = t_mid + (const.TEMPERATURE_SPREAD / 2)

    # Plot pressure and temperature vs time
    ##################################
    fig, ax1 = plt.subplots(figsize=(15, 9))
    color = "red"
    plt.ylim(pres_max, pres_min)
    ax1.plot(time, pressure_kpa, "r.", markersize=1, label="Pressure")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Pressure (kPa)", color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.grid(axis="x")
    ax1.legend(loc="upper left")

    ax2 = ax1.twinx()
    color = "blue"
    plt.ylim(tptr_min, tptr_max)
    ax2.plot(time, temperature_degc, "b.", markersize=1, label="Temperature")
    ax2.set_ylabel("Temperature (°C)", color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    ax2.legend(loc="upper right")

    plt.gcf().autofmt_xdate()
    fig.suptitle(filename)
    fig.tight_layout()

    basename = filename.stem
    pathname = filename.parent
    fig.savefig(f"{pathname}/{basename}.png", dpi=200)

    plt.show()


if __name__ == "__main__":
    main()
