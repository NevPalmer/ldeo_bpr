#!/usr/bin/env python3  # noqa: D100

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import register_matplotlib_converters

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

    args = parser.parse_args()
    filename = args.infile

    csv_datatypes = {
        "Pressure": "float",
        "Temperature": "float",
    }
    df = pd.read_csv(
        filename,
        header=None,
        # skiprows=86_400,  # one day = 86,400
        # nrows=86_000,
        names=("Time", "Pressure", "Temperature"),
        parse_dates=["Time"],
        dtype=csv_datatypes,
    )
    time = df["Time"]
    pressure_kpa = df["Pressure"] / 1000
    temperature_degc = df["Temperature"]

    pres_median = pressure_kpa.median()
    # pres_median = 10216.258
    pres_spread = 25
    pres_min = pres_median - (pres_spread / 2)
    pres_max = pres_median + (pres_spread / 2)

    tptr_median = temperature_degc.median()
    # tptr_median = 5.7
    tptr_spread = 2
    tptr_min = tptr_median - (tptr_spread / 2)
    tptr_max = tptr_median + (tptr_spread / 2)

    # Plot voltage and current vs time
    ##################################
    fig, ax1 = plt.subplots()
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
    plt.show()


if __name__ == "__main__":
    main()
