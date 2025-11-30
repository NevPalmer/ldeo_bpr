"""Read raw realtime serial data from LDEO type APG data loggers.

A script for reading and converting raw realtime serial data from LDEO type
APG data loggers.
"""

import sys

import numpy as np
import serial

from ldeo_bpr import arg_parser
from ldeo_bpr import constants as const
from ldeo_bpr.logger import Logger
from ldeo_bpr.paros import Paros


def main():
    """The first function run when this script is run directly."""
    # Retrieve CLI parameters.
    args = arg_parser.parse_args_apg_read_serial()

    # Default serail timeout for serial connection.
    serial_timeout = 0.1

    # Read in parameters from command line
    args = arg_parser.parse_args_apg_read_serial()

    # Read Paros transducer coefficients from .ini file.
    paros = Paros.from_file(filename=args.apgini, paros_sn=args.snapg)

    # Read APG logger configuration parameters from .ini file.
    logger = Logger.from_file(
        filename=args.loggerini, logger_version=args.loggerversion
    )

    try:
        with serial.Serial(
            port=args.serport,
            baudrate=args.serbaud,
            parity=args.serparity,
            stopbits=args.serstop,
            bytesize=args.serbytesize,
            timeout=serial_timeout,
        ) as ser:
            print(f"connected to: {ser.portstr} at {ser.baudrate}bps.")

            while True:
                raw_record = []
                raw_record_len = 15
                next_byte = ser.read(1)
                while next_byte != b"":  # next_byte will be '' after serial.timeout
                    raw_record.append(next_byte)
                    if len(raw_record) > raw_record_len:
                        raw_record = raw_record[-raw_record_len:]
                    next_byte = ser.read(1)

                if raw_record:
                    raw_ticks = b"".join(raw_record[0:5])
                    raw_tempr = b"".join(raw_record[5:10])
                    raw_press = b"".join(raw_record[10:15])
                    int_ticks = int.from_bytes(raw_ticks, byteorder="big")
                    int_tempr = int.from_bytes(raw_tempr, byteorder="big")
                    int_press = int.from_bytes(raw_press, byteorder="big")
                    print()
                    print(f"raw_ticks: {int_ticks:11} | 0x{int_ticks:010X}")
                    print(f"raw_tempr: {int_tempr:11} | 0x{int_tempr:010X}")
                    print(f"raw_press: {int_press:11} | 0x{int_press:010X}")

                    time_secs = int_ticks / 1000.0
                    hours, remainder = divmod(time_secs, 3600)
                    minutes, seconds = divmod(remainder, 60)
                    print(
                        f"Time since start: "
                        f"{int(hours):02}:{int(minutes):02}:{seconds:06.3f} ",
                        end="",
                    )
                    print(f"({time_secs:.3f} secs)   ")

                    # Temperature period (usec)
                    tmptr_period_usec = (
                        int_tempr / (logger.tp_fctr) + logger.tp_cnst
                    ) / (logger.clock_freq)

                    # Pressure period (usec)
                    presr_period_usec = (
                        int_press / (logger.pp_fctr) + logger.pp_cnst
                    ) / (logger.clock_freq)

                    # Create temperature
                    coef_uv = tmptr_period_usec - paros.u
                    temperature = np.polyval(paros.y, coef_uv)

                    # Create pressure
                    coef_cv = np.polyval(paros.c, coef_uv)
                    coef_dv = np.polyval(paros.d, coef_uv)
                    coef_t0 = np.polyval(paros.t, coef_uv)

                    facts = 1 - (coef_t0**2) / (presr_period_usec**2)
                    pressure_psi = (
                        coef_cv * facts * (1 - coef_dv * facts)
                    )  # pressure in PSIA
                    # Convert pressure units
                    pressure_kpa = pressure_psi * const.PSIA_TO_PASCAL / 1000.0
                    # Convert temperature units
                    temperature_f = temperature * 9 / 5 + 32

                    print(f"Temperature: {temperature:.3f}°C ", end="")
                    print(f"/ {temperature_f:.3f}°F   ")
                    print(f"Pressure: {pressure_kpa:.4f} kPa ", end="")
                    print(f"/ {pressure_psi:.5f} PSI")

    except serial.SerialException as error:
        print(error)
        sys.exit()
    except ValueError as error:
        print(error)
        sys.exit()


##########################
if __name__ == "__main__":
    main()
