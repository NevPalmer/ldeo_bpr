#!/usr/bin/env python3

# Version 20200818
# by Neville Palmer, GNS Science
# 2020/08/18 Start development based on apg_read.py

import os
import sys
from subprocess import call
import argparse
import configparser
import serial
import re
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from numexpr import evaluate


def main():
    # Default values and choices for reading params from command line.
    paros_ini = './ParosAPG.ini'
    logger_ini = './APGlogger.ini'
    logger_versions = ['CSAC2013','Seascan2018']
    ser_port = 'COM11'
    ser_baud = 9600
    ser_parity = 'N'
    ser_stop = 1
    ser_bytesize = 8
    ser_timeout = 30

    # Other constants.
    # Pressure conversion factor from PSIA to Pascal.
    press_conv_fctr = 6.894757293168e3

    # Read in parameters from command line
    helpdesc = (
        f'Reads a raw APG data file and outputs decimated pressure data.'
        f'Two .ini files are required which must contain configuration values '
        f'for the specific Paroscientific pressure transducer and the correct '
        f'version of APG logger board used.')
    parser = argparse.ArgumentParser(description = helpdesc)
    parser.add_argument("-p","--parosini", 
        help=f'Full path and filename for Paros APG configuration settings. '
            f'Default: "{paros_ini}"',
        default=paros_ini)
    parser.add_argument("-s","--snparos", 
        help=f'Serial number of the Paroscientific APG used. This must correspond '
            f'to the serial number of an entry in the parosini file.',
        required=True)
    parser.add_argument("-l","--loggerini", 
        help=f'Full path and filename for the APG logger board configuration '
            f'settings. Default: "{logger_ini}"',
        default=logger_ini)
    parser.add_argument("-v","--version", 
        help=f'Specify the version/firmware of the APG logger board used.',
        choices=logger_versions,
        required=True)
    parser.add_argument("--serport", 
        help=f'Serial port name. Default: "{ser_port}"',
        default=ser_port)
    parser.add_argument("--serbaud", 
        help=f'Serial baud rate. Default: "{ser_baud}"',
        default=ser_baud)
    parser.add_argument("--serparity", 
        help=f'Serial parity. Default: "{ser_parity}"',
        default=ser_parity)
    parser.add_argument("--serstop", 
        help=f'Serial stop bit name. Default: "{ser_stop}"',
        default=ser_stop)
    parser.add_argument("--serbytesize", 
        help=f'Serial byte size. Default: "{ser_bytesize}"',
        default=ser_bytesize)
    args = parser.parse_args()
    
    # Translate argparse parameters (except for window times).
    paros_ini = os.path.normpath(args.parosini)
    paros_sn = args.snparos.strip()
    logger_ini = os.path.normpath(args.loggerini)
    logger_version = args.version.strip()
    ser_port = args.serport.strip()
    ser_parity = args.serparity.strip()

    # Read Paros transducer coefficients into a dict of lists from ini file.
    paros_coefs = ('U', 'Y', 'C', 'D', 'T')
    paros_cfg = configparser.ConfigParser()
    paros_cfg.read(paros_ini)
    paros = {}
    for coef in paros_coefs:
        paros[coef] = (paros_cfg.get(paros_sn,coef).split(','))
        paros[coef] = [float(x) for x in paros[coef][::-1]]
    paros['Y'].append(0.0)
    
    # Read APG logger configuration parameters from ini file.
    logger_cfg = configparser.ConfigParser()
    logger_cfg.read(logger_ini)
    clock_freq = logger_cfg.getint(logger_version,'clock_freq')
    TP_fctr = evaluate(logger_cfg.get(logger_version,'TP_fctr')).item()
    TP_cnst = logger_cfg.getfloat(logger_version,'TP_cnst')
    PP_fctr = evaluate(logger_cfg.get(logger_version,'PP_fctr')).item()
    PP_cnst = logger_cfg.getfloat(logger_version,'PP_cnst')
    
    
    ser = serial.Serial(
        port=ser_port,
        baudrate=ser_baud,
        parity=ser_parity,
        stopbits=ser_stop,
        bytesize=ser_bytesize,
        timeout=ser_timeout
    )
    print("connected to: " + ser.portstr)
    
    while True:
        raw_ticks = ser.read(5)
        raw_tempr = ser.read(5)
        raw_press = ser.read(5)
        int_ticks = int.from_bytes(raw_ticks, byteorder='big')
        int_tempr = int.from_bytes(raw_tempr, byteorder='big')
        int_press = int.from_bytes(raw_press, byteorder='big')
        print(f'raw_ticks: {int_ticks} | {raw_ticks}   ',end='')
        print(f'raw_tempr: {int_tempr} | {raw_tempr}   ',end='')
        print(f'raw_press: {int_press} | {raw_press}')

        time_secs = int_ticks / 1000.0
        dt_time = dt.timedelta(seconds = time_secs)
        hours, remainder = divmod(time_secs, 3600)
        minutes, seconds = divmod(remainder, 60)
        print (f'Time since start: {int(hours):02}:{int(minutes):02}:{seconds:06.3f} ',end='')
        print(f'({time_secs:.3f} secs)   ',end='')
        
        # Temperature period (usec)
        TP = (int_tempr/(TP_fctr)+TP_cnst)/(clock_freq)

        # Pressure period (usec)
        PP = (int_press/(PP_fctr)+PP_cnst)/(clock_freq)

        # Create temperature
        Uv = TP-paros['U'][0]
        temperature = np.polyval(paros['Y'], Uv)
        
        # Create pressure
        Cv = np.polyval(paros['C'], Uv)
        Dv = np.polyval(paros['D'], Uv)
        T0 = np.polyval(paros['T'], Uv)
        
        facts = 1 - (T0**2)/(PP**2) 
        pressure_psi = Cv * facts *(1-Dv*facts) #pressure in PSIA
        pressure_kPa = pressure_psi * press_conv_fctr / 1000.0 #Convert pressure units
        
        print(f'Temperature: {temperature:.3f}°C   ',end='')
        print(f'Pressure: {pressure_kPa:.4f} kPa ',end='')
        print(f'/ {pressure_psi:.5f} PSI')
            
    ser.close()

    return


##########################
if __name__ == '__main__':
    main()