#!/usr/bin/env python3

# Version 20181207
# by Neville Palmer, GNS Science
# 2018/11/17 Start development
# 2018/11/23 SeaScan board format added.
# 2018/11/24 Smoothing and decimation added.
# 2018/11/25 Added ini files for reading in transducer & logger parameters.
# 2018/11/26 Resolved issue with iir decimation (ensure q<13 and n=5).
# 2018/11/28 Add iteration routine for decimation.
# 2018/12/02 Add command line parameters.
# 2018/12/07 Add output tic count summary file.
# 2018/12/14 Add ability to calc clock drift for CSAC logger.
# 2020/08/26 Minor edits for plotting, output file and debug outputs.
# 2020/11/30 Added parameter to plot/save time as either seconds or date-time.
#            Added linear interpolation for temperature upsampling and included
#            temperature in same decimation loop as pressure.

import os
import sys
import argparse
import configparser
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
    logger_versions = ['CSAC2013', 'Seascan2018']
    clk_start = '2000-01-01_00:00:00'  # 'YYYY-MM-DD_hh:mm:ss'
    out_filename = ''
    tmptr_smth_fctr = 1
    decmt_intvl = 0

    # Other constants.
    # Pressure conversion factor from PSIA to Pascal.
    press_conv_fctr = 6.894757293168e3

    # Read in parameters from command line
    helpdesc = (
        'Reads a raw APG data file and outputs decimated pressure data.'
        'Two .ini files are required which must contain configuration values '
        'for the specific Paroscientific pressure transducer and the correct '
        'version of APG logger board used.')
    parser = argparse.ArgumentParser(description=helpdesc)
    parser.add_argument("-i", "--infile",
                        help='Full path and filename of raw APG input file.',
                        required=True)
    parser.add_argument("-p", "--parosini",
                        help='Full path and filename for Paros APG '
                        f'configuration settings. Default: "{paros_ini}"',
                        default=paros_ini)
    parser.add_argument("-s", "--snparos",
                        help='Serial number of the Paroscientific APG used. '
                        'This must correspond to the serial number of an '
                        'entry in the parosini file.',
                        required=True)
    parser.add_argument("-l", "--loggerini",
                        help=f'Full path and filename for the APG logger '
                        f'board configuration settings. '
                        f'Default: "{logger_ini}"',
                        default=logger_ini)
    parser.add_argument("-v", "--version",
                        help='Specify the version/firmware of the APG logger '
                        'board used.',
                        choices=logger_versions,
                        required=True)
    parser.add_argument("-d", "--decimate",
                        help=f'Sample interval in seconds for pressure '
                        f'decimation. Zero for no decimation. '
                        f'Value must equal a single digit integer of seconds '
                        f'or minutes or a multiple of 5 or 10.'
                        f'Default: "{decmt_intvl}"',
                        type=int,
                        default=decmt_intvl)
    parser.add_argument("-t", "--tempsmth",
                        help=f'Temperature smoothing factor (must be an odd '
                        f'integer). 10001 gives strong smoothing. '
                        f'Default: "{tmptr_smth_fctr}"',
                        type=int,
                        default=tmptr_smth_fctr)
    parser.add_argument("-c", "--clkstart",
                        help=f'Precise date and time when the logger clock '
                        f'was started. Default: "{clk_start}"',
                        default=clk_start)
    parser.add_argument("-b", "--beginwndw",
                        help='Date and time to begin data extraction. '
                        'Assumes beginning of file if omitted. '
                        'Format: "YYYY-MM-DD_hh:mm:ss.s"')
    parser.add_argument("-e", "--endwndw",
                        help='Date and time to end data extraction. Assumes '
                        'end of file if omitted. '
                        'Format: "YYYY-MM-DD_hh:mm:ss.s"')
    parser.add_argument("-g", "--gpssynctime",
                        help='Precise date and time from GPS clock for '
                        'syncronising end time. No clock drift adjustment is '
                        'made if omitted. Format: "YYYY-DDD_hh:mm:ss"')
    parser.add_argument("-y", "--synctickcount",
                        help='The hexidecimal tick count that corresponds to '
                        'GPSSYNCTIME. If GPSSYNCTIME is specified and '
                        'SYNCTICKCOUNT is omitted, then it is assumed that an '
                        'artificial frequency was inserted precisely '
                        'at GPSSYNCTIME. This parameter is ignored if '
                        'GPSSYNCTIME is omitted. Format: "0xHHHHHHHHHH"',
                        type=lambda x: int(x, 0))
    parser.add_argument("-o", "--outfile",
                        help='Full path and filename for output file. No file '
                        'will be generated if not specified.',
                        default=out_filename)
    parser.add_argument("-x", "--xaxis",
                        help='Specify the scale format for the X-axis of the '
                        'plot to be displayed as either (s)econds or (d)ate-time.',
                        choices=['s','d'],
                        default='s')
    args = parser.parse_args()

    # Translate argparse parameters (except for window times).
    apg_filename = os.path.normpath(args.infile)
    paros_ini = os.path.normpath(args.parosini)
    paros_sn = args.snparos.strip()
    logger_ini = os.path.normpath(args.loggerini)
    logger_version = args.version.strip()
    decmt_intvl = args.decimate
    tmptr_smth_fctr = args.tempsmth
    clk_start = re.sub('[-: _/]', '_', args.clkstart)
    clk_start_dt = dt.datetime.strptime(clk_start, '%Y_%m_%d_%H_%M_%S')
    if args.gpssynctime:
        gpssynctime = re.sub('[-: _/]', '_', args.gpssynctime)
        gpssync_dt = dt.datetime.strptime(gpssynctime, '%Y_%j_%H_%M_%S')
    sync_tick_count = args.synctickcount
    if args.outfile:
        out_filename = os.path.normpath(args.outfile)
    xaxis_format = args.xaxis.strip()

    # Read Paros transducer coefficients into a dict of lists from ini file.
    paros_coefs = ('U', 'Y', 'C', 'D', 'T')
    paros_cfg = configparser.ConfigParser()
    paros_cfg.read(paros_ini)
    paros = {}
    for coef in paros_coefs:
        paros[coef] = (paros_cfg.get(paros_sn, coef).split(','))
        paros[coef] = [float(x) for x in paros[coef][::-1]]
    paros['Y'].append(0.0)

    # Read APG logger configuration parameters from ini file.
    logger_cfg = configparser.ConfigParser()
    logger_cfg.read(logger_ini)
    clock_freq = logger_cfg.getint(logger_version, 'clock_freq')
    sample_epoch = logger_cfg.getfloat(logger_version, 'epoch')
    rec_len = logger_cfg.getint(logger_version, 'rec_len')
    head_len = logger_cfg.getint(logger_version, 'head_len')
    smpls_per_rec = logger_cfg.getint(logger_version, 'smpls_per_rec')
    TP_fctr = evaluate(logger_cfg.get(logger_version, 'TP_fctr')).item()
    TP_cnst = logger_cfg.getfloat(logger_version, 'TP_cnst')
    PP_fctr = evaluate(logger_cfg.get(logger_version, 'PP_fctr')).item()
    PP_cnst = logger_cfg.getfloat(logger_version, 'PP_cnst')

    rec_fmt = (logger_cfg.get(logger_version, 'rec_fmt').split(','))
    rec_fmt = tuple([float(x) for x in rec_fmt])
    fmt_fields = {}
    fmt_fields['tic'] = logger_cfg.getint(logger_version, 'tic_field')
    fmt_fields['tptr'] = logger_cfg.getint(logger_version, 'temperature_field')
    fmt_fields['pcore'] = logger_cfg.getint(logger_version, 'pcore_field')
    fmt_fields['pn'] = logger_cfg.getint(logger_version, 'pn_field')
    if (abs(fmt_fields['pcore']-fmt_fields['pn']+1) != smpls_per_rec):
        sys.exit(f"The number of samples per record ({smpls_per_rec}) "
                 f"for logger {logger_version}, does "
                 f"not match the number of Pcore and Pn fields in the "
                 f"record format (Pcore through Pn inclusive = "
                 f"{fmt_fields['pcore']-fmt_fields['pn']+1}), "
                 f"as provided in file {logger_ini}.")
    rec_fmt_bits = int(sum(map(abs, rec_fmt)))
    if (rec_fmt_bits != (rec_len*8)):
        sys.exit(f"The total number of bits ({rec_fmt_bits}) given by the "
                 f"record format {rec_fmt} "
                 f"does not equal the record length ({rec_len} bytes x 8), "
                 f"as provided in the file {logger_ini}.")

    # Calculate duration and end time of raw file.
    try:
        fsize = os.path.getsize(apg_filename)  # file size in bytes
    except(FileNotFoundError):
        sys.exit(f'Raw APG file "{apg_filename}" does not exist.')
    print(f'Filesize = {fsize} bytes')
    nrecs = int((fsize - head_len) / rec_len)-1
    print(f'Number of records = {nrecs:d}')
    file_duration_days = nrecs*smpls_per_rec*sample_epoch/3600/24
    print(f'File duration = {file_duration_days} days')
    clk_end_dt = clk_start_dt + dt.timedelta(days=file_duration_days)
    print(f'Time of last sample = {clk_end_dt}')

    # Clock drift
    if (args.gpssynctime is not None):
        drift = clockdrift(clk_start_dt, gpssync_dt, sync_tick_count,
                           rec_fmt[fmt_fields['tic']])
        print(f'Clock drift in milliseconds (logged time minus actual '
              f'GPS time) = {drift}')

    # Calculate window for extracting data.
    if args.beginwndw is None:
        wndw_start_dt = clk_start_dt
    else:
        wndw_start = re.sub('[-: _/]', '_', args.beginwndw)
        wndw_start_dt = dt.datetime.strptime(wndw_start,
                                             '%Y_%m_%d_%H_%M_%S.%f')
    if args.endwndw is None:
        wndw_end_dt = clk_end_dt
    else:
        wndw_end = re.sub('[-: _/]', '_', args.endwndw)
        wndw_end_dt = dt.datetime.strptime(wndw_end, '%Y_%m_%d_%H_%M_%S.%f')
    wndw_start_secs = ((wndw_start_dt - clk_start_dt).total_seconds())
    wndw_start_days = wndw_start_secs / (24*3600)
    print(f'wndw_start_dt = {wndw_start_dt}')
    print(f'wndw_start_days = {wndw_start_days}')
    wndw_len_secs = ((wndw_end_dt - wndw_start_dt).total_seconds())
    wndw_len_days = wndw_len_secs / (24*3600)
    print(f'wndw_len_days = {wndw_len_days}')
    rec_start = int(wndw_start_secs / sample_epoch / smpls_per_rec)
    print(f'rec_start = {rec_start}')
    nrecs_want = int(wndw_len_secs / sample_epoch / smpls_per_rec)+1

    # Extract records from APG file for the time window specified.
    start_byte = head_len + rec_start * rec_len
    records = extractrecords(apg_filename, start_byte, nrecs_want,
                             rec_len, rec_fmt)

    # Assign extracted records to field names.
    last_field = len(rec_fmt) - 1
    pcore = last_field - fmt_fields['pcore']
    pn = last_field - fmt_fields['pn']
    tptr = last_field - fmt_fields['tptr']
    tic = last_field - fmt_fields['tic']
    if pn > pcore:
        press_raw = records[:, pcore:pn+1]
    else:
        press_raw = records[:, pcore:pn-1:-1]
    press_raw = np.cumsum(press_raw, axis=1)
    press_raw = press_raw.reshape((nrecs_want*smpls_per_rec))
    temp_raw = (records[:, tptr])
    ticks = (records[:, tic])

    # Write a summary file showing out of sync tic counts and exit.
    # Uncomment below for troubleshooting.
    '''
    #nom_time = np.arange(64,64+72*np.alen(ticks),72,dtype=int) #Seascan
    nom_time = np.arange(0,80*np.alen(ticks),80,dtype=int) #CSAC
    time_diff = ticks - nom_time
    ticks = np.column_stack((ticks,nom_time,time_diff))
    summary_ticks = [ticks[0]]
    for i in range(1,len(time_diff)):
        if (time_diff[i] != time_diff[i-1]):
            print(f'{time_diff[i]} : {time_diff[i-1]} : '
                  f'{time_diff[i] != time_diff[i-1]}')
            summary_ticks.append(ticks[i])
    np.savetxt('summary_ticks.txt', summary_ticks, fmt = '%d')
    np.savetxt('ticks.txt', ticks, fmt = '%d')
    sys.exit()
    '''

    time_start = rec_start * sample_epoch * smpls_per_rec
    time_end = (rec_start + len(ticks)) * sample_epoch * smpls_per_rec
    seconds_t = np.linspace(time_start, time_end, num=len(ticks),
                            endpoint=False)
    seconds_p = np.linspace(time_start, time_end,
                            num=len(ticks)*smpls_per_rec,
                            endpoint=False)

    # Temperature period (usec)
    TP = (temp_raw/(TP_fctr)+TP_cnst)/(clock_freq)
    # Uncomment one of the lines below to hold temperature fixed
    # TP.fill(5.8224) #Fixed temp for 140344
    # TP.fill(5.7900) #Fixed temp for 140346
    # TP.fill(5.7875) #Fixed temp of +3.06°C for 140339
    # TP.fill(5.753) #Fixed temp of +2.27°C for 140338
    # TP.fill(5.8475) #Fixed temp of +2.55°C for 136309

    # Pressure period (usec)
    PP = (press_raw/(PP_fctr)+PP_cnst)/(clock_freq)

    # Apply smoothing filter to temperature before calculating pressure.
    # This eliminates significant noise from the pressure values.
    if tmptr_smth_fctr >= 5:
        print('Applying temperature smoothing filter.')
        TP = sig.savgol_filter(TP, tmptr_smth_fctr, 3, axis=0, mode='mirror')

    # Calculate temperature array
    print('Calculating temperatures and pressures.')
    Uv = TP-paros['U'][0]
    temperature = np.polyval(paros['Y'], Uv)
    # Upsample temperatures to match frequency of pressure samples by
    #  linear interpolation.
    Uv_expnd = np.interp(seconds_p, seconds_t, Uv)
    temperature_upsmpld = np.polyval(paros['Y'], Uv_expnd)
    
    # Calculate pressure array
    Cv = np.polyval(paros['C'], Uv_expnd)
    Dv = np.polyval(paros['D'], Uv_expnd)
    T0 = np.polyval(paros['T'], Uv_expnd)

    facts = 1 - (T0**2)/(PP**2)
    pressure = Cv * facts * (1 - Dv*facts)  # pressure in PSIA
    pressure = pressure * press_conv_fctr  # Convert pressure units

    # Decimate results
    # To produce sensible decimation results when ftype='iir',
    #    ensure n=5 and iterate using downsampling factor (q) <= 13.
    pressure_dcmtd = []
    if decmt_intvl > 0:
        print(f'Decimating results to {decmt_intvl} second epochs by '
              f'iteration, using factors;', end='', flush=True)
        pressure_dcmtd = pressure
        temperature_dcmtd = temperature_upsmpld
        seconds_dcmtd = seconds_p
        # First decimate to whole seconds
        decmt_intvl_pre = int(1/sample_epoch)
        first = True
        while decmt_intvl_pre > 1:
            if decmt_intvl_pre % 5 == 0:
                decmt_fctr = 5
            else:
                decmt_fctr = decmt_intvl_pre
            if not first:
                print(' :', end='')
            first = False
            print(f' {decmt_fctr}', end='', flush=True)
            pressure_dcmtd = sig.decimate(pressure_dcmtd, decmt_fctr,
                                          n=5, ftype='iir')
            temperature_dcmtd = sig.decimate(temperature_dcmtd, decmt_fctr,
                                          n=5, ftype='iir')
            seconds_dcmtd = seconds_dcmtd[::decmt_fctr]
            decmt_intvl_pre = decmt_intvl_pre // decmt_fctr
        # Now decimate to number of whole seconds requested
        while decmt_intvl > 1:
            if decmt_intvl % 10 == 0:
                decmt_fctr = 10
            elif decmt_intvl % 6 == 0:
                decmt_fctr = 6
            elif decmt_intvl % 5 == 0:
                decmt_fctr = 5
            elif decmt_intvl < 10 and decmt_intvl % 1 == 0:
                decmt_fctr = decmt_intvl
            else:
                sys.exit('\nDecimation failed! The interval specified must be '
                         'a single digit number of minutes or seconds, or be '
                         'divisible by 5 or 10.')
            print(f' : {decmt_fctr}', end='', flush=True)
            pressure_dcmtd = sig.decimate(pressure_dcmtd, decmt_fctr,
                                          n=5, ftype='iir')
            temperature_dcmtd = sig.decimate(temperature_dcmtd, decmt_fctr,
                                          n=5, ftype='iir')
            seconds_dcmtd = seconds_dcmtd[::decmt_fctr]
            decmt_intvl = decmt_intvl // decmt_fctr
        print()

#    print (ticks)
#    print(temperature)
#    print(pressure)
#    [print(press,end=', ') for press in pressure]

    if len(pressure_dcmtd) != 0:
        time = seconds_dcmtd
        pressure_out = pressure_dcmtd
        temperature_out = temperature_dcmtd
    else:
        time = seconds_p
        pressure_out = pressure
        temperature_out = temperature_upsmpld
    time_p = seconds_p
    time_t = seconds_t

    # Convert timestamp from seconds to datetime if specified.
    xlabel = 'Seconds'
    if xaxis_format == 'd':
        print('Calculating a timestamp array from the seconds array.')
        time = [clk_start_dt + dt.timedelta
                (seconds=sec) for sec in time]
        if len(pressure_dcmtd) != 0:
            time_p = [clk_start_dt + dt.timedelta
                    (seconds=sec) for sec in time_p]
            time_t = [clk_start_dt + dt.timedelta
                    (seconds=sec) for sec in time_t]
        xlabel = 'Date-Time'

    results = np.column_stack((time, pressure_out, temperature_out))

    # Save results to file
    if out_filename:
        print('Saving results to file.')
        np.savetxt(out_filename, results, fmt='%s,%f,%f')

    # Set min-max values for plot Y-axes
    print('Generating plot.')
    p_min = np.min(pressure)
    p_max = np.max(pressure)
    p_range = p_max - p_min
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
    p_min = p_min - p_min % intvl - intvl
    p_max = p_max - p_max % intvl + 2*intvl

    t_min = np.min(temperature)
    t_max = np.max(temperature)
    t_range = t_max-t_min
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
    t_min = t_min - t_min % intvl - intvl
    t_max = t_max - t_max % intvl + 2*intvl

    # Plot Results
    fig, ax1 = plt.subplots()
    color = 'red'
    plt.ylim(p_min, p_max)
    if len(pressure_dcmtd) != 0:
        ax1.plot(time_p, pressure, color='pink', marker='.',
                 markersize=1.0, linestyle='')
    ax1.plot(results[:, 0], results[:, 1], color=color, marker='.',
             markersize=2.0, linestyle='solid', linewidth=0.5)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel('Pressure (Pa)', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(axis='x')

    ax2 = ax1.twinx()
    color = 'blue'
    plt.ylim(t_min, t_max)
    ax2.set_ylabel('Temperature (°C)', color=color)
    if len(pressure_dcmtd) != 0:
        ax2.plot(time_t, temperature, color='lightblue', marker='.',
                 markersize=1.0, linestyle='')
    ax2.plot(results[:, 0], results[:, 2], color=color, marker='.',
             markersize=2.0, linestyle='solid', linewidth=0.5)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.suptitle(f'{apg_filename}\nStart: {wndw_start_dt}  -  '
                 f'End: {wndw_end_dt}')
    fig.tight_layout()
    # Uncomment the line below if timestamp format is datetime.
    # plt.gcf().autofmt_xdate()
    plt.show()

    return


###########################################################################
def extractrecords(apg_filename, start_byte, nrecs_want, rec_len, rec_fmt):
    '''
    Extracts binary records from a raw APG data logger file.
    '''
    print('Extracting raw records:', end='')
    with open(apg_filename, 'rb') as apgfile:
        apgfile.seek(start_byte, os.SEEK_CUR)
        records = []
        for i in range(0, nrecs_want):
            record_bin = apgfile.read(rec_len)

            # Print record as a string of Hex and/or Bin values.
            # Uncomment below for troubleshooting.
            '''
            hex_str = ''
            bin_str = ''
            for ch in record_bin:
                #hex_str += hex(ch)+' '
                bin_str += f'{ch:08b}'
            #print(hex_str)
            cum_rec_fmt = np.cumsum(list(map(abs,rec_fmt)))
            new_bin_str = ''
            for i in range( len(bin_str) ):
                if i in cum_rec_fmt:
                    new_bin_str = new_bin_str + ' '
                new_bin_str = new_bin_str + bin_str[i]
            print(new_bin_str)
            '''

            # Split record into array of ints defined as groups of bits
            # by rec_fmt.
            record_int = int.from_bytes(record_bin, byteorder='big',
                                        signed=False)
            record = []
            for signed_bit_len in reversed(rec_fmt):
                bit_len = int(abs(signed_bit_len))
                # Read right most bit_len bits
                field = record_int & (2**bit_len-1)
                # Check for sign bit and convert as a 2s-compliment negative.
                if ((signed_bit_len < 0) and (field & (2**(bit_len-1)))):
                    full_bit_len = bit_len + (bit_len % 8)
                    field = field | ((2**full_bit_len-1) ^ (2**bit_len-1))
                    field = field.to_bytes(full_bit_len//8, byteorder='big',
                                           signed=False)
                    field = int.from_bytes(field, byteorder='big', signed=True)
                record.append(field)
                # Shift to right bit_len bits
                record_int = record_int >> (bit_len)

            records.append(record)
            # Print a "." every 10000 records to indicate script is running.
            if i % 10000 == 0:
                print('.', end='', flush=True)
        print('\n')
    return(np.array(records))


#####################################################################
def clockdrift(clk_start_dt, gpssync_dt, sync_tick_count, tick_bits):
    '''
    Calculates the clock drift using one of two methods.
    The expected number of tick counts between clk_start_dt and gpssync_dt
    will be calculated. The size of the tick record in tick_bits is used to
    determine when the tick count 'rolls over' to zero.
    If sync_tick_count is provided, the difference between this and the
    expected value gives the drift in ticks (milliseconds).
    If sync_tick_count is not present, it is assumed that a fixed frequency has
    been injected into the raw data precisely at gpssync_dt. The precise tick
    count when this frequency starts is detected in the data and this value
    is used in place of sync_tick_count.
    '''
    if (sync_tick_count is None):
        sys.exit('The code for time sync by frequency injection has not been '
                 'written yet.')

    millisecs_logged = ((gpssync_dt - clk_start_dt).total_seconds()) * 1000
    # Number of ticks until rollover and restart at zero.
    tick_rollover = 2**tick_bits  # 1 tick = 1 millisecond
    # APG logger clock offset relative to GPS time (positive = APG > GPS)
    clk_drift_at_end = int(sync_tick_count -
                           (millisecs_logged % tick_rollover))

    return(clk_drift_at_end)


##########################
if __name__ == '__main__':
    main()
