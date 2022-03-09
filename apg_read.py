#!/home/nevillep/miniconda3/bin/python

# Version 20220307
# by Neville Palmer, GNS Science
# 2018/11/17 Start development

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
    apg_ini = './ParosAPG.ini'
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
    parser.add_argument("-a", "--apgini",
                        help='Full path and filename for Paros APG '
                        f'configuration settings. Default: "{apg_ini}"',
                        default=apg_ini)
    parser.add_argument("-s", "--snapg",
                        help='Serial number of the Paroscientific APG used. '
                        'This must correspond to the serial number of an '
                        'entry in the apgini file.',
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
                        help=f'Required sample interval in seconds for '
                        f'pressure decimation. Zero for no decimation. '
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
    parser.add_argument("-f", "--fmttime",
                        help='Specify the format to be used for presenting time '
                        'in outputs and plots, to be displayed as either '
                        '(s)econds or (d)ate-time.',
                        choices=['s','d'],
                        default='s')
    parser.add_argument("-p", "--plot",
                        help='Generate and display a timeline plot.',
                        action='store_true')
    parser.add_argument("-r", "--rawbinary",
                        help='Output raw binary data to stdout.',
                        action='store_true')
    args = parser.parse_args()

    # Translate argparse parameters (except for window times).
    apg_filename = os.path.normpath(args.infile)
    apg_ini = os.path.normpath(args.apgini)
    apg_sn = args.snapg.strip()
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
    time_format = args.fmttime.strip()
    plot_flag = args.plot
    rawbin_flag = args.rawbinary

    # Read Paros transducer coefficients into a dict of lists from ini file.
    paros_coefs = ('U', 'Y', 'C', 'D', 'T')
    paros_cfg = configparser.ConfigParser()
    paros_cfg.read(apg_ini)
    paros = {}
    for coef in paros_coefs:
        paros[coef] = (paros_cfg.get(apg_sn, coef).split(','))
        paros[coef] = [float(x) for x in paros[coef][::-1]]
    paros['Y'].append(0.0)

    # Read APG logger configuration parameters from ini file.
    logger_cfg = configparser.ConfigParser()
    logger_cfg.read(logger_ini)
    logger = {}
    logger['head_len'] = logger_cfg.getint(logger_version, 'head_len')
    logger['rec_len']  = logger_cfg.getint(logger_version, 'rec_len')
    logger['smpls_per_rec']  = logger_cfg.getint(logger_version, 'smpls_per_rec')
    logger['sample_epoch']  = logger_cfg.getfloat(logger_version, 'epoch')
    logger['record_epoch']  = logger['sample_epoch']  * logger['smpls_per_rec'] 
    logger['clock_freq']  = logger_cfg.getint(logger_version, 'clock_freq')
    logger['TP_fctr']  = evaluate(logger_cfg.get(logger_version, 'TP_fctr')).item()
    logger['TP_cnst']  = logger_cfg.getfloat(logger_version, 'TP_cnst')
    logger['PP_fctr']  = evaluate(logger_cfg.get(logger_version, 'PP_fctr')).item()
    logger['PP_cnst']  = logger_cfg.getfloat(logger_version, 'PP_cnst')
    logger['timing']  = logger_cfg.get(logger_version, 'timing')

    logger['rec_fmt']  = (logger_cfg.get(logger_version, 'rec_fmt').split(','))
    logger['rec_fmt']  = tuple([float(x) for x in logger['rec_fmt'] ])
    fmt_field = {}
    fmt_field['tic'] = logger_cfg.getint(logger_version, 'tic_field')
    fmt_field['tptr'] = logger_cfg.getint(logger_version, 'temperature_field')
    fmt_field['pcore'] = logger_cfg.getint(logger_version, 'pcore_field')
    fmt_field['pn'] = logger_cfg.getint(logger_version, 'pn_field')
    if (abs(fmt_field['pcore']-fmt_field['pn']+1) != logger['smpls_per_rec'] ):
        sys.exit(f"The number of samples per record ({logger['smpls_per_rec'] }) "
                 f"for logger {logger_version}, does "
                 f"not match the number of Pcore and Pn fields in the "
                 f"record format (Pcore through Pn inclusive = "
                 f"{fmt_field['pcore']-fmt_field['pn']+1}), "
                 f"as provided in file {logger_ini}.")
    rec_fmt_bits = int(sum(map(abs, logger['rec_fmt'] )))
    if (rec_fmt_bits != (logger['rec_len'] *8)):
        sys.exit(f"The total number of bits ({rec_fmt_bits}) given by the "
                 f"record format {logger['rec_fmt'] } "
                 f"does not equal the record length ({logger['rec_len'] } bytes x 8), "
                 f"as provided in the file {logger_ini}.")
    logger['fmt_field'] = fmt_field
    logger['tic_bit_len'] = logger['rec_fmt'] [fmt_field['tic']]

    # Assign names to column numbers of raw data array.
    # Note that raw data array columns are reverse order to raw binary.
    last_field = len(logger['rec_fmt'] ) - 1
    pcore_col = last_field - fmt_field['pcore']
    pn_col = last_field - fmt_field['pn']
    tptr_col = last_field - fmt_field['tptr']


    # Calculate duration and end time of raw file.
    try:
        fsize = os.path.getsize(apg_filename)  # file size in bytes
    except(FileNotFoundError):
        sys.exit(f'Raw APG file "{apg_filename}" does not exist.')
    print(f'Filesize = {fsize} bytes')
    nrecs = int((fsize - logger['head_len'] ) / logger['rec_len'] )-1
    print(f'Number of records = {nrecs:d}')
    file_duration_days = nrecs*logger['record_epoch'] /3600000/24
    print(f'File duration = {file_duration_days} days')
    clk_end_dt = clk_start_dt + dt.timedelta(days=file_duration_days)
    print(f'Time of last sample = {clk_end_dt}')

    # Clock drift
    if (args.gpssynctime is not None):
        drift = clockdrift(apg_filename, logger, clk_start_dt, gpssync_dt, sync_tick_count)
        print(f'Clock drift in milliseconds (logged time minus actual '
              f'GPS time) = {drift}')

    # Calculate window for extracting data.
    if args.beginwndw is None:
        wndw_begin_dt = clk_start_dt
    else:
        wndw_begin = re.sub('[-: _/]', '_', args.beginwndw)
        wndw_begin_dt = dt.datetime.strptime(wndw_begin,
                                             '%Y_%m_%d_%H_%M_%S.%f')
    if args.endwndw is None:
        wndw_end_dt = clk_end_dt
    else:
        wndw_end = re.sub('[-: _/]', '_', args.endwndw)
        wndw_end_dt = dt.datetime.strptime(wndw_end, '%Y_%m_%d_%H_%M_%S.%f')
    wndw_begin_secs = ((wndw_begin_dt - clk_start_dt).total_seconds())
    wndw_begin_days = wndw_begin_secs / (24*3600)
    print(f'Data window beginning as date-time = {wndw_begin_dt}')
    print(f'Data window beginning in days since clock start: {wndw_begin_days}')
    wndw_len_secs = ((wndw_end_dt - wndw_begin_dt).total_seconds())
    wndw_len_days = wndw_len_secs / (24*3600)
    print(f'Data window length in days: {wndw_len_days}')
    rec_begin = int(wndw_begin_secs*1000 / logger['record_epoch'] )
    print(f'First record of data window: {rec_begin}')
    nrecs_want = int(wndw_len_secs*1000 / logger['record_epoch'] )+1

    # Extract records from APG file for the time window specified.
    print('Extracting raw records:\n', end='')
    records = extractrecords(apg_filename, logger, nrecs_want, rec_begin,
                             rawbin_flag)

    # Save raw records to file as integers with tick rollover removed.
    ''' Comment this line for troubleshooting.
    np.savetxt('raw_records_no-rollover.txt', records, fmt='%d',
            header='',
            comments='')
    #'''


    # Create an array for each raw observable (pressure, temperature, ticks)
    if pn_col > pcore_col:
        press_raw = records[:, pcore_col:pn_col+1]
    else:
        press_raw = records[:, pcore_col:pn_col-1:-1]
    press_raw = np.cumsum(press_raw, axis=1)
    press_raw = press_raw.reshape((nrecs_want*logger['smpls_per_rec'] ))
    temp_raw = (records[:, tptr_col])
    ticks = (records[:, last_field - fmt_field['tic']])  # ticks are equivalent to milliseconds


    actual_end_tick = ticks[-1] + (logger['smpls_per_rec']  - 1) * logger['sample_epoch'] 
    actual_begin_tick = ticks[0]


    # Write a summary file showing out of sync tic counts and exit.
    ''' Comment this line for troubleshooting.
    millisecs_t = np.linspace(nominal_begin_tick, nominal_end_tick,
                              num=nrecs_want,
                              endpoint=True)
    millisecs_p = np.linspace(nominal_begin_tick,
                              nominal_end_tick + logger['record_epoch'] ,
                              num=nrecs_want*logger['smpls_per_rec'] ,
                              endpoint=False)
    time_diff = ticks - millisecs_t
    ticks = np.column_stack((ticks,millisecs_t,time_diff))
    summary_ticks = [ticks[0]]
    for i in range(1,len(time_diff)):
        if (time_diff[i] != time_diff[i-1]):
            #print(f'{time_diff[i]} : {time_diff[i-1]} : '
            #      f'{time_diff[i] != time_diff[i-1]}')
            summary_ticks.append(ticks[i])
    np.savetxt('summary_ticks.txt', summary_ticks, fmt='%d',
                header='actual,nominal,difference',
                comments='')
    np.savetxt('ticks.txt', ticks, fmt='%d',
                header='actual,nominal,difference',
                comments='')
    #sys.exit()
    #'''

    # If the nominal tick count and actual recorded tick count are not precisely
    # aligned then use linear interpolation to generate a precisely periodic
    # record of raw temperature and pressure values.
    nominal_begin_tick = rec_begin * logger['record_epoch'] 
    nominal_end_tick = (rec_begin + nrecs_want -1) * logger['record_epoch'] 

    # print(f'nom_begin: {nominal_begin_tick}\nact_begin: {actual_begin_tick}\n'
    #       f'nom_end: {nominal_end_tick}\nact_end: {actual_end_tick}\n')
    if (actual_begin_tick == nominal_begin_tick 
        and actual_end_tick == nominal_end_tick):
        millisecs_t = np.linspace(nominal_begin_tick,
                                  nominal_end_tick,
                                  num=nrecs_want,
                                  endpoint=True)
        millisecs_p = np.linspace(nominal_begin_tick,
                                  nominal_end_tick + logger['record_epoch'] ,
                                  num=nrecs_want*logger['smpls_per_rec'] ,
                                  endpoint=False)
    else:
        # Determine first tick count to achieve fixed period epochs.
        final_begin_tick = (actual_begin_tick + (logger['record_epoch']  - 
                            actual_begin_tick % logger['record_epoch'] ))
        # Determine final tick count to achieve fixed period epochs.
        final_end_tick = (actual_end_tick - (actual_end_tick % logger['record_epoch'] ))
        epoch_count = int((final_end_tick-final_begin_tick) / logger['record_epoch'] ) + 1
        millisecs_t  = np.linspace(final_begin_tick,
                                   final_end_tick,
                                   num=epoch_count,
                                   endpoint=True)
        millisecs_p  = np.linspace(final_begin_tick,
                                   final_end_tick + logger['record_epoch'] ,
                                   num=epoch_count*logger['smpls_per_rec'] ,
                                   endpoint=False)
        ticks_p = np.linspace(actual_begin_tick,
                                  actual_end_tick + logger['record_epoch'] ,
                                  num=nrecs_want*logger['smpls_per_rec'] ,
                                  endpoint=False)
        # Interpolate to generate fixed period observation epochs.
        temp_raw = np.interp(millisecs_t, ticks, temp_raw)
        press_raw = np.interp(millisecs_p, ticks_p, press_raw)

    # Temperature period (usec)
    TP = (temp_raw/(logger['TP_fctr'] )+logger['TP_cnst'] )/(logger['clock_freq'] )
    # Uncomment one of the lines below to hold temperature fixed
    # TP.fill(5.8224) #Fixed temp for 140344
    # TP.fill(5.7900) #Fixed temp for 140346
    # TP.fill(5.7875) #Fixed temp of +3.06°C for 140339
    # TP.fill(5.753) #Fixed temp of +2.27°C for 140338
    # TP.fill(5.8475) #Fixed temp of +2.55°C for 136309

    # Pressure period (usec)
    PP = (press_raw/(logger['PP_fctr'] )+logger['PP_cnst'] )/(logger['clock_freq'] )

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
    Uv_expnd = np.interp(millisecs_p, millisecs_t, Uv)
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
        millisecs_dcmtd = millisecs_p
        # First decimate to whole seconds
        sample_freq = int(1000/logger['sample_epoch'] )
        decmt_intvl_pre = sample_freq
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
            millisecs_dcmtd = millisecs_dcmtd[::decmt_fctr]
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
            millisecs_dcmtd = millisecs_dcmtd[::decmt_fctr]
            decmt_intvl = decmt_intvl // decmt_fctr
        print()

    if len(pressure_dcmtd) != 0:
        time = millisecs_dcmtd / 1000
        pressure_out = pressure_dcmtd
        temperature_out = temperature_dcmtd
    else:
        time = millisecs_p / 1000
        pressure_out = pressure
        temperature_out = temperature_upsmpld
    time_p = millisecs_p / 1000
    time_t = millisecs_t / 1000

    # Convert timestamp from seconds to datetime if specified.
    xlabel = 'Seconds'
    if time_format == 'd':
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

    if plot_flag:
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
        #'''
        if len(pressure_dcmtd) != 0:
            ax1.plot(time_p, pressure, color='pink', marker='.',
                    markersize=1.0, linestyle='')
        #'''
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
        #'''
        if len(pressure_dcmtd) != 0:
            ax2.plot(time_t, temperature, color='lightblue', marker='.',
                    markersize=1.0, linestyle='')
        #'''
        ax2.plot(results[:, 0], results[:, 2], color=color, marker='.',
                markersize=2.0, linestyle='solid', linewidth=0.5)
        ax2.tick_params(axis='y', labelcolor=color)
        fig.suptitle(f'{apg_filename}\nBegin: {wndw_begin_dt}  -  '
                    f'End: {wndw_end_dt}')
        fig.tight_layout()
        # Uncomment the line below if timestamp format is datetime.
        # plt.gcf().autofmt_xdate()
        plt.show()

    return


###########################################################################
def extractrecords(apg_filename, logger, nrecs_want, rec_begin, rawbin_flag):
    '''
    Extracts binary records from a raw APG data logger file.
    '''
    with open(apg_filename, 'rb') as apgfile:
        begin_byte = logger['head_len']  + rec_begin * logger['rec_len']
        apgfile.seek(begin_byte , os.SEEK_CUR)
        records = []
        for i in range(0, nrecs_want):
            record_bin = apgfile.read(logger['rec_len'] )

            # Print record as a string of Hex and/or Bin values to stdout.
            if rawbin_flag:
                hex_str = ''
                bin_str = ''
                for ch in record_bin:
                    #hex_str += hex(ch)+' '
                    bin_str += f'{ch:08b}'
                #print(hex_str)
                cum_rec_fmt = np.cumsum(list(map(abs,logger['rec_fmt'] )))
                new_bin_str = ''
                for i in range( len(bin_str) ):
                    if i in cum_rec_fmt:
                        new_bin_str = new_bin_str + ' '
                    new_bin_str = new_bin_str + bin_str[i]
                print(new_bin_str)

            # Split record into array of ints defined as groups of bits
            # by logger['rec_fmt'] .
            record_int = int.from_bytes(record_bin, byteorder='big',
                                        signed=False)
            record = []
            for signed_bit_len in reversed(logger['rec_fmt'] ):
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
        print()

        records = np.array(records)

        # Save raw records to file as integers without tick rollover removed.
        ''' Comment this line for troubleshooting.
        np.savetxt('raw_records_rollover.txt', records, fmt='%d',
                header='',
                comments='')
        #'''


        # Shift the tick count if necessary, so that it relates to the first sample
        # in each record (instead of the last).
        if logger['timing']  == 'first':
            first_tic = 0
        elif logger['timing']  == 'last':
            first_tic = int((logger['smpls_per_rec']  - 1) * logger['sample_epoch'] )
        last_field = len(logger['rec_fmt'] ) - 1
        # ticks are equivalent to milliseconds
        ticks = (records[:, last_field - logger['fmt_field']['tic']] - first_tic)


        # Remove tick count rollovers and make actual ticks continuously increasing.
        nominal_begin_tick = rec_begin * logger['record_epoch'] 
        nominal_end_tick = (rec_begin + nrecs_want -1) * logger['record_epoch'] 
        # print(f'Tick field length (in bits): {tic_field_len}')
        rollover_period = 2 ** logger['tic_bit_len']  # in millisec
        # print(f'Rollover length (in millisec/ticks): {rollover_period}')
        # The number of rollovers prior to the beginning of the specified data window.
        nom_rollovers_begin = int(nominal_begin_tick / rollover_period)
        nom_rollover_balance = nominal_begin_tick % rollover_period
        # Does the nominal rollover count align with the actual count.
        # (Within 1e7 millisec or approx 166 minutes.)
        if abs(ticks[0] - nom_rollover_balance) < 1e7:
            actl_rollovers_begin = nom_rollovers_begin
        elif ticks[0] - nom_rollover_balance > 1e7:
            actl_rollovers_begin = nom_rollovers_begin - 1
        else:
            actl_rollovers_begin = nom_rollovers_begin + 1

        # {rollovers} contains the index of the first record after each rollover.
        rollovers = np.where(ticks[:-1] > ticks[1:])[0] + 1
        #print(f'Index values of rollovers within ticks array: {rollovers}')
        cumtv_rollovers = actl_rollovers_begin
        #print(f'Count of cumulative rollovers at beginning of window: {cumtv_rollovers}')
        if cumtv_rollovers != 0:
            if rollovers.size == 0:
                ticks = ticks + rollover_period * cumtv_rollovers
            else:
                ticks[0:rollovers[0]] = (ticks[0:rollovers[0]] 
                                        + rollover_period * cumtv_rollovers)

        for idx, rollover in np.ndenumerate(rollovers):
            if rollover == rollovers[-1]:
                nxt_rollover = nrecs_want
            else:
                nxt_rollover = rollovers[idx[0]+1]
            #print(rollover, nxt_rollover)
            # If the tick count does not rollover cleanly and takes more than one
            # record to reset to zero, then for first record after the rollover 
            # calc as the previous cumulative tick count plus a std record period.
            if nxt_rollover - rollover == 1:
                ticks[rollover] = ticks[rollover-1] + logger['record_epoch'] 
            else:
                cumtv_rollovers = cumtv_rollovers + 1
                ticks[rollover:nxt_rollover] = (ticks[rollover:nxt_rollover] 
                                            + rollover_period * cumtv_rollovers)

    records[:, -1] = ticks

    return(records)


#####################################################################
def clockdrift(apg_filename, logger, clk_start_dt, gpssync_dt, sync_tick_count):
    '''
    Calculates the clock drift using one of two methods.
    The expected number of tick counts between clk_start_dt and gpssync_dt
    will be calculated. The size of the tick record in logger['tic_bit_len'] is
    used to determine when the tick count 'rolls over' to zero.
    If sync_tick_count is provided, the difference between this and the
    expected value gives the drift in ticks (milliseconds).
    If sync_tick_count is not present, it is assumed that a fixed frequency has
    been injected into the raw data precisely at gpssync_dt. The precise tick
    count when this frequency starts is detected in the data and this value
    is used in place of sync_tick_count.
    '''
    
    millisecs_logged = ((gpssync_dt - clk_start_dt).total_seconds()) * 1000

    if (sync_tick_count is None):
        # Assign names to column numbers of raw data array.
        # Note that raw data array columns are reverse order to raw binary.
        last_field = len(logger['rec_fmt'] ) - 1
        tick_col = last_field - logger['fmt_field']['tic']
        pcore_col = last_field - logger['fmt_field']['pcore']
        pn_col = last_field - logger['fmt_field']['pn']
        tptr_col = last_field - logger['fmt_field']['tptr']

        # Window for identifying sync_tick_count is +/-5 minutes long.
        sync_wndw_secs = 5 * 60
        # GPS sync time (gpssync_dt) is mid point of  window for sync.
        wndw_begin_secs = ((gpssync_dt - clk_start_dt).total_seconds())
        wndw_begin_secs = wndw_begin_secs - (sync_wndw_secs)
        rec_begin = int(wndw_begin_secs*1000 / logger['record_epoch'] )
        nrecs_want = int(sync_wndw_secs*2000 / logger['record_epoch'] )+1
        
        sync_records = extractrecords(apg_filename, logger, nrecs_want, rec_begin,
                             False)
        # Save raw records to file as integers with tick rollover removed.
        #''' Comment this line for troubleshooting.
        np.savetxt('raw_sync_records.txt', sync_records, fmt='%d',
                header='',
                comments='')
        #'''
        
        # Identify the start of the record block where the pressure values 
        # start changing again (ie This is where frequency injection for time 
        # sync occurs).
        # Identify all consecutive row pairs where p_core changes.
        pcore_diff = np.diff(sync_records[:,pcore_col])
        pcore_diff_row = (np.where(pcore_diff!=0)[0])+1
        # Select the final instance where p_core starts changing
        # This exculdes any single noise values occuring before actual frequency
        # injection.
        diff_row_increments = np.diff(pcore_diff_row)
        x = np.where(diff_row_increments>1)[0]+1
        x = np.insert(x,0,0)
        poss_sync_row = pcore_diff_row[x]-1
        print("poss_sync_row: ",poss_sync_row)

        # For each poss_sync_row check:
        #   - immed prev row has all Pn values as zero,
        #      (Indicates two consecutive p_core values to be identical although
        #       surrounding values continue changing)
        #   - immed next row does not have all Pn values as zero.
        #       (Indicates noise value occuring before actual frequency injection.)
        # It is possible for two consecutive p_core values
        # to be identical although surrounding values continue changing.
        for n in range(np.size(poss_sync_row)-1, -1, -1):
            sync_row = poss_sync_row[n]
            prev_sync_row = sync_row - 1
            prev_block = sync_records[prev_sync_row,:]
            next_sync_row = sync_row + 1
            next_block = sync_records[next_sync_row,:]

            if pn_col > pcore_col:
                prev_pn = prev_block[pcore_col+1:pn_col+1]
                next_pn = next_block[pcore_col+1:pn_col+1]
            else:
                prev_pn = prev_block[pcore_col-1:pn_col-1:-1]
                next_pn = next_block[pcore_col-1:pn_col-1:-1]
            
            prev_nonzero = (np.where(prev_pn!=0)[0])
            next_nonzero = (np.where(next_pn!=0)[0])
            if not prev_nonzero.any() and next_nonzero.any():
                sync_row = poss_sync_row[n]
                break

        try:
            sync_block = sync_records[sync_row,:]
        except UnboundLocalError:
            sys.exit(f'Unable to determine clock drift.\n'
                f'The period {gpssync_dt} +/-{sync_wndw_secs} seconds does not include a \n'
                f'frequency injection for syncing to.')

        if pn_col > pcore_col:
            pn = sync_block[pcore_col+1:pn_col+1]
        else:
            pn = sync_block[pcore_col-1:pn_col-1:-1]

        nonzero = (np.where(pn!=0)[0])
        if nonzero.any():
            i = len(pn)+1 - nonzero.size
        else:
            i = len(pn)+1

        sync_tick_count = sync_block[tick_col] + (i * logger['sample_epoch'])
        sync_tick_count = int(sync_tick_count)

    else:
        # Number of ticks until rollover and restart at zero.
        tick_rollover = 2**logger['tic_bit_len']  # 1 tick = 1 millisecond
        millisecs_logged = millisecs_logged % tick_rollover

    # APG logger clock offset relative to GPS time (positive = APG > GPS)
    clk_drift_at_end = int(sync_tick_count - (millisecs_logged))

    return(clk_drift_at_end)


##########################
if __name__ == '__main__':
    main()
