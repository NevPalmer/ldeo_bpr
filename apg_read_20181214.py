#!/usr/bin/env python3

# Version 20181207
# by Neville Palmer, GNS Science
# 2018/11/17 Start development
# 2018/11/23 SeaScan board format added.
# 2018/11/24 Smoothing and decimation added.
# 2018/11/25 Added ini files for reading in transducer & logger parameters.
# 2018/11/26 Resolved issue with iir decimation (ensure q<13 and n=5).
#               Iterate for desired decimation.
# 2018/11/28 Add iteration routine for decimation.
# 2018/12/02 Add command line parameters.
# 2018/12/07 Add output tic count summary file.
# 2018/12/14 Add ability to calc clock drift for CSAC logger.

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
    logger_versions = ['CSAC2013','Seascan2018']
    clk_start =  '2000-01-01_00:00:00' # 'YYYY-MM-DD_hh:mm:ss'
    out_filename = './apgOut.txt'
    tmptr_smth_fctr = 1
    decmt_intvl = 0

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
    parser.add_argument("-i","--infile", 
        help=f'Full path and filename of raw APG input file.',
        required=True)
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
    parser.add_argument("-d","--decimate", 
        help=f'Sample interval in seconds for pressure decimation. Zero for no '
            f'decimation. Default: "{decmt_intvl}"',
        type=int,
        default=decmt_intvl)
    parser.add_argument("-t","--tempsmth", 
        help=f'Temperature smoothing factor (must be an odd integer). 10001 '
            f'gives strong smoothing. Default: "{tmptr_smth_fctr}"',
        type=int,
        default=tmptr_smth_fctr)
    parser.add_argument("-c","--clkstart", 
        help=f'Precise date and time when the logger clock was started. '
            f'Default: "{clk_start}"',
        default=clk_start)
    parser.add_argument("-b","--beginwndw", 
        help=f'Date and time to begin data extraction. Assumes beginning of '
            f'file if omitted. Format: "YYYY-MM-DD_hh:mm:ss.s"')
    parser.add_argument("-e","--endwndw", 
        help=f'Date and time to end data extraction. Assumes end of file if '
            f'omitted. Fotmat: "YYYY-MM-DD_hh:mm:ss.s"')
    parser.add_argument("-g","--gpssynctime", 
        help=f'Precise date and time from GPS clock for syncronising end time. '
            f'No clock drift adjustment is made if omitted. '
            f'Fotmat: "YYYY-DDD_hh:mm:ss"')
    parser.add_argument("-y","--synctickcount", 
        help=f'The hexidecimal tick count that corresponds to GPSSYNCTIME. '
            f'If GPSSYNCTIME is specified and SYNCTICKCOUNT is omitted, then '
            f'it is assumed that an artificial frequency was inserted precisely '
            f'at GPSSYNCTIME. This parameter is ignored if GPSSYNCTIME is omitted. '
            f'Fotmat: "0xHHHHHHHHHH"',
        type=lambda x: int(x,0))
    parser.add_argument("-o","--outfile", 
        help=f'Full path and filename for output file. Default: "{out_filename}"',
        default=out_filename)
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
    if (args.gpssynctime != None):
        gpssynctime = re.sub('[-: _/]', '_', args.gpssynctime)
        gpssync_dt = dt.datetime.strptime(gpssynctime, '%Y_%j_%H_%M_%S')
    sync_tick_count = args.synctickcount
    out_filename = os.path.normpath(args.outfile)

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
    sample_epoch = logger_cfg.getfloat(logger_version,'epoch')
    rec_len = logger_cfg.getint(logger_version,'rec_len')
    head_len = logger_cfg.getint(logger_version,'head_len')
    smpls_per_rec = logger_cfg.getint(logger_version,'smpls_per_rec')
    TP_fctr = evaluate(logger_cfg.get(logger_version,'TP_fctr')).item()
    TP_cnst = logger_cfg.getfloat(logger_version,'TP_cnst')
    PP_fctr = evaluate(logger_cfg.get(logger_version,'PP_fctr')).item()
    PP_cnst = logger_cfg.getfloat(logger_version,'PP_cnst')
    
    rec_fmt = (logger_cfg.get(logger_version,'rec_fmt').split(','))
    rec_fmt = tuple([float(x) for x in rec_fmt])
    fmt_fields = {}
    fmt_fields['tic'] = logger_cfg.getint(logger_version,'tic_field')
    fmt_fields['tptr'] = logger_cfg.getint(logger_version,'temperature_field')
    fmt_fields['pcore'] = logger_cfg.getint(logger_version,'pcore_field')
    fmt_fields['pn'] = logger_cfg.getint(logger_version,'pn_field')
    if (abs(fmt_fields['pcore']-fmt_fields['pn']+1) != smpls_per_rec):
        sys.exit(f"The number of samples per record ({smpls_per_rec}) "
                 f"for logger {logger_version}, does "
                 f"not match the number of Pcore and Pn fields in the "
                 f"record format (Pcore through Pn inclusive = "
                 f"{fmt_fields['pcore']-fmt_fields['pn']+1}), "
                 f"as provided in file {logger_ini}.")
    rec_fmt_bits = int(sum(map(abs,rec_fmt)))
    if (rec_fmt_bits != (rec_len*8)):
        sys.exit(f"The total number of bits ({rec_fmt_bits}) given by the "
                 f"record format {rec_fmt} "
                 f"does not equal the record length ({rec_len} bytes x 8), "
                 f"as provided in the file {logger_ini}.")

    # Calculate duration and end time of raw file.
    try:
        fsize = os.path.getsize(apg_filename) # file size in bytes
    except(FileNotFoundError):
        sys.exit(f'Raw APG file "{apg_filename}" does not exist.')
    print(f'Filesize = {fsize} bytes')
    nrecs = int((fsize - head_len) / rec_len)
    print(f'Number of records = {nrecs:d}')
    file_duration_days = nrecs*smpls_per_rec*sample_epoch/3600/24
    print(f'File duration = {file_duration_days} days')
    clk_end_dt = clk_start_dt + dt.timedelta(days=file_duration_days)
    print(f'Time of last sample = {clk_end_dt}')
    
    # Clock drift
    if (args.gpssynctime != None):
        drift = clockdrift(clk_start_dt, gpssync_dt, sync_tick_count, rec_fmt[fmt_fields['tic']])
        print(f'Clock drift in milliseconds = {drift}')


    # Calculate window for extracting data.
    if args.beginwndw == None:
        wndw_start_dt = clk_start_dt
    else:
        wndw_start = re.sub('[-: _/]', '_', args.beginwndw)
        wndw_start_dt = dt.datetime.strptime(wndw_start, '%Y_%m_%d_%H_%M_%S.%f')
    if args.endwndw == None:
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
    records = extractrecords(apg_filename, start_byte, nrecs_want, rec_len, rec_fmt)

    # Assign extracted records to field names.    
    last_field = len(rec_fmt) - 1
    pcore = last_field - fmt_fields['pcore']
    pn = last_field - fmt_fields['pn']
    tptr = last_field - fmt_fields['tptr']
    tic = last_field - fmt_fields['tic']
    if pn > pcore:
        press_raw = records[:,pcore:pn+1]
    else:
        press_raw = records[:,pcore:pn-1:-1]
    press_raw = np.cumsum(press_raw,axis=1)
    press_raw = press_raw.reshape((nrecs_want*smpls_per_rec))
    temp_raw = (records[:,tptr:tptr+1])
    ticks = (records[:,tic])

    # Write a summary file showing out of sync tic counts and exit.
    # Uncomment below for troubleshooting.
    '''
    nom_time = np.arange(64,64+72*np.alen(ticks),72,dtype=int)
    time_diff = ticks - nom_time
    ticks = np.column_stack((ticks,nom_time,time_diff))
    summary_ticks = [ticks[0]]
    for i in range(1,len(time_diff)):
        if (time_diff[i] != time_diff[i-1]):
            print(f'{time_diff[i]} : {time_diff[i-1]} : {time_diff[i] != time_diff[i-1]}')
            summary_ticks.append(ticks[i])
    np.savetxt('summary_ticks.txt', summary_ticks, fmt = '%d')
    np.savetxt('ticks.txt', ticks, fmt = '%d')
    sys.exit()
    '''
    
    time_start = rec_start * sample_epoch * smpls_per_rec
    time_end = (rec_start + len(ticks)) * sample_epoch * smpls_per_rec
    seconds_t = np.linspace(time_start,time_end,num=len(ticks),endpoint=False)
    seconds_p = np.linspace(time_start,time_end,num=len(ticks)*smpls_per_rec,endpoint=False)

    # Temperature period (usec)
    TP = (temp_raw/(TP_fctr)+TP_cnst)/(clock_freq)
    # Pressure period (usec)
    PP = (press_raw/(PP_fctr)+PP_cnst)/(clock_freq)

    # Apply smoothing filter to temperature before calculating pressure.
    # This eliminates significant noise from the pressure values.
    if tmptr_smth_fctr >= 5:
        TP = sig.savgol_filter(TP, tmptr_smth_fctr, 3, axis=0, mode='mirror')
        
    # Create temperature
    Uv = TP-paros['U'][0]
    temperature = np.polyval(paros['Y'], Uv)
    
    if decmt_intvl != 0:
        temperature = temperature.reshape((nrecs_want))
        temperature_resmpld = sig.resample_poly(
                temperature, 
                smpls_per_rec,
                decmt_intvl * int(1 / sample_epoch))
    else:
        temperature_resmpld = temperature * np.ones((1, smpls_per_rec), dtype=int)
        temperature_resmpld = temperature_resmpld.reshape((nrecs_want * smpls_per_rec))
        '''
        temperature_resmpld = sig.resample_poly(
                temperature, 
                smpls_per_rec,
                1)'''
        
    Uv_expnd = Uv * np.ones((1, smpls_per_rec), dtype=int)
    Uv_expnd = Uv_expnd.reshape((nrecs_want * smpls_per_rec))

    # Create pressure
    Cv = np.polyval(paros['C'], Uv_expnd)
    Dv = np.polyval(paros['D'], Uv_expnd)
    T0 = np.polyval(paros['T'], Uv_expnd)
    
    facts = 1 - (T0**2)/(PP**2) 
    pressure = Cv * facts *(1-Dv*facts) #pressure in PSIA
    pressure = pressure * press_conv_fctr #Convert pressure units
    
    # To produce sensible decimation results when ftype='iir', 
    #    ensure n=5 and iterate using downsampling factor (q) <= 13.
    pressure_dcmtd = []
    if decmt_intvl > 0:
        pressure_dcmtd = pressure
        seconds_dcmtd = seconds_p
        # First decimate to whole seconds
        decmt_intvl_pre = int(1/sample_epoch)
        while decmt_intvl_pre > 1:
            if decmt_intvl_pre % 5 == 0:
                decmt_fctr = 5
            else:
                decmt_fctr = decmt_intvl_pre
            pressure_dcmtd = sig.decimate(pressure_dcmtd,decmt_fctr,n=5,ftype='iir')
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
                sys.exit('The pressure decimation factor must be a single digit '
                    'number of minutes or seconds, or be divisible by 5 or 10.')
            pressure_dcmtd = sig.decimate(pressure_dcmtd,decmt_fctr,n=5,ftype='iir')
            seconds_dcmtd = seconds_dcmtd[::decmt_fctr]
            decmt_intvl = decmt_intvl // decmt_fctr
    
#    print (ticks)
#    print(temperature)
#    print(pressure)
#    [print(press,end=', ') for press in pressure]
    if len(pressure_dcmtd) != 0:
        times_dt = [clk_start_dt + dt.timedelta(seconds = sec) for sec in seconds_dcmtd]
        results = np.column_stack((times_dt, pressure_dcmtd, temperature_resmpld))
    else:
        times_dt = [clk_start_dt + dt.timedelta(seconds = sec) for sec in seconds_p]
        results = np.column_stack((times_dt, pressure, temperature_resmpld))
    np.savetxt(out_filename, results, fmt = '%s,%f,%f')
    
    # Set min-max values for plot Y-axes
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
    p_min = p_min - p_min%intvl - intvl
    p_max = p_max - p_max%intvl + 2*intvl
    
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
    
    fig, ax1 = plt.subplots()
    color = 'red'
    plt.ylim(p_min,p_max)
    ax1.plot(seconds_p, pressure, color='red')
    if len(pressure_dcmtd) != 0:
        ax1.plot(seconds_dcmtd, pressure_dcmtd, color='yellow',marker='')
    ax1.set_xlabel('seconds')
    ax1.set_ylabel('Pressure (Pa)', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(axis='x')

    ax2 = ax1.twinx()
    color = 'blue'
    plt.ylim(t_min,t_max)
    ax2.set_ylabel('Temperature (°C)', color=color)
    ax2.plot(seconds_t, temperature, color=color)
    if len(pressure_dcmtd) != 0:
        ax2.plot(seconds_dcmtd, temperature_resmpld, color='green',marker='')
    ax2.tick_params(axis='y', labelcolor=color)
    fig.suptitle(f'Start: {wndw_start_dt}\nEnd:   {wndw_end_dt}')
    fig.tight_layout()

    plt.show()
    
    return

###########################################################################
def extractrecords(apg_filename, start_byte, nrecs_want, rec_len, rec_fmt):
    '''
    Extracts binary records from a raw APG data logger file.
    '''
    with open(apg_filename, 'rb') as apgfile:
        apgfile.seek(start_byte, os.SEEK_CUR)
        records = []
        for i in range(0,nrecs_want):
            record_bin = apgfile.read(rec_len)
            
            
            # Print record as a string of Hex and/or Bin values.
            # Uncomment below for troubleshooting.
            '''
            hex_str = ''
            bin_str = ''
            for ch in record_bin:
                #hex_str += hex(ch)+' '
                bin_str += f'{ch:08b}'
            #print( hex_str )
            cum_rec_fmt = np.cumsum(list(map(abs,rec_fmt)))
            new_bin_str = ''
            for i in range( len(bin_str) ):
                if i in cum_rec_fmt:
                    new_bin_str = new_bin_str + ' '
                new_bin_str = new_bin_str + bin_str[i]
            print( new_bin_str )
            '''
                        
            # Split record into array of ints defined as groups of bits by rec_fmt
            record_int = int.from_bytes( record_bin, byteorder='big', signed=False )
            record = []
            for signed_bit_len in reversed(rec_fmt):
                bit_len = int(abs(signed_bit_len))
                field = record_int & (2**bit_len-1) # Read right most bit_len bits
                # Check for sign bit and convert field as a 2s-compliment negative.
                if ((signed_bit_len < 0) and (field & (2**(bit_len-1)))):
                    full_bit_len = bit_len + (bit_len % 8)
                    field = field | ((2**full_bit_len-1) ^ (2**bit_len-1))
                    field = field.to_bytes(full_bit_len//8, byteorder='big', signed=False)
                    field = int.from_bytes(field, byteorder='big', signed=True)
                record.append(field)
                record_int = record_int >> (bit_len) # Shift to right bit_len bits
                            
            records.append(record)
            # Print a "." every 1000 records to indicate script is running.
            if i%1000 == 0:
                print('.',end='')
        print(f'\n')
    return(np.array(records))


#####################################################################
def clockdrift(clk_start_dt, gpssync_dt, sync_tick_count, tick_bits):
    '''
    Calculates the clock drift using one of two methods.
    The expected number of tick counts between clk_start_dt and gpssync_dt
    will be calculated. The size of the tick record in tick_bits is used to 
    determine when the tick count 'rolls over' to zero.
    If sync_tick_count is provided, the difference between this and the expected 
    value gives the drift in ticks (milliseconds).
    If sync_tick_count is not present, it is assumed that a fixed frequency has
    been injected into the raw data precisely at gpssync_dt. The precise tick
    count when this frequency starts is detected in the data and this value 
    is used in place of sync_tick_count.
    '''
    if (sync_tick_count == None):
        sys.exit('The code for time sync by frequency injection has not been written yet.')
    
    millisecs_logged = ((gpssync_dt - clk_start_dt).total_seconds()) * 1000
    # Number of ticks until rollover and restart at zero.
    tick_rollover = 2**tick_bits # 1 tick = 1 millisecond
    # APG logger clock offset relative to GPS time (positive = APG > GPS)
    clk_drift_at_end = int(sync_tick_count - (millisecs_logged % tick_rollover))

    return(clk_drift_at_end)
    

##########################
if __name__ == '__main__':
    main()