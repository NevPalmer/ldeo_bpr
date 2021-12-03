#! /bin/sh -
'eval' 'exec python3 -tt -- "$0" "$@"'

# Version 20181128
# by Neville Palmer, GNS Science
# 2018/11/17 Start dev
# 2018/11/23 SeaScan board format added.
# 2018/11/24 Smoothing and decimation added.
# 2018/11/25 Added ini files for reading in transducer & logger parameters.
# 2018/11/26 Resolved issue with iir decimation (ensure q<13 and n=5).
#               Iterate for desired decimation.
# 2018/11/28 Add iteration routine for decimation.

import os
import sys
from configparser import ConfigParser
from numexpr import evaluate
import datetime as dt
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as sig

paros_ini = 'ParosAPG.ini'
apglogger_ini = 'APGlogger.ini'
'''
###change below for BPR used,
# Site TX17-2 used instrument UTIG-4
paros_sn = '136308'
logger_type = 'CSAC2013'
apg_filename = 'd:/UTIG_BPR_download/BPR_UTIG-4_TX17-2.apg';
clk_start =  '2017-06-25 22:34:00' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2017-07-01 00:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2017-07-01 06:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'

# Site TX17-4 used instrument UTIG-1,  136309
paros_sn = '136309'
logger_type = 'CSAC2013'
apg_filename = 'd:/UTIG_BPR_download/BPR_UTIG-1_TX17-4.apg';
clk_start =  '2017-06-26 00:36:30' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2017-07-01 00:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2017-07-01 12:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'


# SeaScan APG logger test 20181121
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181121.apg';
clk_start =  '2018-11-21 02:20:30' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
#wndw_start = '2018-11-21 02:20:30.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-21 04:30:50.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#This is the moment when 50kHz frequency is inserted for stop timing.
#wndw_start = '2018-11-21 04:30:50.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-21 04:31:10.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#Section of 50kHz frequency inserted for stop timing.
#wndw_start = '2018-11-21 04:30:52.6' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-21 04:31:20.0' # 'YYYY-MM-DD hh:mm:ss.ss'

# SeaScan APG logger test 20181123
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181123.apg';
clk_start =  '2018-11-23 02:50:00' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
#wndw_start = '2018-11-27 01:53:30.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-27 01:54:30.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_start = '2018-11-24 00:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-24 00:10:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'

# SeaScan APG logger test 20181127
paros_sn = '140338'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140338_20181127.apg';
clk_start =  '2018-11-27 05:25:00' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2018-11-27 06:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-27 06:10:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_start = '2018-11-27 19:39:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-27 19:45:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'

# SeaScan APG logger test 20181127
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181127.apg';
clk_start =  '2018-11-27 20:51:00' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2018-11-27 21:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-27 21:30:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_start = '2018-11-27 21:40:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-27 21:44:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'

# SeaScan APG logger test 20181127-2
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181127-2.apg';
clk_start =  '2018-11-27 23:52:00' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2018-11-27 23:57:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-28 00:03:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_start = '2018-11-28 01:42:05.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-28 01:46:4.0' # 'YYYY-MM-DD hh:mm:ss.ss'

# SeaScan APG logger test 20181127-3
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181127-3.apg';
clk_start =  '2018-11-28 01:40:00' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
#wndw_start = '2018-11-28 01:40:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-28 01:40:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_start = '2018-11-28 01:42:01.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-28 01:46:30.0' # 'YYYY-MM-DD hh:mm:ss.ss'

# SeaScan APG logger test 20181127-4
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181127-4.apg';
clk_start =  '2018-11-28 02:24:30' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2018-11-28 02:24:30.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-28 02:38:30.0' # 'YYYY-MM-DD hh:mm:ss.ss'
'''
# SeaScan APG logger test 20181127-5
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181127-5.apg';
clk_start =  '2018-11-28 03:05:00' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2018-11-28 03:09:59.5' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-28 03:10:00.5' # 'YYYY-MM-DD hh:mm:ss.ss'
'''
# SeaScan APG logger test 20181128-4
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181128-4.apg';
clk_start =  '2018-11-28 23:35:30' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2018-11-28 23:37:59.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-28 23:38:01.0' # 'YYYY-MM-DD hh:mm:ss.ss'
'''
def main():
    # Pressure conversion factor from PSIA to Pascal.
    press_conv_fctr = 6.894757293168e3
    # Temperature smoothing factor (must be an odd integer)
    # 10001 gives strong smoothing
    tmptr_smth_fctr = 1
    # Required sample interval in seconds for pressure decimation.
    # Zero for no decimation.
    decmt_intvl = 0
    
    # Read Paros transducer coefficients into a dict of lists from ini file.
    paros_coefs = ('U', 'Y', 'C', 'D', 'T')
    paros_cfg = ConfigParser()
    paros_cfg.read(paros_ini)
    paros = {}
    for coef in paros_coefs:
        paros[coef] = (paros_cfg.get(paros_sn,coef).split(','))
        paros[coef] = [float(x) for x in paros[coef][::-1]]
    paros['Y'].append(0.0)
    
    # Read APG logger configuration parameters from ini file.
    logger_cfg = ConfigParser()
    logger_cfg.read(apglogger_ini)
    clock_freq = logger_cfg.getint(logger_type,'clock_freq')
    sample_epoch = logger_cfg.getfloat(logger_type,'epoch')
    rec_len = logger_cfg.getint(logger_type,'rec_len')
    head_len = logger_cfg.getint(logger_type,'head_len')
    smpls_per_rec = logger_cfg.getint(logger_type,'smpls_per_rec')
    print(logger_cfg.get(logger_type,'TP_fctr'))
    TP_fctr = evaluate(logger_cfg.get(logger_type,'TP_fctr')).item()
    TP_cnst = logger_cfg.getfloat(logger_type,'TP_cnst')
    PP_fctr = evaluate(logger_cfg.get(logger_type,'PP_fctr')).item()
    PP_cnst = logger_cfg.getfloat(logger_type,'PP_cnst')
    
    rec_fmt = (logger_cfg.get(logger_type,'rec_fmt').split(','))
    rec_fmt = tuple([float(x) for x in rec_fmt])
    fmt_fields = {}
    fmt_fields['tic'] = logger_cfg.getint(logger_type,'tic_field')
    fmt_fields['tptr'] = logger_cfg.getint(logger_type,'temperature_field')
    fmt_fields['pcore'] = logger_cfg.getint(logger_type,'pcore_field')
    fmt_fields['pn'] = logger_cfg.getint(logger_type,'pn_field')
    if (abs(fmt_fields['pcore']-fmt_fields['pn']+1) != smpls_per_rec):
        sys.exit(f"The number of samples per record ({smpls_per_rec}) "
                 f"for logger {logger_type}, does "
                 f"not match the number of Pcore and Pn fields in the "
                 f"record format (Pcore through Pn inclusive = "
                 f"{fmt_fields['pcore']-fmt_fields['pn']+1}), "
                 f"as provided in file {apglogger_ini}.")
    rec_fmt_bits = int(sum(map(abs,rec_fmt)))
    if (rec_fmt_bits != (rec_len*8)):
        sys.exit(f"The total number of bits ({rec_fmt_bits}) given by the "
                 f"record format {rec_fmt} "
                 f"does not equal the record length ({rec_len} bytes x 8), "
                 f"as provided in the file {apglogger_ini}.")


    clk_start_dt = dt.datetime.strptime(clk_start, '%Y-%m-%d %H:%M:%S')
    wndw_start_dt = dt.datetime.strptime(wndw_start, '%Y-%m-%d %H:%M:%S.%f')
    wndw_end_dt = dt.datetime.strptime(wndw_end, '%Y-%m-%d %H:%M:%S.%f')
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
    
    fsize = os.path.getsize(apg_filename) # file size in bytes
    print(f'Filesize = {fsize} bytes')
    nrecs = int((fsize - head_len) / rec_len)
    print(f'Number of records = {nrecs:d}')
    file_duration_days = nrecs*smpls_per_rec*sample_epoch/3600/24
    print(f'File duration = {file_duration_days} days')
    clk_end_dt = clk_start_dt + dt.timedelta(days=file_duration_days)
    print(f'Time of last record = {clk_end_dt}')

    with open(apg_filename, 'rb') as apgfile:
        apgfile.seek(head_len + rec_start * rec_len, os.SEEK_CUR)
        records = []
        for i in range(0,nrecs_want):
            record_bin = apgfile.read(rec_len)
            
            
            # Print record as string of Hex and/or Bin values uncomment below for troubleshooting
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
        records = np.array(records)
        print(f'\n')
        
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
            TP_offset = np.mean(TP)
            TP = sig.savgol_filter(TP-TP_offset, tmptr_smth_fctr, 3, axis=0)
            TP = TP + TP_offset
        
        # Create temperature
        Uv = TP-paros['U'][0];
        temperature = np.polyval(paros['Y'], Uv)
        temperature = temperature.reshape((nrecs_want))
        Uv_expnd = Uv * np.ones((1,smpls_per_rec),dtype=int)
        Uv_expnd = Uv_expnd.reshape((nrecs_want*smpls_per_rec))

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
            decmt_intvl_pre = 1/sample_epoch
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
        
#        print (ticks)
#        print(temperature)
#        print(pressure)
#        [print(press,end=', ') for press in pressure]
        
        # Set min-max values for plot Y-axes
        p_min = np.min(pressure)
        p_max = np.max(pressure)
        p_range = p_max-p_min
        if p_range == 0:
            intvl = 10
        else:
            int_log = int(np.log10(p_range))
            if 10**int_log/p_range >= 0.5:
                intvl = 10**int_log/10
            elif 10**int_log/p_range >= 0.2:
                intvl = 10**int_log/5
            elif 10**int_log/p_range >= 0.1:
                intvl = 10**int_log/2
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
            if 10**int_log/t_range >= 0.5:
                intvl = 10**int_log/10
            elif 10**int_log/t_range >= 0.2:
                intvl = 10**int_log/5
            elif 10**int_log/t_range >= 0.1:
                intvl = 10**int_log/2
            else:
                intvl = 10**int_log
        t_min = t_min - t_min%intvl - intvl
        t_max = t_max - t_max%intvl + 2*intvl
        
        fig, ax1 = plt.subplots()
        color = 'red'
        plt.ylim(p_min,p_max)
        ax1.plot(seconds_p, pressure, color='red')
        if len(pressure_dcmtd) != 0:
            ax1.plot(seconds_dcmtd, pressure_dcmtd, color='yellow')
        ax1.set_xlabel('seconds')
        ax1.set_ylabel('Pressure (Pa)', color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(axis='x')

        ax2 = ax1.twinx()
        color = 'blue'
        plt.ylim(t_min,t_max)
        ax2.set_ylabel('Temperature (°C)', color=color)
        ax2.plot(seconds_t, temperature, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        fig.suptitle(f'Start: {wndw_start}\nEnd:   {wndw_end}')
        fig.tight_layout()

        plt.show()
    
    return

if __name__ == '__main__':
    main()