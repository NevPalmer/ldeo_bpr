#! /bin/sh -
'eval' 'exec python3 -tt -- "$0" "$@"'

# Version 20181126
# by Neville Palmer, GNS Science
# 2018/11/17 Start dev
# 2018/11/23 SeaScan board format added.
# 2018/11/24 Smoothing and decimation added.
# 2018/11/25 Added ini files for reading in transducer & logger parameters.
# 2018/11/26 Resolved issue with iir decimation (ensure q<13 and n=5).
#               Iterate for desired decimation.

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
wndw_end =   '2017-07-01 12:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
'''
# Site TX17-4 used instrument UTIG-1,  136309
paros_sn = '136309'
logger_type = 'CSAC2013'
apg_filename = 'd:/UTIG_BPR_download/BPR_UTIG-1_TX17-4.apg';
clk_start =  '2017-06-26 00:36:30' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2017-07-01 00:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2017-07-01 12:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
'''
paros_sn = '140340'
logger_type = 'Seascan2018'
apg_filename = 'd:/UTIG_BPR_download/Test_paros140340_20181121.apg';
clk_start =  '2018-11-21 01:46:30' # 'YYYY-MM-DD hh:mm:ss'
#Times to extract data for.
wndw_start = '2018-11-21 01:46:30.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2018-11-21 04:20:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_start = '2018-11-21 02:30:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
#wndw_end =   '2018-11-21 02:30:01.0' # 'YYYY-MM-DD hh:mm:ss.ss'
'''

def main():
    # Pressure conversion factor from PSIA to Pascal.
    press_conv_fctr = 6.894757293168e3
    # Temperature smoothing factor (must be an odd integer)
    # 10001 gives strong smoothing
    tmptr_smth_fctr = 1001
    
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
            
            '''
            # Print record as string of Hex values
            str = ''
            for ch in record_bin:
                str += hex(ch)+" "
            print( str )
            '''
                        
            record_int = int.from_bytes( record_bin, byteorder='big', signed=False )
            # Split record into array of ints defined as groups of bits by rec_fmt
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
        millisec_t = ((ticks-64) / 1 )
        x = np.linspace(0,nrecs_want,nrecs_want+1)
        xinterp = np.linspace(0,nrecs_want-1/smpls_per_rec,nrecs_want*smpls_per_rec)
        ticks = np.append(ticks, (2*ticks[-1] - ticks[-2])) #extrapolate final record
        ticks = np.interp(xinterp, x, ticks)
        ticks = np.array([int(tick) for tick in ticks])
        millisec_p = ((ticks-64) / 1 )
        
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
        
        P_offset = np.mean(pressure)
        pressure_dec10 = sig.decimate(pressure,10,n=5,ftype='iir')
        millisec_dec10 = millisec_p[::10]
        pressure_dec100 = sig.decimate(pressure_dec10,10,n=5,ftype='iir')
        millisec_dec100 = millisec_p[::100]
        pressure_dec1000 = sig.decimate(pressure_dec100,10,n=5,ftype='iir')
        millisec_dec1000 = millisec_p[::1000]
        pressure_dec6000 = sig.decimate(pressure_dec1000,6,n=5,ftype='iir')
        millisec_dec6000 = millisec_p[::6000]

        '''
        pressure_B = pressure - P_offset
        pressure_B = sig.decimate(pressure_B,10,n=10,ftype='fir')
        pressure_B = sig.decimate(pressure_B,10,n=10,ftype='fir')
        pressure_B = sig.decimate(pressure_B,10,n=10,ftype='fir')
        pressure_B = sig.decimate(pressure_B,6,n=10,ftype='fir')
        pressure_B = pressure_B + P_offset
        millisec_pB = millisec_p[::6000]
        '''
        
#        print (ticks)
#        print(temperature)
#        print(pressure)
#        [print(press,end=', ') for press in pressure]
        
        # Set min-max values for plot Y-axes
        p_min = np.min(pressure)
        p_max = np.max(pressure)
        p_range = p_max-p_min
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
        
#        p_min = (np.mean(pressure)//1000)*1000 - 10000
#        t_min = (np.mean(temperature)//0.1)*0.1 - 0.5
        
        fig, ax1 = plt.subplots()
        color = 'red'
        plt.ylim(p_min,p_max)
        ax1.plot(millisec_p, pressure, color='red')
        ax1.plot(millisec_dec10, pressure_dec10, color='yellow')
        ax1.plot(millisec_dec100, pressure_dec100, color='green')
        ax1.plot(millisec_dec1000, pressure_dec1000, color='magenta')
        ax1.plot(millisec_dec6000, pressure_dec6000, color='cyan')
        ax1.set_xlabel('milliseconds')
        ax1.set_ylabel('Pressure (Pa)', color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(axis='x')

        ax2 = ax1.twinx()
        color = 'blue'
        plt.ylim(t_min,t_max)
        ax2.set_ylabel('Temperature (°C)', color=color)
        ax2.plot(millisec_t, temperature, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        fig.suptitle(f'Start: {wndw_start}\nEnd:   {wndw_end}')
        fig.tight_layout()

        plt.show()
    
    return

if __name__ == '__main__':
    main()