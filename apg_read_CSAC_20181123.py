#! /bin/sh -
'eval' 'exec python3 -tt -- "$0" "$@"'

# Version 20181123 CSAC
# by Neville Palmer, GNS Science
# 2018/11/17 Start dev

# APG data format is:
#   File header, 512 bytes (16 x 32 byte records).
#   Sequential 32-byte (256-bit) records, each containing one tic count,
#   one temperature sample and eight pressure samples.
#   1.    tic count    34 bits unsigned int
#   2.    Temp         40 bits unsigned int
#   3.    unused        2 bits
#   4-10. P8-P2        20 bits signed int each (7 difference values)
#   11.   Pcore        40 bits unsigned int
# Note: The eight pressure values are unpacked as:
#   [Pcore Pcore+P2 Pcore+P2+P3 Pcore+P2+P3+P4... Pcore+P2+..+P8]
#   P2-P8 are signed (2s compliment) integers.

import os
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

###change below for each BPR, as well as coefficients
'''
# TX17-1 used UTIG-2, APG serial number 136863
apg_filename = 'd:/UTIG_BPR_download/BPR_UTIG-2_TX17-1.apg';
U0 = 5.818855
Y = [-3986.358, -10512.00, 0]
C = [-31892.15,142.9143,87999.78]
D = [0.040572, 0]
T = [30.18210,1.340158, 54.26334, 139.7080,0]
temp_clock=50  # Mhz
clk_start =  '2017-06-26 01:04:00' # 'YYYY-MM-DD hh:mm:ss'

# TX17-2 used UTIG-4, APG serial number 236308
apg_filename = 'd:/UTIG_BPR_download/BPR_UTIG-4_TX17-2.apg';
U0 = 5.826349
Y = [-3960.002, -11216.04, 0]
C = [-28133.56,649.2896,87462.27]
D = [0.038439, 0]
T = [30.31328,1.82197, 58.6033, 139.7210,0]
temp_clock=50  # Mhz
clk_start =  '2017-06-25 22:34:00' # 'YYYY-MM-DD hh:mm:ss'

# TX17-3 used UTIG-3, APG serial number 136305
apg_filename = 'd:/UTIG_BPR_download/BPR_UTIG-3_TX17-3.apg';
U0 = 5.774811
Y = [-3993.698, -11177.66, 0]
C = [-24594.27,-99.44320,79270.41]
D = [0.039887, 0]
T = [30.06436,0.852899, 53.6318, 93.96366,0]
temp_clock=50  # Mhz
clk_start =  '2017-06-25 21:17:00' # 'YYYY-MM-DD hh:mm:ss'
'''
# TX17-4 used UTIG-1, APG serial number 136309
apg_filename = 'd:/UTIG_BPR_download/BPR_UTIG-1_TX17-4.apg';
U0 = 5.848144
Y = [-3962.094, -11170.40, 0]
C = [-29930.79,664.1395,95597.41]
D = [0.039290, 0]
T = [30.01267,1.200911, 55.91407, 144.6710,0]
temp_clock=50  # Mhz
clk_start =  '2017-06-26 00:36:30' # 'YYYY-MM-DD hh:mm:ss'

#Times to extract data for.
wndw_start = '2017-07-01 00:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'
wndw_end =   '2017-07-11 00:00:00.0' # 'YYYY-MM-DD hh:mm:ss.ss'

#CSAC Board
epoch_logged = 0.01 # seconds per sample as logged
epoch_decmtd = 60 # required seconds per sample after decimating
rec_len = 32 # record length in bytes
head_size = 16 # header size in records
head_len = head_size * rec_len
smpls_per_rec = 8 # number of pressure samples per record
# Format of records as a tuple of integer numbers indicating bits per field.
# A negative indicates a signed (2s compliment) value.
# Values are: 34 Tick count, 40 Temperature, 2 free, 7 x -20 P2-P8, 40 P Core
rec_fmt = ( 34, 40, 2, -20, -20, -20, -20, -20, -20, -20, 40 )

def main():
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
    rec_start = int(wndw_start_secs / epoch_logged / smpls_per_rec)
    print(f'rec_start = {rec_start}')
    nrecs_want = int(wndw_len_secs / epoch_logged / smpls_per_rec)+1
    
    fsize = os.path.getsize(apg_filename) # file size in bytes
    print(f'Filesize = {fsize} bytes')
    nrecs = int((fsize - head_len) / rec_len)
    print(f'Number of records = {nrecs:d}')
    file_duration_days = nrecs*smpls_per_rec*epoch_logged/3600/24
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
                bit_len = abs(signed_bit_len)
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
            if i%1000 == 0:
                print('.',end='')
        records = np.array(records)
        print(f'\n')
        
        press_raw = records[:,0:8]
#        print(press_raw)
        press_raw = np.cumsum(press_raw,axis=1)
        press_raw = press_raw.reshape((nrecs_want*smpls_per_rec))
        
        temp_raw = (records[:,9:10])
#        print(temp_raw)
        
        ticks = (records[:,10])
        millisec_t = ((ticks-64) / 1 )
        x = np.linspace(0,nrecs_want,nrecs_want+1)
        xinterp = np.linspace(0,nrecs_want-1/smpls_per_rec,nrecs_want*smpls_per_rec)
        ticks = np.append(ticks, (2*ticks[-1] - ticks[-2])) #extrapolate final record
        ticks = np.interp(xinterp, x, ticks)
        ticks = np.array([int(tick) for tick in ticks])
        millisec_p = ((ticks-64) / 1 )
        
        # Temperature period (usec)
        TP = (temp_raw/(2**28)+1)/(temp_clock)
        # Pressure period (usec)
        PP = (press_raw/(2**25)+1)/(temp_clock)
        
        # Create temperature
        U = TP-U0;
        temperature = np.polyval([Y[2], Y[1], Y[0], 0], U)
        temperature = temperature.reshape((nrecs_want))
        U_expnd = U * np.ones((1,smpls_per_rec),dtype=int)
        U_expnd = U_expnd.reshape((nrecs_want*smpls_per_rec))

        # Create pressure
        Cv = np.polyval([C[2], C[1], C[0]], U_expnd)
        Dv = np.polyval([D[1], D[0]], U_expnd)
        T0 = np.polyval([T[4], T[3], T[2], T[1], T[0]], U_expnd)
        
        facts = 1 - (T0**2)/(PP**2) 
        P_psia = Cv * facts *(1-Dv*facts)
        P_pa = P_psia * 6.894757293168e3
        
#        print (ticks)
#        print(temperature)
#        print(P_pa)
#        [print(press,end=', ') for press in P_pa]
        
        # Set minimum values for plot Y-axes
        p_min = (np.mean(P_pa)//1000)*1000 - 10000
        t_min = (np.mean(temperature)//0.1)*0.1 - 0.5
        
        fig, ax1 = plt.subplots()
        color = 'tab:red'
        plt.ylim(p_min,p_min+20000)
        ax1.set_xlabel('milliseconds')
        ax1.set_ylabel('Pressure (Pa)', color=color)
        ax1.plot(millisec_p, P_pa, color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax2 = ax1.twinx()
        color = 'tab:blue'
        plt.ylim(t_min,t_min+1.0)
        ax2.set_ylabel('Temperature (°C)', color=color)
        ax2.plot(millisec_t, temperature, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax1.grid(axis='x')
        fig.suptitle(f'Start: {wndw_start}\nEnd:   {wndw_end}')
        fig.tight_layout()
        plt.show()

    
    return


if __name__ == '__main__':
    main()