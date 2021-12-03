#!/usr/bin/env python3

import datetime as dt

def main():
    # Length of tic count field in bits.
    tic_bits = 34

    # UTIG2 - TX17-1
    clk_start = '2017-177 01:04:00' # 'YYYY-DOY hh:mm:ss'
    clk_end = '2018-281 05:27:30' # 'YYYY-DOY hh:mm:ss'
    tick_count_end = 0x017037C909
    clk_offset_at_end = clock_os(clk_start, clk_end, tick_count_end, tic_bits)
    print(f'Clock offset in milliseconds = {clk_offset_at_end}')

    # UTIG4 - TX17-2
    clk_start = '2017-176 22:34:00' # 'YYYY-DOY hh:mm:ss'
    clk_end = '2018-281 03:25:30' # 'YYYY-DOY hh:mm:ss'
    tick_count_end = 0x0170516C1D
    clk_offset_at_end = clock_os(clk_start, clk_end, tick_count_end, tic_bits)
    print(f'Clock offset in milliseconds = {clk_offset_at_end}')

    # UTIG3 - TX17-3
    clk_start = '2017-176 21:17:00' # 'YYYY-DOY hh:mm:ss'
    clk_end = '2018-281 06:18:30' # 'YYYY-DOY hh:mm:ss'
    tick_count_end = 0x0171364DB8
    clk_offset_at_end = clock_os(clk_start, clk_end, tick_count_end, tic_bits)
    print(f'Clock offset in milliseconds = {clk_offset_at_end}')

    # UTIG1 - TX17-4
    clk_start = '2017-177 00:36:30' # 'YYYY-DOY hh:mm:ss'
    clk_end = '2018-284 02:17:00' # 'YYYY-DOY hh:mm:ss'
    tick_count_end = 0x017f15a1ab
    clk_offset_at_end = clock_os(clk_start, clk_end, tick_count_end, tic_bits)
    print(f'Clock offset in milliseconds = {clk_offset_at_end}')


def clock_os(clk_start, clk_end, tick_count_end, tic_bits):
    clk_start_dt = dt.datetime.strptime(clk_start, '%Y-%j %H:%M:%S')
    clk_end_dt = dt.datetime.strptime(clk_end, '%Y-%j %H:%M:%S')
    seconds_logged = int((clk_end_dt - clk_start_dt).total_seconds())
    # Number of ticks untill rollover and restart at zero.
    tick_rollover = 2**tic_bits # 1 tick = 1 millisecond
    millisecs_logged = (seconds_logged * 1000)

    # APG logger clock offset relative to GPS time (positive = APG > GPS)
    clk_offset_at_end = tick_count_end - (millisecs_logged % tick_rollover)
    return(clk_offset_at_end)
    
if __name__ == '__main__':
    main()