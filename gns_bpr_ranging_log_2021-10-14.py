#!/usr/bin/env python3

# Version 20200818
# by Neville Palmer, GNS Science
# 2020/08/18 Start development based on apg_read.py
# 2021/10/13 Add TCP-NMEA stream

import os
import sys
import argparse
import configparser
import csv
import serial
import time
import re
import datetime as dt
import numpy as np
import socket
import io

# Other constants.
sound_speed = 1500  # Speed of sound in water (typical 1450 to 1570 m/sec)
turn_around_time = 12.5 # Delay in ms for reply from BPR transducer.


def main():
    now = dt.datetime.now()
    timestamp = now.strftime('%Y-%m-%d_%H-%M')
    # Default values and choices for reading params from command line.
    out_filename = f'./bpr_rangelog_{timestamp}.csv'
    etech_ser_port = 'COM1'
    etech_ser_baud = 9600
    etech_ser_parity = 'N'
    etech_ser_stop = 1
    etech_ser_bytesize = 8
    etech_turnTime = 12.5
    etech_sndSpd = 1500
    nmea_host = 'das.tan.niwa.co.nz'
    nmea_port = 3000
    
    COLUMNS = ('utcTime', 'rangeTime', 'range', 'lat', 'latDec', 'lon', 'lonDec', 'qlty', 'noSats', 'hdop', 'htAmsl', 'htAmslUnit', 'geiodSep', 'geiodSepUnit', 'cog', 'sogKt', 'sogKm', 'turnTime', 'sndSpd', 'tx', 'rx')
    DISPLAY_COLS = ('utcTime', 'rangeTime', 'range', 'lat', 'lon', 'cog', 'sogKt')

    # Read in parameters from command line
    helpdesc = (
        'Communicate with Edgetech 8011M acoustic deck box via serial.\n'
        'Send and receive ranging commands and store ranges.')
    parser = argparse.ArgumentParser(description=helpdesc)
    parser.add_argument("--outfile",
                        help=f'Full path and name of output CSV file. Default: "{out_filename}"',
                        default=out_filename)
    parser.add_argument("--serport",
                        help=f'Serial port name. Default: "{etech_ser_port}"',
                        default=etech_ser_port)
    parser.add_argument("--serbaud",
                        help=f'Serial baud rate. Default: {etech_ser_baud}',
                        default=etech_ser_baud)
    parser.add_argument("--serparity",
                        help=f'Serial parity. Default: "{etech_ser_parity}"',
                        default=etech_ser_parity)
    parser.add_argument("--serstop",
                        help=f'Serial stop bit name. Default: {etech_ser_stop}',
                        default=etech_ser_stop)
    parser.add_argument("--serbytesize",
                        help=f'Serial byte size. Default: {etech_ser_bytesize}',
                        default=etech_ser_bytesize)
    parser.add_argument("--host",
                        help=f'Host address for TCP NMEA stream. Default: {nmea_host}',
                        default=nmea_host)
    parser.add_argument("--port",
                        help=f'Host port for TCP NMEA stream. Default: {nmea_port}',
                        default=nmea_port)
    args = parser.parse_args()

    # Translate argparse parameters.
    out_filename = os.path.normpath(args.outfile)
    etech_ser_port = args.serport.strip()
    etech_ser_baud = int(args.serbaud)
    etech_ser_parity = args.serparity.strip()
    etech_ser_stop = int(args.serstop)
    etech_ser_bytesize = int(args.serbytesize)
    nmea_host = args.host.strip()
    nmea_port = int(args.port)

    # Open TCP port to NMEA server.
    try:
        nmea_sckt = socket.socket()
        nmea_sckt.connect((nmea_host, nmea_port))
    except ConnectionRefusedError as error:
        print(f'Connection refused while attempting to connect to NMEA '
              f'server {nmea_host} on port {nmea_port}.')
        sys.exit()
    except OSError as error:
        print(error)
        print(f'while attempting to connect to NMEA server {nmea_host}.')
        sys.exit()

    # Open serial port to EdgeTech 8011M deckbox.
    with serial.Serial(port=etech_ser_port,
                       baudrate=etech_ser_baud,
                       parity=etech_ser_parity,
                       stopbits=etech_ser_stop,
                       bytesize=etech_ser_bytesize,
                       timeout=0.01) as ser:
        print("Connected to EgeTech deckbox: " + ser.portstr)
        
        print(','.join(DISPLAY_COLS))
        with open(out_filename, 'w') as csvfile:
            dw = csv.DictWriter(csvfile, delimiter=',', fieldnames=COLUMNS)
            dw.writeheader()

        while True:
            nmeaPGGA, nmeaPVTG = getNMEA(nmea_sckt)
            response, flag = getResponse(ser)
            response = response.decode('UTF-8').strip().split(' ')
            if (response[0] == 'RNG:'):
                print('==================================================================')
                #print(f'Response Flag: {flag}')
                #print(response)
                range = {}
                range['turnTime'] = etech_turnTime
                range['sndSpd'] = etech_sndSpd
                range['tx'] = float(response[3])
                range['rx'] = float(response[6])
                range['rangeTime'] = float(response[9])
                range['range'] = (range['rangeTime'] /2) * range['sndSpd']
                

                #print(nmeaPGGA)
                #print(nmeaPVTG)

                nmea = NMEAtoDict(nmeaPGGA, nmeaPVTG)
                range |= nmea
            
                values = []
                for key in DISPLAY_COLS:
                    values.append(range[key])
                print(str(values).strip('[]'))

                with open(out_filename, 'a+', newline='') as csvfile:
                    dw = csv.DictWriter(csvfile, delimiter=',', fieldnames=COLUMNS)
                    dw.writerow(range)

                #print(range)
            
    nmea_sckt.close()

    return


def getResponse(ser):
    response = []
    byte_2 = b''
    byte_1 = b''
    byte_0 = ser.read(1)
    while byte_0 != b'':  # next_byte will be '' after serial.timeout
        response.append(byte_0)
        # byte_2: '*' indicates success, '#' indicates error
        if (byte_2 + byte_1 + byte_0) in (b'*\r\n', b'#\r\n'):
            # print('Command completed.')
            flag = b''.join(response[-3:-2])
            response = b''.join(response[0:-5])
            return response, flag
        byte_2 = byte_1
        byte_1 = byte_0
        byte_0 = ser.read(1)
    flag = b'T'
    response = b''.join(response)

    return response, flag


def sendCommand(ser, command):
    ser.write(command)
    # print(f'Sent command: {command}')

    return


def getNMEA(nmea_sckt):
    def getLine(nmea_sckt):
        data = nmea_sckt.recv(2048).decode('utf-8')
        nmeaLine = io.StringIO(data).read().strip()
        if not chksum_nmea(nmeaLine):
            return(false)
        nmeaItem = re.match(r'\$(.*)\*', nmeaLine)[1].split(',')
        return(nmeaItem)

    nmeaPGGA = getLine(nmea_sckt)
    while nmeaPGGA[0][1:] != 'PGGA':
        nmeaPGGA = getLine(nmea_sckt)
        
    nmeaPVTG = getLine(nmea_sckt)
    while nmeaPVTG[0][1:] != 'PVTG':
        nmeaPVTG = getLine(nmea_sckt)
        
    return(nmeaPGGA, nmeaPVTG)
    

def NMEAtoDict(nmeaPGGA, nmeaPVTG):
    nmea = {}
    nmea['utcTime'] = f"{nmeaPGGA[1][0:2]}:{nmeaPGGA[1][2:4]}:{nmeaPGGA[1][4:9]}"
    nmea['lat'] = f"{nmeaPGGA[2][0:2]}\u00b0{nmeaPGGA[2][2:11]}\'{nmeaPGGA[3]}"
    nmea['latDec'] = int(nmeaPGGA[2][0:2]) + float(nmeaPGGA[2][2:11]) / 60
    nmea['lon'] = f"{nmeaPGGA[4][0:3]}\u00b0{nmeaPGGA[4][3:12]}\'{nmeaPGGA[3]}"
    nmea['lonDec'] = int(nmeaPGGA[4][0:3]) + float(nmeaPGGA[4][3:12]) / 60
    nmea['qlty'] = fixQlty(int(nmeaPGGA[6]))
    nmea['noSats'] = int(nmeaPGGA[7])
    nmea['hdop'] = float(nmeaPGGA[8])
    nmea['htAmsl'] = float(nmeaPGGA[9])
    nmea['htAmslUnit'] = nmeaPGGA[10].lower()
    nmea['geiodSep'] = float(nmeaPGGA[11])
    nmea['geiodSepUnit'] = nmeaPGGA[12].lower()
    
    nmea['cog'] = float(nmeaPVTG[1])
    nmea['sogKt'] = float(nmeaPVTG[5])
    nmea['sogKm'] = float(nmeaPVTG[7])
                
    return(nmea)


def fixQlty(n):
    return [
        'Invalid',
        'GPS fix - (Standard Positioning Service)',
        'DGPS fix',
        'PPS fix',
        'Real Time Kinematic',
        'Float RTK',
        'Estimated (dead reckoning)',
        'Manual input mode',
        'Simulation mode'
    ][n]


def chksum_nmea(sentence):
    sentenceMatch = re.match(r'\$(.*)\*(.{2})', sentence)
    chksumdata = sentenceMatch[1]
    cksum = int(sentenceMatch[2], 16)  # convert from hex string
    csum = 0

    # For each char in chksumdata, XOR against the previous XOR'd char.
    # The final XOR of the last char will be our checksum to verify
    # against the checksum we sliced off the NMEA sentence.
    for c in chksumdata:
        # XOR'ing value of csum against the next char in line
        # and storing the new XOR value in csum
        csum ^= ord(c)

    # Do we have a validated sentence?
    if csum == cksum:
        return True
    else:
        return False


##########################
if __name__ == '__main__':
    main()
