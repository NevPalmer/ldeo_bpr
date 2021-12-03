#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 16:51:37 2018

@author: nev
"""
import sys
import socket
import io
import re


def Main():
    #host = '192.168.37.189'
    #port = 50000
    host = 'das.tan.niwa.co.nz'
    port = 3000

    try:
        nmea_sckt = socket.socket()
        nmea_sckt.connect((host, port))
    except ConnectionRefusedError as error:
        print(f'Connection refused while attempting to connect to NMEA '
              f'server {host} on port {port}.')
        sys.exit()
    except OSError as error:
        print(error)
        print(f'while attempting to connect to NMEA server {host}.')
        sys.exit()

    for n in range(4):
        print(f'Iteration: {n}')
        data = nmea_sckt.recv(2048).decode('utf-8')
        buf = io.StringIO(data)
        line = buf.readline().strip()
        while line:
            nmealine = re.match(r'\$(.*)\*', line)[1]
            nmeaItem = (nmealine.split(','))

            if nmeaItem[0][1:] == 'PGGA':
                print(line)
                print(nmeaItem)
                print('Type: {}'.format(nmeaItem[0][1:]))
                print('UTC time: {}:{}:{}'.format(
                    nmeaItem[1][0:2], nmeaItem[1][2:4], nmeaItem[1][4:9]))
                print('Lat:  {}\u00b0{}\'{}'.format(
                    nmeaItem[2][0:2], nmeaItem[2][2:11], nmeaItem[3]))
                print('Lon: {}\u00b0{}\'{}'.format(
                    nmeaItem[4][0:3], nmeaItem[4][3:12], nmeaItem[5]))
                if nmeaItem[6] != '':
                    print('Fix quality: {}'.format(fixQlty(int(nmeaItem[6]))))
                    if nmeaItem[6] != '0':
                        print('No Sats: {:02d}'.format(int(nmeaItem[7])))
                        try:
                            print('HDOP: {:.1f}'.format(float(nmeaItem[8])))
                        except:
                            print('HDOP: not reported')
                        print('Ht AMSL: {:.1f} {}'.format(
                            float(nmeaItem[9]), nmeaItem[10].lower()))
                        print('Geoid sep: {:.1f} {}'.format(
                            float(nmeaItem[11]), nmeaItem[12].lower()))
                else:
                    print('Fix quality: not reported')
                print('Checksum: {} {}'.format(
                    hex(int(line[-2:], 16)),
                    chksum_nmea(line)
                ))

            elif nmeaItem[0][1:] == 'PVTG':
                print(line)
                print(nmeaItem)
                print('Type: {}'.format(nmeaItem[0][1:]))
                if nmeaItem[1] != '':
                    print(f'COG: {float(nmeaItem[1]):.1f}\u00b0 True')
                    print(f'SOG: {float(nmeaItem[5]):.1f} knots')
                    print(f'SOG: {float(nmeaItem[7]):.1f} km/h')
                else:
                    print('No Track Made Good reported.')
                print('Checksum: {} {}'.format(
                    hex(int(line[-2:], 16)),
                    chksum_nmea(line)
                ))
            line = buf.readline().strip()

        print("="*100)

    nmea_sckt.close()


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


if __name__ == '__main__':
    Main()
