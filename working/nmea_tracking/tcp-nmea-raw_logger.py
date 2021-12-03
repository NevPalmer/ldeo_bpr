#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 15 Oct 2021

@author: nev
"""
import sys
import socket
import io
import re
import datetime as dt


def Main():
    #host = '192.168.37.189'
    #port = 50000
    host = 'das.tan.niwa.co.nz'
    gps_port = 3000
    #gyro_port = 3006
    
    gps_sckt = getTCPsocket(host, gps_port)
    #gyro_sckt = getTCPsocket(host, gyro_port)

    while True:
        print(f'  ==>> {dt.datetime.now()}')
        streamNMEA(gps_sckt)
        #streamNMEA(gyro_sckt)

    gps_sckt.close()
    gyro_sckt.close()

def getTCPsocket(host, port):
    #try:
    nmea_sckt = socket.socket()
    nmea_sckt.connect((host, port))
    return(nmea_sckt)

def streamNMEA(sckt):
    data = sckt.recv(2048).decode('utf-8')
    buf = io.StringIO(data)
    line = buf.readline().strip()
    while line:
        print(line)
        line = buf.readline().strip()



if __name__ == '__main__':
    Main()
