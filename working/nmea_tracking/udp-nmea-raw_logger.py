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
    ip_addr = '0.0.0.0'
    gps = 6000
    dgps = 6061
    gyro = 6006
    windsonic_port = 6088
    windsonic_stbd = 6089
    
    gps_sckt = getUDPsocket(ip_addr, gps)
    dgps_sckt = getUDPsocket(ip_addr, dgps)
    gyro_sckt = getUDPsocket(ip_addr, gyro)
    wind_port_sckt = getUDPsocket(ip_addr, windsonic_port)
    wind_stbd_sckt = getUDPsocket(ip_addr, windsonic_stbd)

    while True:
        streamNMEA(gps_sckt, gps)
        #streamNMEA(dgps_sckt, dgps)
        #streamNMEA(gyro_sckt, gyro)
        #streamNMEA(wind_port_sckt, windsonic_port)
        #streamNMEA(wind_stbd_sckt, windsonic_stbd)

    gps_sckt.close()
    dgps_sckt.close()
    gyro_sckt.close()
    wind_port_sckt.close()
    wind_stbd_sckt.close()

def getUDPsocket(ip_addr, port):
    nmea_sckt = socket.socket(socket.AF_INET, # Internet
                              socket.SOCK_DGRAM) # UDP
    nmea_sckt.bind((ip_addr, port))
    nmea_sckt.settimeout(0.001)
    return(nmea_sckt)

def streamNMEA(sckt, port):
    try:
        data, addr = sckt.recvfrom(2048)
        data = data.decode('utf-8')
        buf = io.StringIO(data)
        line = buf.readline().strip()
        while line:
            print(f'  ==>> {dt.datetime.now()}')
            print(f'{port}: {line}')
            line = buf.readline().strip()
    except socket.timeout:
        pass


if __name__ == '__main__':
    Main()
