#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 19:49:21 2024

@author: zhao
"""

import numpy as np
from pyginvc.libs import geotools as gt
import argparse

def gen_fautgeom(origin, leng, width, dep, strike, dip, ss, ds, op):
    '''
    convert to faultgeom format
    
    origin : [lat, lon] of the center of the fault patch
    leng   : length of the fault patch
    width  : width of the fault patch
    dep    : depth of the center of the fault patch
    srike  : strike angle in degree
    dip    : dip angle in degree
    ss     : strike slip in millimeter
    ds     : dip slip in millimeter
    op     : tensile slip in millimeter
    '''
    
    rdip    = np.deg2rad(dip)
    rstrk   = np.deg2rad(strike)
    pt1N    = -leng/2
    pt1E    = width/2
    pt2N    =  leng/2
    pt2E    = width/2
    pt1E    = pt1E*np.cos(rdip)
    pt2E    = pt2E*np.cos(rdip)
    pt1Z    = dep-width/2*np.sin(rdip)
    pt2Z    = dep-width/2*np.sin(rdip)
    depth   = pt1Z

    pt1New  = complex(pt1N, pt1E)*np.exp(complex(0, rstrk))
    pt2New  = complex(pt2N, pt2E)*np.exp(complex(0, rstrk))

    slatlon = gt.utm2llh([pt1New.real, pt1New.imag], origin)
    elatlon = gt.utm2llh([pt2New.real, pt2New.imag], origin)
    geom    = np.array([width, depth, dip+180, slatlon[0,0], slatlon[0,1], elatlon[0,0], elatlon[0,1], ss, ds, op])
    return geom

def read_geom_param(fspfile):
    with open(fspfile, 'r') as fid:
        lines = fid.readlines()
    for line in lines:
        if line.find('Dx') != -1:
            length = float(line.split()[5])
            width  = float(line.split()[9])
        if line.find('STRK')!=-1:
            strike = float(line.split()[5])
            dip    = float(line.split()[8])
    return length, width, strike, dip

def main():
    parser = argparse.ArgumentParser(description="Convert USGS coseismic slip model in fsp file to faultgeom file.")
    parser.add_argument('--fsp_file', type=str, required=True)
    args      = parser.parse_args()

    fspfile   = args.fsp_file
    dat       = np.genfromtxt(fspfile, comments="%")
    faultgeom = np.zeros((len(dat),10))
    length, width, strike, dip = read_geom_param(fspfile)
    for i in range(len(dat)):
        origin = dat[i,0:2]
        depth  = dat[i,4]
        slip   = dat[i,5]
        rake   = dat[i,6]
        radrak = np.deg2rad(rake)
        ss     = slip*np.cos(radrak)*1e3
        ds     = slip*np.sin(radrak)*1e3
        op     = 0.0
    
        geom   = gen_fautgeom(origin, length, width, depth, strike, dip, ss, ds, op)
        faultgeom[i] = geom
    
    header  = " width_km   depth_km    dip_deg   lat0_deg   lon0_deg   lat1_deg   lon1_deg  slp_strk_mm  slp_ddip_mm  slp_tens_mm"
    fmt     = 7*'%11.3f'+3*'%13.3f'
    np.savetxt('output.faultgeom', faultgeom, fmt=fmt, header=header)


if __name__ == '__main__':
    main()
