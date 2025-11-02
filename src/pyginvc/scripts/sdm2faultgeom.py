#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 19:49:21 2024

@author: zhao
"""

import numpy as np
import pyginvc.libs.geotools as gt
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
    pt2Z    = dep+width/2*np.sin(rdip)
    depth   = pt1Z

    pt1New  = complex(pt1N, pt1E)*np.exp(complex(0, rstrk))
    pt2New  = complex(pt2N, pt2E)*np.exp(complex(0, rstrk))

    slatlon = gt.utm2llh([pt1New.real, pt1New.imag], origin)
    elatlon = gt.utm2llh([pt2New.real, pt2New.imag], origin)
    geom    = np.array([width, depth, dip+180, slatlon[0,0], slatlon[0,1], elatlon[0,0], elatlon[0,1], ss, ds, op])
    return geom

def main():
    parser = argparse.ArgumentParser(description="Convert SDM file to faultgeom file.")
    parser.add_argument('--sdm_file', type=str, required=True)
    parser.add_argument('--ss_sign', type=float, default=1.0, required=False)
    parser.add_argument('--ds_sign', type=float, default=-1.0, required=False)
    args = parser.parse_args()

    sdmfile   = args.sdm_file
    dat       = np.genfromtxt(sdmfile, skip_header=0)
    faultgeom = np.zeros((len(dat),10))
    ss_sign   = args.ss_sign
    ds_sign   = args.ds_sign
    for i in range(len(dat)):
        length = dat[i,6]
        width  = dat[i,5]
        strike = dat[i,10]
        dip    = dat[i,11]
        depth  = dat[i,2]
        origin = dat[i,0:2]
        ss     =-dat[i,7]*1e3*ss_sign
        ds     = dat[i,8]*1e3*ds_sign
        op     = 0.0
    
        geom   = gen_fautgeom(origin, length, width, depth, strike, dip, ss, ds, op)
        faultgeom[i] = geom

    if len(faultgeom)>2 and abs(faultgeom[0,1]-faultgeom[1,1])>1.0:
        idx    = np.arange(0,len(faultgeom))
        nseg   = len(np.where(abs(faultgeom[:,1]-faultgeom[0,1])<0.1)[0])
        ndip   = len(faultgeom)//nseg
        print(nseg, ndip, len(idx))
        idx    = idx.reshape((nseg, ndip)).T.flatten()
        faultgeom = faultgeom[idx]
    
    header  = " width_km   depth_km    dip_deg   lat0_deg   lon0_deg   lat1_deg   lon1_deg  slp_strk_mm  slp_ddip_mm  slp_tens_mm"
    fmt     = 7*'%11.3f'+3*'%13.3f'
    np.savetxt('output.faultgeom', faultgeom, fmt=fmt, header=header)


if __name__ == '__main__':
    main()
