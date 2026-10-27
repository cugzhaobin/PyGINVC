#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 00:26:56 2020

@author: zhao
"""

import h5py
import numpy as np
import geotools as gt
import math, sys


def mesh2trigeom(meshfile, origin):
    '''
    Convert trelis output mesh file to faultgeom format

    Input:
        meshfile = output of trelis/cubit in hdf5 format
        origin   = [lat, lon]
    Output:
    '''
    h5       = h5py.File(meshfile, 'r')
    e        = h5['coordx'][()]
    n        = h5['coordy'][()]
    u        = h5['coordz'][()]
    vertex   = np.vstack((e,u,n)).T
    mapelem  = h5['connect1'][()]

    with open('element', 'w') as efid:
        efid.write('# idx1   idx2   idx3   slip_ss_mm  slip_ds_mm   slip_ts_mm\n')
        for i in range(len(mapelem)):
            efid.write('{} {} {} {} {} {}\n'.format(mapelem[i,0], mapelem[i,1], mapelem[i,2], 0.0, 1.0, 0.0))
    
    print(origin)
    with open('vertex', 'w') as vfid:
        vfid.write('# Idx   Lat    Lon    Depth(km)\n')
        for i in range(len(u)):
            llh = gt.utm2llh([n[i], e[i]], origin)
            print(n[i], e[i], origin)
            vfid.write('{:5d} {:12.4f} {:12.4f} {:12.4f}\n'.format(i+1, llh[0,0], llh[0,1], -u[i]))


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('trelis2faultgeom.py <mesh file> <ref lat> <ref lon>')
        sys.exit()

    meshfile = sys.argv[1]
    lat      = float(sys.argv[2])
    lon      = float(sys.argv[3])
    origin   = [lat, lon]
    mesh2trigeom(meshfile, origin)
