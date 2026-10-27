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


def mesh2faultgeom(meshfile, origin):
    '''
    Convert trelis output mesh file to faultgeom format

    Input:
        meshfile = output of trelis/cubit in hdf5 format
        origin   = [lat, lon]
    Output:
    '''
    lon      = float(sys.argv[2])
    lat      = float(sys.argv[3])
    h5       = h5py.File(meshfile, 'r')
    e        = h5['coordx'][()]
    n        = h5['coordy'][()]
    u        = h5['coordz'][()]
    vertex   = np.vstack((e,u,n)).T
    mapelem  = h5['connect1'][()]-1
    origin   = [lat, lon]
    elements = []

    for i in range(len(mapelem)):
        patch = []
        for j in mapelem[i]:
            vertex = np.array([e[j], n[j], u[j]])
            patch.append(vertex)
        elements.append(patch)
    
    
    for i in range(len(elements)):
        patch = elements[i]
        v1    = patch[0]
        v2    = patch[1]
        v3    = patch[2]
        v4    = patch[3]
        length1 = np.sqrt((v1[0] - v2[0])**2 + (v1[1] - v2[1])**2)
        length2 = np.sqrt((v3[0] - v4[0])**2 + (v3[1] - v4[1])**2)
        length  = (length1+length2)/2
        width1  = np.sqrt((v1[0]-v4[0])**2 + (v1[1]-v4[1])**2 + (v1[2]-v4[2])**2)
        width2  = np.sqrt((v2[0]-v3[0])**2 + (v2[1]-v3[1])**2 + (v2[2]-v3[2])**2)
        width   = (width1+width2)/2

        # dip angle of the fault
        dh1   = v3[2]-v2[2]
        dh2   = v4[2]-v1[2]
        dh    = (dh1+dh2)/2
        dip  = np.rad2deg(math.asin(dh/width))
        dip  = np.abs(dip)+180
    
        llh1 = gt.utm2llh([v1[1], v1[0]], origin)
        llh2 = gt.utm2llh([v2[1], v2[0]], origin)
    
        fmt  = 3*"%10.3f"+4*"%10.4f"+3*"%15.3f"
        print(fmt %(width, np.abs(v1[2]), dip, llh2[0,0], llh2[0,1], llh1[0,0], llh1[0,1], 1.0, 1.0, 0.0))


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('trelis2faultgeom.py <mesh file> <ref lon> <ref lat>')
        sys.exit()

    meshfile = sys.argv[1]
    lon      = float(sys.argv[2])
    lat      = float(sys.argv[3])
    origin   = [lat, lon]
    mesh2faultgeom(meshfile, origin)
