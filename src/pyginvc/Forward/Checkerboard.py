#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. August 28, 2016
# Mod by Zhao Bin, Dec. 7, 2018. replace print() with print; with open as fid

import numpy as np
import random, sys

def sync_gpsfile(gpsfile, mean, std):
    '''
    Generate GPS data file.

    Input:
        gpsfile = GPS data file in GMT format
        mean    = mean value
        std     = standard 
    '''
    data = np.loadtxt(gpsfile, usecols=[0,1,2,3,4,5,6], comments='#')
    with open('gps_syn.gmtvec', 'w') as fid:
        for i in range(len(data)):
            sige = random.normalvariate(mean, std)
            sign = random.normalvariate(mean, std)
            fid.write("%10.4f\t%10.4f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f\n" %(data[i,0], data[i,1], data[i,2], data[i,3], sige, sign, 0.0))


def make_checkerborad(faultgeom, Ns, Nd, ss, dd):
    '''
    Create checkerborad input for subfaults
    The example is Nd=4; Ns=6; dd=2; ss=3				   
     _________________________						   
     |	 |   |   |	 |   |   |					   
     | 1 | 1 | 1 | 0 | 0 | 0 |                                             
     |___|___|___|___|___|___|                                             
     |	 |   |   |	 |   |   |                                         
     | 1 | 1 | 1 | 0 | 0 | 0 |                                             
     |___|___|___|___|___|___|                                             
     |	 |   |   |	 |   |   |					   
     | 0 | 0 | 0 | 1 | 1 | 1 |                                             
     |___|___|___|___|___|___|                                             
     |	 |   |   |	 |   |   |                                         
     | 0 | 0 | 0 | 1 | 1 | 1 |                                             
     |___|___|___|___|___|___|                                             
  									   
    '''

#   fg = flt.LoadFaultGeom(faultgeom)
    fg = np.loadtxt(faultgeom)

    if len(fg) != Ns*Nd:
        print('please check FaultGeom file or [Ns, Nd]')
        sys.exit()

    filename = "FaultGeom_"+str(ss)+"by"+str(dd)

    with open(filename, 'w') as fid:
        for ii in range(0, Nd):
            id  = np.ceil((ii+1)/dd)
            md  = np.mod(id,2)
            for jj in range(0, Ns):
                ns  = np.ceil((jj+1)/ss)
                ms  = np.mod(ns, 2)
                index = ii*Ns + jj
                if md == ms:
                    fid.write("%f %f %f %f %f %f %f %f %f %f\n" %(fg[index,0], fg[index,1], fg[index,2], fg[index,3], fg[index,4], fg[index,5], fg[index,6], 0, 0, 0))
                else:
                    fid.write("%f %f %f %f %f %f %f %f %f %f\n" %(fg[index,0], fg[index,1], fg[index,2], fg[index,3], fg[index,4], fg[index,5], fg[index,6], 1000, 0, 0))
