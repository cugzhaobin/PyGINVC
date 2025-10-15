#!/usr/bin/env python
import numpy as np
from numpy import sin, cos
import argparse


def main(args):
    fltfile = args.faultfile
    angle   = args.angle
    beta    = np.deg2rad(angle)
    R       = np.array([[cos(beta), -sin(beta)],
                        [sin(beta), cos(beta)]])
    dat  =  np.genfromtxt(fltfile)
    
    if dat.ndim == 1:
        dat = np.expand_dims(dat, axis=0)
    if dat.shape[1] == 10:
        slip = dat[:,[7,8]]
        out  = R.dot(slip.T).T
        out2 = np.column_stack((dat[:,0:7], out, dat[:,9]))
        np.savetxt('new.faultgeom', out2, fmt='%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.3f %10.3f %10.3f')
    elif dat.shape[1] == 6:
        slip = dat[:,[3,4]]
        out  = R.dot(slip.T).T
        rake = np.rad2deg(np.arctan2(out[:,1], out[:,0]))
        out2 = np.column_stack((dat[:,0:3], out, dat[:,5]))
        np.savetxt('slip_new.tri', out2, fmt='%4d %4d %4d %10.3f %10.3f %10.3f')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to VTK format.")
    parser.add_argument('--faultfile', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    parser.add_argument('--angle', default=0, help='degree', type=float)
    args = parser.parse_args()
    main(args)
