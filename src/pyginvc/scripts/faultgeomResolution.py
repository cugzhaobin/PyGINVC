#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to VTK format
# Rewritten by Zhao Bin, April 12, 2021

import sys, argparse
import h5py
import numpy as np
from   scipy.linalg import inv, pinv

def main():
    parser = argparse.ArgumentParser(description="Compute model resolution matrix and output to FaultGeom format.")
    parser.add_argument('--faultfile', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    parser.add_argument('--green', required=True, type=str)
    parser.add_argument('--smooth', required=True, type=str)
    args    = parser.parse_args()

    fltfile = args.faultfile
    fltgeom = np.loadtxt(args.faultfile)
    G       = np.empty(0)
    Lap     = np.empty(0)
    with h5py.File(args.green) as h5:
        G   = h5['G'][()]
    with h5py.File(args.smooth) as h5:
        Lap = h5['lap'][()]

    Gg = (G.T).dot(G) + (Lap.T).dot(Lap)*0.0
    Gg = pinv(Gg).dot(G.T)
    R  = np.diag(Gg.dot(G))
    R  = R.reshape((Gg.shape[0]//3,3))
    newout = np.column_stack((fltgeom[:,0:7], R))
    np.savetxt('aaa', newout)

if __name__ == '__main__':
    main()
