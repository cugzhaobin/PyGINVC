#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to input file of the VISCO2PT5D software
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Ocb 4, 2016 
# Rewritten by Zhao Bin, April 6, 2021

from pyginvc.Geometry.Fault import Fault
import sys, argparse
import numpy as np

def main(args):
    fltfile = args.faultfile
    dat     = np.genfromtxt(fltfile)
    flt     = Fault(fltfile, len(dat), 1, False)
    flt.FaultGeom2VISCO2PT5D()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to VISCO2.5D input file.")
    parser.add_argument('--faultfile', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    args = parser.parse_args()
    main(args)
