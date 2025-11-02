#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to input file of the VISCO1D software
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Ocb 4, 2016 
# Rewritten by Zhao Bin, April 6, 2021

from pyginvc.Geometry.Fault import Fault
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to VISCO1D input file.")
    parser.add_argument('--faultfile', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    args    = parser.parse_args()

    fltfile = args.faultfile
    flt     = Fault(fltfile, 1, 1, False)
    flt.FaultGeom2VISCO1D()

if __name__ == '__main__':
    main()
