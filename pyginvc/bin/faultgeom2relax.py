#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to input file of the RELAX software
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Ocb 4, 2016 
# Rewritten by Zhao Bin, April 6, 2021

from pyginvc.Geometry.Fault import Fault
import sys, argparse

def main(args):
    fltfile = args.fault_file
    flt     = Fault(fltfile, 1, 1, False)
    flt.FaultGeom2RELAX()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to RELAX input file.")
    parser.add_argument('--fault_file', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    args = parser.parse_args()
    main(args)
