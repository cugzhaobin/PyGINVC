#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to VTK format
# Rewritten by Zhao Bin, April 12, 2021

from pyginvc.Geometry.Fault import Fault
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to VTK format.")
    parser.add_argument('--faultfile', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    parser.add_argument('--origin', default=[], help='[lat, lon]', nargs=2, type=float)
    parser.add_argument('--scale', default=1.0, type=float)
    args    = parser.parse_args()
    fltfile = args.faultfile
    origin  = args.origin
    scale   = args.scale
    flt     = Fault(fltfile, 1, 1, False, origin=origin)
    flt.FaultGeom2VTK(scale=1)

if __name__ == '__main__':
    main()
