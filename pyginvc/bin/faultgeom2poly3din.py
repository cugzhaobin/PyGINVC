#!/usr/bin/env python
# This program is used to create input file for POLY3D program
# This program is suitable for rectangualr dislocation
# The Input file for this script is FaultGeom, which contain coseismic rupture
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Oct 28, 2016 

# import libs
import numpy as np
import sys, os
import argparse
import poly3d
import pyginvc.libs.geotools as gt

def main(args):
    # get the input parameter
    slipmodel         = args.faultfile
    origin            = args.origin
    print(origin)

    # get uniqe vertex, elem index, disp_element and origin
    vertex, elem, disp = poly3d.vertex_element_from_faultgeom(slipmodel, origin)

    bctype = [args.bctype for _ in range(len(disp))]
    # open and creat input file for POLY3D
    poly3d.gen_poly3d_input(vertex, elem, disp, bctype=bctype)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to Poly3D input file.")
    parser.add_argument('--faultfile', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    parser.add_argument('--origin', default=[], help='[lat, lon]', nargs=2, type=float)
    parser.add_argument('--bctype', default='bbb', help='bbb/ttt', type=str)
    args = parser.parse_args()
    main(args)
