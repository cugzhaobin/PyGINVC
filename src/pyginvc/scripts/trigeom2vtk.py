#!/usr/bin/env python
# Convert triangular element file, output of IDSMODEL/GINVPY to VTK format
# Rewritten by Zhao Bin, April 29, 2021

from pyginvc.Geometry.Triangle import Triangle
import sys, argparse

def main(args):
    vertex  = args.vertex
    element = args.element
    origin  = args.origin
    flt     = Triangle(vertex, element, origin=origin)
    flt.Triangle2VTK()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert Triangle file to VTK format.")
    parser.add_argument('--vertex', type=str, required=True, help='index lat lon dep(up positive)')
    parser.add_argument('--element', type=str, required=True, help='indx1 indx2 indx3 strike-slip dip-slip open-slip')
    parser.add_argument('--origin', default=[], help='[lat, lon]', nargs=2, type=float)
    args = parser.parse_args()
    main(args)
