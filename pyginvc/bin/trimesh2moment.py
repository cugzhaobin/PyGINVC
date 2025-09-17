#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to input file of the Coulomb software
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Ocb 4, 2016 
# Rewritten by Zhao Bin, April 6, 2021

from pyginvc.Geometry.Triangle import Triangle
import sys, argparse

def main(args):
    print(args)
    vertex       = args.vertex
    element      = args.element
    shearmodulus = args.shear_modulus
    flt          = Triangle(vertex, element)
    Mo, Mw       = flt.Moment(shearmodulus=shearmodulus)
    print(Mo, Mw)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute moment and magnitude from triangular dislocation model.")
    parser.add_argument('--vertex', type=str, required=True, help='')
    parser.add_argument('--element', type=str, required=True, help='idx1 idx2 idx3 strike_slip dip_slip ex_slip in mm')
    parser.add_argument('--shear_modulus', type=float, required=False, default=3e10, help='shear modulus in Pa')
    args = parser.parse_args()
    main(args)
