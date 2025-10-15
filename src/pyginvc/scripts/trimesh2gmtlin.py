#!/usr/bin/env python
# This script depends on PyUnicycle.Geometry
from pyginvc.Geometry.Triangle import Triangle
import argparse
import numpy as np

def main(args):
    vertexfile  = args.vertex
    elementfile = args.element
    coordtype   = args.coordtype

    tri = Triangle(vertexfile, elementfile)

    if coordtype == 'neu':
        x  = tri.vertex_enu
    if coordtype == 'llh':
        x  = tri.vertex_llh
    v  = tri.element-1

    with open('ds_slip.gmtlin', 'w') as fid:
        fid.write('# export from PyUnicycle\n')
        for i in range(len(v)):
            if coordtype == 'enu' or coordtype == 'enu':
                enu = tri.vertex_enu[tri.element[i]-1]
            elif coordtype == 'llh':
                enu = tri.vertex_llh[tri.element[i]-1]
                enu = enu[:,[1,0,2]]
                enu[:,2] = -enu[:,2]
            fid.write('> -Z%f\n' %(tri.slip[i,1]))
            fid.write('%f %f %f\n' %(enu[0,0], enu[0,1], enu[0,2]))
            fid.write('%f %f %f\n' %(enu[1,0], enu[1,1], enu[1,2]))
            fid.write('%f %f %f\n' %(enu[2,0], enu[2,1], enu[2,2]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert vertex/element file to gmtlin")
    parser.add_argument('--vertex', type=str, required=True, help='')
    parser.add_argument('--element', type=str, required=True, help='')
    parser.add_argument('--coordtype', type=str, required=True, help='llh/enu')
    args = parser.parse_args()
    main(args)
