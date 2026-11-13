#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA. July 10, 2017
# plot 3D slip distribution using Output file from DISMODEL/GINVPY
# MOD by Zhao Bin, add option of NEU/LLH coordinate system
# MOD by Zhao Bin, add option to read aftershock catalog. Sep 12, 2017
# MOD by Zhao Bin, add option to text the index of each patch
# MOD by Zhao Bin, plot surface velocity. Feb 27, 2019.
# MOD by Zhao Bin, change zip() to list(zip()). Sep. 2, 2019.

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
import argparse, os
import matplotlib
from pyginvc.Geometry.Patch import Fault


# Custom colormap
cdict = {
    'red': ((0., 1, 1), (0.10, 1, 1), (0.20, 0, 0), (0.66, 1, 1), (0.89, 1, 1), (1, 0.5, 0.5)),
    'green': ((0., 1, 1), (0.10, 1, 1), (0.20, 0, 0), (0.375, 1, 1), (0.64, 1, 1), (0.91, 0, 0), (1, 0, 0)),
    'blue': ((0., 1, 1), (0.15, 1, 1), (0.20, 1, 1), (0.34, 1, 1), (0.65, 0, 0), (1, 0, 0))
}
whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet', cdict, 256)

def main():
    parser = argparse.ArgumentParser(description="Plot 3D slip distribution.")
    parser.add_argument('--faultfile', type=str, required=True, help='Fault geometry file')
    parser.add_argument('--aftershock', type=str, required=False, default="", help='Aftershock data file')
    parser.add_argument('--coordtype', type=str, required=False, default='llh', help='Coordinate system: llh/enu')
    parser.add_argument('--showindex', type=bool, required=False, default=False, help='Show patch indices: False/True')
    parser.add_argument('--azimuth', type=float, required=False, default=-20, help='Azimuth in degrees')
    parser.add_argument('--elevation', type=float, required=False, default=40, help='Elevation in degrees')
    parser.add_argument('--sarfile', type=str, required=False, default="", help='InSAR file')
    parser.add_argument('--gpsfile', type=str, required=False, default="", help='GPS file')
    parser.add_argument('--savefig', type=str, required=False, default="", help='Save the figure to this path (e.g. faultgeom.png)')
    parser.add_argument('--scale', type=float, required=False, default=1000, help='scale for vector')
    parser.add_argument('--nomap', action='store_true', help='Disable the cartopy map view')
    parser.add_argument('--showslip', action='store_true', help='Draw slip vectors on each patch')
    parser.add_argument('--slipscale', type=float, default=0.6, help='longest slip arrow as a fraction of the median patch size')
    parser.add_argument('--sliphead', type=float, default=0.18, help='slip arrow-head length as a fraction of the arrow')
    parser.add_argument('--slipheadangle', type=float, default=15.0, help='slip arrow-head half-angle in degrees')
    args = parser.parse_args()

    faultfile = args.faultfile
    
    afsfile   = args.aftershock
    gpsfile   = args.gpsfile
    sarfile   = args.sarfile
    coordtype = args.coordtype
    showtext  = args.showindex
    azim      = args.azimuth
    elev      = args.elevation
    scale     = args.scale

    # The geometry file is already discretised, so doSubFault=False and the
    # number of patches is simply the number of rows in the file.
    ngeom     = len(np.atleast_2d(np.genfromtxt(faultfile, comments='#')))
    flt       = Fault(faultfile, 1, ngeom, False)
    flt.load_fault()
    flt.plot_faultgeom(afsfile=afsfile, gpsfile=gpsfile, sarfile=sarfile,
                       coordtype=coordtype, show_idx=showtext,
                       show_map=not args.nomap, scale=scale,
                       show_slip=args.showslip, slip_scale=args.slipscale,
                       slip_head=args.sliphead,
                       slip_head_angle=args.slipheadangle,
                       azimuth=azim, elevation=elev,
                       savefig=args.savefig)


if __name__ == '__main__':
    main()
