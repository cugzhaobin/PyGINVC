#!/usr/bin/env python
# Plot 3D triangular element surface with or without slip 
# Written by Zhao Bin, Institute of Seismology, CEA. April 24, 2017.
# Modified by Zhao Bin, adding aftershock catalog as input parameter. April 11, 2018

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
import os, argparse
from pyginvc.Geometry.Triangle import Triangle

def plot_gps(ax, gpsfile, scale=3000):
    """Plot aftershocks if the file exists."""
    if os.path.exists(gpsfile):
        gps = np.genfromtxt(gpsfile)
        ax.quiver(gps[:,0], gps[:,1], np.zeros(len(gps)), gps[:,2], gps[:,3], np.zeros(len(gps)), length=0.001)
    else:
        print('Input gps file does not exist!')


def main():
    parser = argparse.ArgumentParser(description="plot 3D slip distribution.")
    parser.add_argument('--vertex', type=str, required=True, help='index lat lon dep(up positive)')
    parser.add_argument('--element', type=str, required=True, help='indx1 indx2 indx3 strike-slip dip-slip open-slip')
    parser.add_argument('--aftershock', type=str, required=False, default="", help='aftershock data')
    parser.add_argument('--coordtype', type=str, required=False, default='llh', help='llh/enu')
    parser.add_argument('--azimuth', type=float, required=False, default=-40, help='azimuth in degree')
    parser.add_argument('--elevation', type=float, required=False, default=40, help='elevation in degree')
    parser.add_argument('--sarfile', type=str, required=False, default="", help='InSAR file')
    parser.add_argument('--gpsfile', type=str, required=False, default="", help='GPS file')
    parser.add_argument('--savefig', type=bool, required=False, default=False, help='Y/N save figure')
    args = parser.parse_args()

    # get parameters
    vertexfile  = args.vertex
    elementfile = args.element
    afsfile     = args.aftershock
    azimuth     = args.azimuth
    elevation   = args.elevation

    # read in vertex and elements
    tri         = Triangle(vertexfile, elementfile)
    vertex_llh  = tri.vertex_llh
    element     = tri.element

    vertex_llh[:,2] = -np.abs(vertex_llh[:,2])

    # read in slip
    slipmodel = np.genfromtxt(elementfile)
    n,m       = slipmodel.shape
    #slipmodel[:,4] = 0.0
    if m == 6:
        slip  = np.sqrt(slipmodel[:,3]**2+slipmodel[:,4]**2)
        nslip = slip/slip.max()

    fig = plt.figure()
    ax  = fig.add_subplot(projection='3d')
    ax.view_init(elev=elevation, azim=azimuth)
    ax.scatter(vertex_llh[:,1], vertex_llh[:,0], vertex_llh[:,2], s=0.1)
    for i in range(len(element)):
        indx = element[i,:]-1
        lon  = np.array(vertex_llh[indx, 1])
        lat  = np.array(vertex_llh[indx, 0])
        dep  = np.array(vertex_llh[indx, 2])
        verts= [list(zip(lon, lat, dep))]
        subfault = Poly3DCollection(verts)
        if m == 6:
            subfault.set_facecolor(plt.cm.jet(nslip[i]))
        subfault.set_linewidth(1.5)
        ax.add_collection3d(subfault)
    #   ax.text(np.mean(lon), np.mean(lat), np.mean(dep), str(i), color='black', fontsize=5)

    indx1 = element[:,0]-1
    s   = ax.scatter(vertex_llh[indx1,1], vertex_llh[indx1,0], c=slip, cmap=plt.cm.jet, s=1e-5, lw=0)
    cb = plt.colorbar(s, shrink=0.8, pad=-0.07)

    # plot aftershock
    if os.path.exists(afsfile) == False:
        print(' Input aftershock catalog does not exist!')
    else:
        afs = np.genfromtxt(afsfile)
        ax.scatter(afs[:,0], afs[:,1], -afs[:,2])
        
    # plot InSAR data
    if os.path.exists(args.sarfile) == False:
        print(' Input InSAR file does not exist!')
    else:
        sar = np.genfromtxt(args.sarfile)
        s   = ax.scatter(sar[:,0], sar[:,1], np.zeros(len(sar)), c=sar[:,2], cmap=plt.cm.rainbow)
        plt.colorbar(s, shrink=0.8, pad=0.1, label='LOS (mm)')

    # plot GPS data
    if os.path.exists(args.gpsfile):
        plot_gps(ax, args.gpsfile)

    cb.set_label('Total slip (mm)')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth (km)')
    
    if bool(args.savefig) == True:
        plt.savefig('slip3d.pdf', format='pdf')
    plt.show()



if __name__ == '__main__':
    main()
