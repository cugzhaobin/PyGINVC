#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA. July 10, 2017
# plot 3D slip distribution using Ouput file from DISMODEL/GINVPY
# MOD by Zhao Bin, add option of NEU/LLH coordiante system
# MOD by Zhao Bin, add option to read aftershock catalog. Sep 12, 2017
# MOD by Zhao Bin, add option to text the index of each patch
# MOD by Zhao Bin, plot surface velocity. Feb 27, 2019.
# MOD by Zhao Bin, change zip() to list(zip()). Sep. 2, 2019.


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from pyginvc.Geometry.Fault import Fault
import numpy as np
import os, argparse


cdict = {'red': ((0., 1, 1),
                 (0.10, 1, 1),
                 (0.20, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.10, 1, 1),
                   (0.20, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.15, 1, 1),
                  (0.20, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}
whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)

def main():
    parser = argparse.ArgumentParser(description="plot 3D slip distribution.")
    parser.add_argument('--faultfile', type=str, required=True, help='faultgeom format')
    parser.add_argument('--aftershock', type=str, required=False, default="", help='aftershock data')
    parser.add_argument('--coordtype', type=str, required=False, default='llh', help='llh/enu')
    parser.add_argument('--showindex', type=bool, required=False, default=False, help='False/True')
    parser.add_argument('--azimuth', type=float, required=False, default=-40, help='azimuth in degree')
    parser.add_argument('--elevation', type=float, required=False, default=40, help='elevation in degree')
    args = parser.parse_args()
    
    faultfile         = args.faultfile
    afsfile           = args.aftershock
    coordtype         = args.coordtype
    showtext          = args.showindex
    azim              = args.azimuth
    elev              = args.elevation
    
    geom              = np.genfromtxt(faultfile)
    nelem             = len(geom)
    flt               = Fault(faultfile, nelem, 1, False)
#   geom              = flt.LoadFaultGeom(faultfile)
#   [disgeom, origin] = flt.SetOrigin(geom)
#   felem             = flt.FaultGeom2AllVertex(geom, disgeom)
    felem             = flt.felem
    slip              = felem[:,7]/felem[:,7].max()
    
    if coordtype == 'llh':
        x                 = felem[:, [9,12,15,18]]
        y                 = felem[:, [10,13,16,19]]
        z                 = felem[:, [11,14,17,20]]
    elif coordtype == 'enu':
        x                 = felem[:, [21,24,27,30]]
        y                 = felem[:, [22,25,28,31]]
        z                 = felem[:, [23,26,29,32]]
    
    fig               = plt.figure()
    ax                = fig.add_subplot(projection='3d')
    ax.view_init(elev=elev, azim=azim)
    ax.scatter(x, y, z, s=0.0)
    
    for i in range(len(felem)):
        if coordtype == 'llh':
            indx = [9, 12, 15, 18]
            lon  = np.array(felem[i, indx])
            indx = [10, 13, 16, 19]
            lat  = np.array(felem[i, indx])
            indx = [11, 14, 17, 20]
            dep  = np.array(felem[i, indx])
            verts= [list(zip(lon, lat, dep))]
            subfault = Poly3DCollection(verts, facecolors='w', linewidths=1, alpha=0.3)
        elif coordtype == 'enu':
            indx = [21, 24, 27, 30]
            e    = np.array(felem[i, indx])
            indx = [22, 25, 28, 31]
            n    = np.array(felem[i, indx])
            indx = [23, 26, 29, 32]
            u    = np.array(felem[i, indx])
            verts= [list(zip(e, n, u))]
            subfault = Poly3DCollection(verts, facecolors='w', linewidths=1, alpha=0.3)
        subfault.set_facecolor(plt.cm.jet(slip[i]))
        subfault.set_linewidth(5)
        ax.add_collection3d(subfault)
        if coordtype == 'llh' and showtext:
            ax.text(np.mean(lon), np.mean(lat), np.mean(dep), str(i), color='red', fontsize=5)
        elif coordtype == 'neu' and showtext:
            ax.text(np.mean(n), np.mean(e), np.mean(u), str(i), color='red', fontsize=5)
    
    
    # plot aftershock/velocity/displacement
    if os.path.exists(afsfile) == False:
        print(' Input aftershock catalog does not exist!')
    else:
        afs = np.genfromtxt(afsfile)
        if afs.shape[1] >= 8:
            ax.scatter(afs[:,0], afs[:,1], afs[:,2], s=4, marker='s')
            ax.quiver(afs[:,0], afs[:,1], 0.0*afs[:,0], afs[:,2], afs[:,3], 0.0*afs[:,0], length=0.01, arrow_length_ratio=0.03)
        else:
            ax.scatter(afs[:,0], afs[:,1], afs[:,2], s=10, marker='s')
    
    
#   s  = ax.scatter(geom[:,4], geom[:,3],c=felem[:,7]/1000.,cmap=plt.cm.jet, s=0.00001,lw=0)
#   cb = plt.colorbar(s, shrink=0.6, pad=0.1)
#   cb.set_label('Total slip (m)', fontsize=12)
    if coordtype == 'llh':
        ax.set_xlabel('Longitude', fontsize=12, labelpad=12)
        ax.set_ylabel('Latitude', fontsize=12, labelpad=18)
    elif coordtype == 'enu':
        ax.set_xlabel('East (km)')
        ax.set_ylabel('North (km)')
    ax.set_zlabel('Depth (km)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('slip.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('slip.jpg', format='jpg', dpi=600)
    plt.show()


if __name__ == '__main__':
    main()
