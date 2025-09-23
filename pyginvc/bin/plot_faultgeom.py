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

def load_fault_data(faultfile):
    """Load fault geometry data from file."""
    geom  = np.genfromtxt(faultfile)
    nelem = len(geom)
    flt   = Fault(faultfile, nelem, 1, False)
    felem = flt.felem
    return felem

def plot_fault_patches(ax, felem, coordtype, showtext, slip):
    """Plot fault patches on the 3D axis."""
    for i in range(len(felem)):
        if coordtype == 'llh':
            lon = felem[i, [9, 12, 15, 18]]
            lat = felem[i, [10, 13, 16, 19]]
            dep = felem[i, [11, 14, 17, 20]]
            verts = [list(zip(lon, lat, dep))]
            ax.scatter(lon, lat, dep, s=0)
        elif coordtype == 'enu':
            e = felem[i, [21, 24, 27, 30]]
            n = felem[i, [22, 25, 28, 31]]
            u = felem[i, [23, 26, 29, 32]]
            verts = [list(zip(e, n, u))]
            ax.scatter(e, n, u, s=0)

        
        subfault = Poly3DCollection(verts, facecolors='w', linewidths=1, alpha=0.3)
        subfault.set_facecolor(plt.cm.jet(slip[i]))
        subfault.set_linewidth(5)
        ax.add_collection3d(subfault)

        if showtext:
            if coordtype == 'llh':
                ax.text(np.mean(lon), np.mean(lat), np.mean(dep), str(i), color='red', fontsize=5)
            elif coordtype == 'enu':
                ax.text(np.mean(n), np.mean(e), np.mean(u), str(i), color='red', fontsize=5)

def plot_aftershocks(ax, afsfile):
    """Plot aftershocks if the file exists."""
    if os.path.exists(afsfile):
        afs = np.genfromtxt(afsfile)
        ax.scatter(afs[:, 0], afs[:, 1], afs[:, 2], s=4, marker='s')
    else:
        print('Input aftershock catalog does not exist!')
        
def plot_sar(ax, sarfile):
    """Plot aftershocks if the file exists."""
    if os.path.exists(sarfile):
        sar = np.genfromtxt(sarfile)
        s   = ax.scatter(sar[:,0], sar[:,1], np.zeros(len(sar)), c=sar[:,2], cmap=plt.cm.rainbow)
        plt.colorbar(s, shrink=0.8, pad=0.1, label='LOS (mm)')
    else:
        print('Input sar file does not exist!')

def plot_gps(ax, gpsfile, scale=1000):
    """Plot aftershocks if the file exists."""
    if os.path.exists(gpsfile):
        gps = np.genfromtxt(gpsfile)
        ax.quiver(gps[:,0], gps[:,1], np.zeros(len(gps)), gps[:,2], gps[:,3], np.zeros(len(gps)))
    else:
        print('Input gps file does not exist!')

def main(args):
    faultfile = args.faultfile
    afsfile   = args.aftershock
    coordtype = args.coordtype
    showtext  = args.showindex
    azim      = args.azimuth
    elev      = args.elevation
    """
    # Load fault data
    felem = load_fault_data(faultfile)
    slip  = felem[:, 7] / felem[:, 7].max()
    
    if coordtype == 'llh':
        x = felem[:, [9,12,15,18]]
        y = felem[:, [10,13,16,19]]
        z = felem[:, [11,14,17,20]]
    elif coordtype == 'enu':
        x = felem[:, [21,24,27,30]]
        y = felem[:, [22,25,28,31]]
        z = felem[:, [23,26,29,32]]

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.view_init(elev=elev, azim=azim)
    ax.scatter(x, y, z, s=0.0)

    # Plot fault patches
    for i in range(len(felem)):
        if coordtype == 'llh':
            lon = felem[i, [9, 12, 15, 18]]
            lat = felem[i, [10, 13, 16, 19]]
            dep = felem[i, [11, 14, 17, 20]]
            verts = [list(zip(lon, lat, dep))]
            ax.scatter(lon, lat, dep, s=0)
        elif coordtype == 'enu':
            e = felem[i, [21, 24, 27, 30]]
            n = felem[i, [22, 25, 28, 31]]
            u = felem[i, [23, 26, 29, 32]]
            verts = [list(zip(e, n, u))]
            ax.scatter(e, n, u, s=0)

        
        subfault = Poly3DCollection(verts, facecolors='w', linewidths=1, alpha=0.3)
        subfault.set_facecolor(plt.cm.jet(slip[i]))
        subfault.set_linewidth(5)
        ax.add_collection3d(subfault)

        if showtext:
            if coordtype == 'llh':
                ax.text(np.mean(lon), np.mean(lat), np.mean(dep), str(i), color='red', fontsize=5)
            elif coordtype == 'enu':
                ax.text(np.mean(n), np.mean(e), np.mean(u), str(i), color='red', fontsize=5)

    # Plot aftershocks
    plot_aftershocks(ax, afsfile)
    
    # Plot InSAR
    plot_sar(ax, args.sarfile)
    
    # Plot GPS
    plot_gps(ax, args.gpsfile)
    
    # Set labels and save plot
    if coordtype == 'llh':
        ax.set_xlabel('Longitude', fontsize=12, labelpad=12)
        ax.set_ylabel('Latitude', fontsize=12, labelpad=18)
    elif coordtype == 'enu':
        ax.set_xlabel('East (km)')
        ax.set_ylabel('North (km)')
    ax.set_zlabel('Depth (km)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    s  = ax.scatter(felem[:,9], felem[:,10], c=felem[:,7]/1000.,cmap=plt.cm.jet, s=0.00001,lw=0)
    plt.colorbar(s, shrink=0.6, pad=0.1, label='Slip (m)')
    #   cb.set_label('Total slip (m)', fontsize=12)
    
    if args.savefig == True:
        plt.savefig('slip.pdf', format='pdf', bbox_inches='tight')
        plt.savefig('slip.jpg', format='jpg', dpi=600)
    plt.show()"""
    geom  = np.genfromtxt(faultfile)
    nelem = len(geom)
    flt   = Fault(faultfile, nelem, 1, False)
    flt.plot_faultgeom()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot 3D slip distribution.")
    parser.add_argument('--faultfile', type=str, required=True, help='Fault geometry file')
    parser.add_argument('--aftershock', type=str, required=False, default="", help='Aftershock data file')
    parser.add_argument('--coordtype', type=str, required=False, default='llh', help='Coordinate system: llh/enu')
    parser.add_argument('--showindex', type=bool, required=False, default=False, help='Show patch indices: False/True')
    parser.add_argument('--azimuth', type=float, required=False, default=-40, help='Azimuth in degrees')
    parser.add_argument('--elevation', type=float, required=False, default=40, help='Elevation in degrees')
    parser.add_argument('--sarfile', type=str, required=False, default="", help='InSAR file')
    parser.add_argument('--gpsfile', type=str, required=False, default="", help='GPS file')
    parser.add_argument('--savefig', type=bool, required=False, default=False, help='Y/N save figure')
    args = parser.parse_args()
    main(args)
