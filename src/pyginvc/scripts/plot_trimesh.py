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


def tri_coordinate_system(tri_enu):
    '''
    Unit normal, strike and dip vectors of one triangular element.

    tri_enu has columns east/north/up (km).  This mirrors
    pyginvc.Greens.tdcalc.build_tri_coordinate_system, the basis used by the
    triangular Green's functions, where the physical slip vector is

        slip = ss * Vstrike + ds * Vdip  (+ opening * Vnorm)

    so the arrows plotted here point in the same direction as the slip used in
    the forward/inverse model.  Returns None for a degenerate element.
    '''
    vnorm = np.cross(tri_enu[1] - tri_enu[0], tri_enu[2] - tri_enu[0])
    norm  = np.linalg.norm(vnorm)
    if norm == 0:
        return None
    vnorm = vnorm / norm

    up      = np.array([0.0, 0.0, 1.0])
    vstrike = np.cross(up, vnorm)
    if np.linalg.norm(vstrike) == 0:
        # Horizontal element: strike is taken as north (or south if the
        # normal points downwards), as in the TD Green's functions.
        vstrike = np.array([0.0, 1.0, 0.0]) * vnorm[2]
        if tri_enu[0][2] > 0:
            vstrike = -vstrike
    vstrike = vstrike / np.linalg.norm(vstrike)
    vdip    = np.cross(vnorm, vstrike)
    return vnorm, vstrike, vdip


def plot_slip_vectors(ax, vertex_llh, vertex_enu, element, slip,
                      slip_scale=0.6, slip_head=0.18, slip_head_angle=15.0):
    '''
    Draw one in-plane slip arrow at the centroid of each triangular element.

    vertex_llh : (n, 3) array, lat/lon/depth with depth already negated for
                 plotting (up positive).
    vertex_enu : (n, 3) array, east/north/up in km.  z is negated internally
                 so this routine can be fed the raw Triangle.vertex_enu.
    element    : (m, 3) array of 1-based vertex indices.
    slip       : (m, 3) array of strike-slip / dip-slip / opening.
    slip_scale : longest arrow as a fraction of the median element size.
    slip_head  : arrow-head length as a fraction of the arrow length.
    slip_head_angle : half-angle (degrees) between each barb and the shaft.
    '''
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    slip = np.atleast_2d(np.asarray(slip, dtype=float))
    ss   = np.nan_to_num(slip[:, 0])
    ds   = np.nan_to_num(slip[:, 1])
    mag  = np.hypot(ss, ds)
    if mag.size == 0 or mag.max() <= 0:
        return

    # ENU vertices with z positive upwards (Triangle stores depth as +down).
    enu = np.array(vertex_enu, dtype=float, copy=True)
    enu[:, 2] = -enu[:, 2]

    centres = []
    sizes   = []
    for i, indx in enumerate(element):
        idx  = indx - 1
        tri  = enu[idx]
        basis = tri_coordinate_system(tri)
        if basis is None:
            continue
        vnorm, vstrike, vdip = basis

        # Physical element size (km): mean edge length.
        edges = (np.linalg.norm(tri[1] - tri[0]),
                 np.linalg.norm(tri[2] - tri[1]),
                 np.linalg.norm(tri[0] - tri[2]))
        size = np.mean(edges)
        if size <= 0 or mag[i] <= 0:
            continue

        vec = ss[i] * vstrike + ds[i] * vdip
        vlen = np.linalg.norm(vec)
        if vlen == 0:
            continue
        v_hat = vec / vlen
        # In-plane direction perpendicular to the slip vector (for the barbs).
        w_hat = np.cross(vnorm, v_hat)
        wnorm = np.linalg.norm(w_hat)
        if wnorm == 0:
            continue
        w_hat = w_hat / wnorm

        # Convert the ENU directions (per km of arrow) to plot coordinates
        # (lon, lat, z=-depth): east/north km -> degrees, up km -> plot z.
        llh      = vertex_llh[idx]
        lat_c    = np.mean(llh[:, 0])
        lon_c    = np.mean(llh[:, 1])
        z_c      = np.mean(llh[:, 2])
        km_per_lon = 111.32 * max(np.cos(np.deg2rad(lat_c)), 1e-3)
        p_v = np.array([v_hat[0] / km_per_lon, v_hat[1] / 111.32, v_hat[2]])
        p_w = np.array([w_hat[0] / km_per_lon, w_hat[1] / 111.32, w_hat[2]])

        centres.append((np.array([lon_c, lat_c, z_c]), mag[i], p_v, p_w))
        sizes.append(size)

    if not centres or not sizes:
        return

    arrow_scale = slip_scale * np.median(sizes) / mag.max()
    head_frac   = float(slip_head)
    head_ang    = np.deg2rad(float(slip_head_angle))

    segments = []
    for centre, m, p_v, p_w in centres:
        length = m * arrow_scale
        if length <= 0:
            continue
        tip  = centre + length * p_v
        head = head_frac * length
        back = tip - head * np.cos(head_ang) * p_v
        side = head * np.sin(head_ang)

        segments.append([centre, tip])
        segments.append([tip, back + side * p_w])
        segments.append([tip, back - side * p_w])

    if segments:
        ax.add_collection3d(
            Line3DCollection(segments, colors='k', linewidths=0.6))


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
    parser.add_argument('--showslip', action='store_true', help='Draw slip vectors on each element')
    parser.add_argument('--slipscale', type=float, default=0.6, help='longest slip arrow as a fraction of the median element size')
    parser.add_argument('--sliphead', type=float, default=0.18, help='slip arrow-head length as a fraction of the arrow')
    parser.add_argument('--slipheadangle', type=float, default=15.0, help='slip arrow-head half-angle in degrees')
    args = parser.parse_args()

    # get parameters
    vertexfile  = args.vertex
    elementfile = args.element
    afsfile     = args.aftershock
    azimuth     = args.azimuth
    elevation   = args.elevation

    # read in vertex and elements
    tri         = Triangle(vertexfile, elementfile)
    tri.load_fault()
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

    # plot slip vectors (one arrow per triangular element)
    if args.showslip:
        plot_slip_vectors(ax, vertex_llh, tri.vertex_enu, element, tri.slip,
                          slip_scale=args.slipscale,
                          slip_head=args.sliphead,
                          slip_head_angle=args.slipheadangle)

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
