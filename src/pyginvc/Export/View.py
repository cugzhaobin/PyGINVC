#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 4 2016

import h5py, sys
import numpy as np
import matplotlib.pyplot as plt

from pyginvc.libs import geotools as gt
from pyginvc.Geometry.Fault import Fault
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_obs_mod(solutionfile, scale=100):
    '''
    Plot observation and prediction for GPS and InSAR data

    Parameters
    ----------
    solutionfile : string
        HDF5 file output from PyGINV.Inverse
    scale : float, optional
        scale of vector. The default is 100.

    Returns
    -------
    None.

    '''

    with h5py.File(solutionfile, 'r') as h5:
        #
        # For GPS stations
        #
        if 'obs' not in h5.keys():
            print('No group named obs.')
            sys.exit()
        else:
            if 'llh_gps' not in h5['obs'].keys():
                print('No dataset named llh_gps.')
                sys.exit()
        llh_gps = h5['obs']['llh_gps'][()]
        origin  = h5['flt']['origin'][()]
        enu     = np.zeros((len(llh_gps),2)) 
        for i in range(len(llh_gps)):
            enu[i] = gt.llh2localxy(llh_gps[i], origin)
        if len(llh_gps) == 0:
            print('No GPS data.')
            ndim  = 0
        else:
            ndim  = h5['obs']['ndim'][()]
            d_gps = h5['obs']['d_gps'][()]
            d_gps = d_gps.reshape(len(llh_gps), ndim)
            m_gps = h5['sol']['dhat'][0:ndim*len(llh_gps)]
            m_gps = m_gps.reshape(len(llh_gps), ndim)
            plt.figure(figsize=(10,5))
            plt.subplot(1,2,1)
            plt.quiver(llh_gps[:,1], llh_gps[:,0], d_gps[:,0], d_gps[:,1], color='b', scale=scale)
            plt.quiver(llh_gps[:,1], llh_gps[:,0], m_gps[:,0], m_gps[:,1], color='r', scale=scale)
            if ndim == 3:
                plt.scatter(llh_gps[:,1], llh_gps[:,0], c=d_gps[:,2], cmap=plt.cm.jet, s=50, linewidths=0.1, edgecolors='black')
                plt.scatter(llh_gps[:,1], llh_gps[:,0], c=m_gps[:,2], cmap=plt.cm.jet, s=15, linewidths=0.1, edgecolors='black')
            plt.title('Observation VS Prediction')
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.subplot(1,2,2)
            plt.quiver(llh_gps[:,1], llh_gps[:,0], d_gps[:,0]-m_gps[:,0], d_gps[:,1]-m_gps[:,1], color='b', scale=scale)
            if ndim == 3:
                plt.scatter(llh_gps[:,1], llh_gps[:,0], c=d_gps[:,2]-m_gps[:,2], cmap=plt.cm.jet,linewidths=0.1, edgecolors='black')
            plt.title('Residual')
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.savefig('gps_mod_obs.pdf', format='pdf')
            plt.show()

        # 
        # For InSAR
        #
        if 'llh_sar' in h5['obs'].keys():
            llh_sar = h5['obs']['llh_sar'][()]
            if len(llh_sar) == 0:
                print(' No InSAR points')
            else:
                d_sar = h5['obs']['d_sar'][()]
                m_sar = h5['sol']['dhat'][ndim*len(llh_gps):ndim*len(llh_gps)+len(llh_sar)]
                plt.figure(figsize=(10,5))
                plt.subplot(1,2,1)
                plt.scatter(llh_sar[:,1], llh_sar[:,0], c=d_sar, cmap=plt.cm.jet, s=2)
                plt.xlabel('Longitude')
                plt.ylabel('Latitude')
                plt.title('InSAR LOS Observation')
                cb    = plt.colorbar()
                cb.set_label('Lin changes (mm)')
                plt.subplot(1,2,2)
                plt.scatter(llh_sar[:,1], llh_sar[:,0], c=m_sar, cmap=plt.cm.jet, s=2)
                plt.xlabel('Longitude')
                plt.ylabel('Latitude')
                plt.title('InSAR LOS Prediction ')
                cb    = plt.colorbar()
                cb.set_label('Line of sight changes (mm)')
                plt.savefig('sar_mod_obs.pdf', format='pdf')
                plt.show()


def plot_slip_3d(solutionfile, elevation=40, azimuth=-81, coordtype='llh', show=True):
    '''
    Plot 3D finite slip distribution.

    Parameters
    ----------
    solutionfile : string
        HDF5 file output from PyGINV program
    elevation : float, optional
        look angle
    azimuth : float, optional
        look angle. The default is -81.
    coordtype : string, optional
        'llh/enu'. The default is 'llh'.

    Returns
    -------
    None.

    '''

    with h5py.File(solutionfile, 'r') as h5:
        if 'flt' not in h5.keys():
            print('No fault information for plotting.')
            sys.exit()
        else:
            fig      = plt.figure()
            ax       = fig.add_subplot(111, projection='3d')
            ax.view_init(elev=elevation, azim=azimuth)
            
            if h5['flt/flttype'][()] == 'rectangle':
                geom     = h5['flt/geom'][()]
                dis_geom = h5['flt/dis_geom'][()]
                slip     = h5['sol/slip'][()]
                slip     = slip[0,0:len(geom)*3]
                slip     = slip.reshape(len(geom), 3)
                ts       = np.sqrt(slip[:,0]**2+slip[:,1]**2)
                ts1      = ts/ts.max()
                felem    = Fault.FaultGeom2AllVertex(geom, dis_geom, origin=[])
    
                if coordtype == 'llh':
                    idx  = np.array([9,12,15,18])
                    x    = felem[:,idx]
                    y    = felem[:,idx+1]
                    z    = felem[:,idx+2]
                elif coordtype == 'enu':
                    idx  = [21,24,27,30]
                    x    = felem[:,idx]
                    y    = felem[:,idx+1]
                    z    = felem[:,idx+2]
    
                for i in range(len(felem)):
                    verts= [list(zip(x[i], y[i], z[i]))]
                    subfault = Poly3DCollection(verts, facecolors='w', linewidths=1, alpha=0.3)
                    subfault.set_facecolor(plt.cm.jet(ts1[i]))
                    subfault.set_linewidth(0.05)
                    ax.add_collection3d(subfault)
                s    = ax.scatter(geom[:,4], geom[:,3], c=ts/1000., cmap=plt.cm.jet, s=0.00001, lw=0)
                cb   = plt.colorbar(s, shrink=0.8, pad=-0.07)
                cb.set_label('Total slip(m)')
                if coordtype == 'llh':
                    ax.set_xlabel('Longitude')
                    ax.set_ylabel(':atitude')
                elif coordtype == 'enu':
                    ax.set_xlabel('East (km)')
                    ax.set_ylabel('North (km)')
                ax.set_zlabel('Depth (km)')
                plt.show()
            elif h5['flt/flttype'][()] == 'triangle':
                element  = h5['flt/element'][()]-1
                slip     = h5['sol/slip'][()]
                slip     = slip[0,0:len(element)*3]
                slip     = slip.reshape(len(element), 3)
                ts       = np.sqrt(slip[:,0]**2+slip[:,1]**2)
                ts1      = ts/ts.max()
                print(ts1)

                if coordtype == 'llh':
                    # lat, lon, down
                    vertex_llh = h5['flt/vertex_llh'][()]
                    idx = element[:,0]
                    ax.scatter(vertex_llh[:,1],
                               vertex_llh[:,0],
                              -vertex_llh[:,2],
                               s=0.1)
                    for i in range(len(element)):
                        indx = element[i,:]
                        lat  = np.array(vertex_llh[indx, 0])
                        lon  = np.array(vertex_llh[indx, 1])
                        dep  = np.array(-vertex_llh[indx, 2])
                        verts= [list(zip(lon, lat, dep))]
                        subfault = Poly3DCollection(verts)
                        subfault.set_facecolor(plt.cm.jet(ts1[i]))
                        subfault.set_linewidth(1.5)
                        ax.add_collection3d(subfault)

                    idx = element[:,0]
                    s   = ax.scatter(vertex_llh[idx,1], 
                                     vertex_llh[idx,0], 
                                     c=ts, cmap=plt.cm.jet, s=1e-5, lw=0)
                    cb = plt.colorbar(s, shrink=0.8, pad=-0.07)
                    cb.set_label('Total slip (mm)')
        
                    plt.xlabel('Longitude')
                    plt.ylabel('Latitude')
                    ax.set_zlabel('Depth (km)')
            
                elif coordtype in {'neu', 'enu', 'NEU', 'ENU'}:
                    vertex_enu = h5['flt/vertex_enu'][()]
                    idx = element[:,0]
                    s   = ax.scatter(vertex_enu[idx,1], 
                                     vertex_enu[idx,0], 
                                     -vertex_enu[idx,2], 
                                     c=slip, cmap=plt.cm.jet, s=1e-5, lw=0)
                    
                    for i in range(len(element)):
                        indx = element[i,:]
                        e    = np.array(vertex_enu[indx, 1])
                        n    = np.array(vertex_enu[indx, 1])
                        u    = np.array(vertex_enu[indx, 2])
                        verts= [list(zip(e, n, u))]
                        subfault = Poly3DCollection(verts)
                        subfault.set_facecolor(plt.cm.jet(ts1[i]))
                        subfault.set_linewidth(1.5)
                        ax.add_collection3d(subfault)
                    # plot color bar
                    cb = plt.colorbar(s, shrink=0.8, pad=-0.07)
                    cb.set_label('Total slip (m)')
                    
                    
                    plt.xlabel('East (m)')
                    plt.ylabel('North (m)')
                    ax.set_zlabel('Depth (m)')
                plt.show()
