#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 7 2016
# Revised by Zhao Bin, Apr. 26, 2019. to Object Oriented Code.

import logging
import numpy as np
from scipy import linalg
from pyginvc.Greens.tde import tde
from pyginvc.libs import geotools as gt
from pyginvc.Greens.BaseGreen import BaseGreen

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")


class Meade(BaseGreen):
    '''
    Meade is a class representing triangualr dislocation
    '''
    def __init__(self, flt, data, dict_green):
        '''
        Constructor.
        
        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        '''
        super(Meade, self).__init__(flt, data, dict_green)

    def GenGreens(self, flt, data, dict_green):
        '''
        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'

        '''
        
        llh_gps       = data.llh_gps
        llh_lev       = data.llh_lev
        llh_sar       = data.llh_sar
        ndim          = data.ndim
        unit          = data.unit
        
        origin        = flt.origin
        node          = flt.vertex_enu
        element       = flt.element
        
        greentype     = dict_green['greentype']
        nu            = dict_green['nu']
        mu            = float(dict_green['modulus'])
        
        ss            = greentype[0]
        ds            = greentype[1]
        op            = greentype[2]
        
        # convert the LLH to local coordinate
        xy_gps        = np.zeros((len(llh_gps),2))
        xy_lev        = np.zeros((len(llh_lev),2))
        xy_sar        = np.zeros((len(llh_sar),2))
        if len(llh_gps) > 0:
            for i in range(len(llh_gps)):
                xy_gps[i,:] = gt.llh2localxy(llh_gps[i], origin)
            
        if len(llh_lev) > 0:
            for i in range(len(llh_lev)):
                xy_lev[i,:] = gt.llh2localxy(llh_lev[i], origin)
            
        if len(llh_sar) > 0:
            for i in range(len(llh_sar)):
                xy_sar[i,:] = gt.llh2localxy(llh_sar[i], origin)
        
        # if we have GPS data
        G_dis = np.array([])
        if len(xy_gps) > 0 and len(element)>0:
            G_dis = self.MakeGGPS_Tridisloc(node, element, xy_gps, mu, nu, ss, ds, op, ndim)
    
            # print the status
            logging.info('Green function for %d GPS stations are computed.' %(len(xy_gps)))
        G = G_dis
                    
        
        # if we have SAR data
        G_sar   = np.array([])
        if len(xy_sar) > 0:
            G_sar = self.MakeGSAR_Tridisloc(unit, node, element, xy_sar, nu, ss, ds, op)
    
            # print the status
            logging.info('Green function for %d InSAR points are computed.' %(len(xy_sar)))
        
        # print the status
        logging.info('Green function for geodetic data are computed using triangular dislocation model.')
        
        self.G     = G
        self.G_sar = G_sar
        if 'gps_ramp' in dict_green.keys() and dict_green['gps_ramp']:
            self.G_gps_ramp = self.MakeGGPSRamp(xy_gps, ndim)
        if 'sar_ramp' in dict_green.keys() and dict_green['sar_ramp']:
            sizes = data.n_sar
            G_sar_ramp = []
            for i in range(sizes):
                start = sum(sizes[:i])
                end   = sum(sizes[:i+1])
                G_sar_ramp.append(self.MakeGSARRamp(xy_sar[start:end]))
            self.G_sar_ramp = linalg.block_diag(*G_sar_ramp)
 

    def MakeGGPS_Tridisloc(self, node, element, xy, mu, nu, ss, ds, op, gdim):
        '''
        Make Green Function using triangular dislocation methods for GPS data
        Written by Zhao Bin of Institute of Seismology @ UC Berkeley Dec 17, 2016
        Modified by Zhao Bin, Institute of Seismology. May 18, 2017. 
        Modified by Zhao Bin, Institute of Seismology. Aug 10, 2017. 
            Fix Bug of calling tde.calc_tri_displacements
    
        NOTICE: input of calc_tri_displacements, downward is positive
        
        Input:
            node    = array n*3, vertices contain x, y, z/ enu in km
            element = array m*3, index triangular
            xy      = array n*2, east and north of GPS station / en
            mu      = float, 3e10
            nu      = float, Possion ratio 0.25
            ss      = int, strike slip
            ds      = int, dip slip
            op      = int, open slip
            gdim    = int, dimension of GPS data 2/3
        
        Output:
            G       = array
        '''

        # number of stations
        nsta      = len(xy)
        # number of elements
        nelem     = len(element)

        element   = element-1
        # whether we used 3D data
        if gdim == 3:
            ndata = nsta*3
        else:
            ndata = nsta*2

        G         = np.zeros((ndata,3*nelem))
        for i in range(len(element)):
            # !!! please make sure this is right !!!
            X =  node[element[i,:], 0]
            Y =  node[element[i,:], 1]
            Z =  node[element[i,:], 2]

            # strike slip
            if ss == 1:
                U = tde.calc_tri_displacements(xy[:,0], xy[:,1], np.zeros(len(xy)), X, Y, Z, nu, -1, 0, 0)
                if gdim == 3:
                    G[:,3*i] = np.vstack((U['x'], U['y'], U['z'])).T.flatten()
                else:
                    G[:,3*i] = np.vstack((U['x'], U['y'])).T.flatten()
            # dip slip
            if ds == 1:
                U = tde.calc_tri_displacements(xy[:,0], xy[:,1], np.zeros(len(xy)), X, Y, Z, nu,  0, 0, -1)
                if gdim == 3:
                    G[:,3*i+1] = np.column_stack((U['x'], U['y'], U['z'])).flatten()
                else:
                    G[:,3*i+1] = np.column_stack((U['x'], U['y'])).flatten()
#               import matplotlib.pyplot as plt
#               from matplotlib.tri import Triangulation
#               tri = Triangulation(node[:,0], node[:,1], element)
#               plt.triplot(tri)
#               plt.scatter(xy[:,0], xy[:,1])
#               plt.quiver(xy[:,0], xy[:,1], U['x'], U['y'])
#               print(U['x'])
#               plt.show()
            # open slip
            if op == 1:
                U = tde.calc_tri_displacements(xy[:,0], xy[:,1], np.zeros(len(xy)), X, Y, Z, nu,  0, 1, 0)
                if gdim == 3:
                    G[:,3*i+2] = np.vstack((U['x'], U['y'], U['z'])).T.flatten()
                else:
                    G[:,3*i+2] = np.vstack((U['x'], U['y'])).T.flatten()
    
        return G
    

    def MakeGSAR_Tridisloc(self, unit, node, element, xy, nu, ss, ds, op):
        '''
        Make Green function for InSAR data and using triangular dislocation
        Written by Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley
        Nov 17. 2016
        
        Input:
            unit    = array, n*3 vector of LOS direction
            node    = array, n*3 vertices of mesh
            element = array, n*3 indices of triangular
            xy      = array, n*2 location of SAR data
            nu      = float, Possion ratio 0.25
            ss      = int, strike slip
            ds      = int, dip slip
            op      = int, open slip
        
        Output:
            G       = array, green function
        
        '''
        # number of SAR data
        nsta  = len(xy)
        # number of elements
        nelem = len(element)
    
        # Descending track unit vector
        if len(unit) == 0:
            lookangle = 23.0
            track = 13.9
            rad_lookangle = np.deg2rad(lookangle)
            rad_track     = np.deg2rad(track)
            unit  = np.array([-np.cos(rad_track)*np.sin(rad_lookangle),
                               np.sin(rad_track)*np.sin(rad_lookangle),
                              -np.cos(rad_lookangle)])
            unit  = -1.0*unit
            
        # init Green function G
        G = np.zeros((nsta, 3*nelem))
    
        # for each element
        for i in range(nelem):
            range_1 = np.zeros((nsta,1))
            range_2 = np.zeros((nsta,1))
            range_3 = np.zeros((nsta,1))
    
            X = node[element[i,:], 0]
            Y = node[element[i,:], 1]
            Z = node[element[i,:], 2]
    
            # strike slip
            if ss == 1:
                ux, uy, uz = tde.calc_tri_displacements(xy[:,0], xy[:,1], np.zeros(len(xy)), X, Y, Z, nu,  1, 0, 0)
                u = np.array([ux, uy, uz])
                for j in range(nelem):
                    range_1[j,:] = unit[j,:].dot(u[j,:])
                    G[:,3*i] = range_1.T
            # dip slip
            if ds == 1:
                ux, uy, uz = tde.calc_tri_displacements(xy[:,0], xy[:,1], np.zeros(len(xy)), X, Y, Z, nu,  0, 1, 0)
                u = np.array([ux, uy, uz])
                for j in range(nelem):
                    range_2[j,:] = unit[j,:].dot(u[j,:])
                    G[:,3*i+1] = range_2.T
            # open slip
            if op == 1:
                ux, uy, uz = tde.calc_tri_displacements(xy[:,0], xy[:,1], np.zeros(len(xy)), X, Y, Z, nu,  0, 0, 1)
                u = np.array([ux, uy, uz])
                for j in range(nelem):
                    range_3[j,:] = unit[j,:].dot(u[j,:])
                    G[:,3*i+2] = range_3.T
    
        return G

