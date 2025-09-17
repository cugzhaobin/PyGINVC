#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 7 2016
# Revised by Zhao Bin, Apr. 26, 2019. to Object Oriented Code.
# Revised by Zhao Bin, Dev. 12, 2021. Rotate Green's function when rake_beta !=0

import numpy as np
import logging
from pyginvc.Greens.tdcalc import TDdispHS
from pyginvc.Greens.BaseGreen import BaseGreen

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class Nikkhoo(BaseGreen):
    '''
    Nikkhoo is a class representing triangualr dislocation
    '''
    def __init__(self, flt, data, dict_green):
        '''
        Constructor.
        
        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        '''
        super(Nikkhoo, self).__init__(flt, data, dict_green)
#       import h5py, os
        
#       greenfile = dict_green['greenfile']   
#       if greenfile == "":
#           self.GenGreens(flt, data, dict_green)
#       elif greenfile == "SAVE":
#           self.GenGreens(flt, data, dict_green)
#           with h5py.File('greenfunc.h5', 'w') as h5:
#               h5.create_dataset('G', data = self.G, compression='gzip')
#               h5.create_dataset('G_sar', data = self.G_sar, compression='gzip')
#       elif os.path.isfile(greenfile):
#           with h5py.File(greenfile, 'r') as h5:
#               self.G     = h5['G'].value[()]
#               self.G_sar = h5['G_sar'].value[()]
#               logging.info('Load Greens function from {}'.format(greenfile))
#       else:
#           self.GenGreens(flt, data, dict_green)
#       self.modulus = float(dict_green['modulus'])
#       return


    def GenGreens(self, flt, data, green_dict):
        '''
        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'

        '''
        from pyginvc.libs import geotools as gt
        llh_gps       = data.llh_gps
        llh_lev       = data.llh_lev
        llh_sar       = data.llh_sar
        ndim          = data.ndim
        unit          = data.unit
        
        origin        = flt.origin
        node          = flt.vertex_enu
        element       = flt.element
        
        greentype     = green_dict['greentype']
        nu            = green_dict['nu']
        mu            = float(green_dict['modulus'])
        
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
                xy_gps[i,:] = gt.llh2utm(llh_gps[i], origin)
            
        if len(llh_lev) > 0:
            for i in range(len(llh_lev)):
                xy_lev[i,:] = gt.llh2localxy(llh_lev[i], origin)
                xy_lev[i,:] = gt.llh2utm(llh_lev[i], origin)
            
        if len(llh_sar) > 0:
            for i in range(len(llh_sar)):
                xy_sar[i,:] = gt.llh2localxy(llh_sar[i], origin)
                xy_sar[i,:] = gt.llh2utm(llh_sar[i], origin)
#       import matplotlib.pyplot as plt
#       from matplotlib.tri import Triangulation
#       tri = Triangulation(node[:,0], node[:,1], element-1)
#       plt.scatter(xy_gps[:,0], xy_gps[:,1])
#       plt.triplot(tri)
#       plt.show()
        
        logging.info('Begin to generating Greens function.')
        # if we have GPS data
        G_dis = np.array([])
        if len(xy_gps) > 0 and len(element)>0:
            logging.info('Begin to compute Greens function for GPS station.')
            G_dis = self.MakeGGPS_Tridisloc(node, element, xy_gps, mu, nu, ss, ds, op, ndim)
    
            # print the status
            logging.info('Green function for %d GPS stations are computed.' %(len(xy_gps)))
        G = G_dis
                    
        
        # if we have SAR data
        G_sar   = np.array([])
        if len(xy_sar) > 0:
            logging.info('Begin to compute Greens function for InSAR pixels.')
            G_sar = self.MakeGSAR_Tridisloc(unit, node, element, xy_sar, mu, nu, ss, ds, op)
    
            # print the status
            logging.info('Green function for %d InSAR points are computed.' %(len(xy_sar)))
        
        # print the status
        logging.info('Green function for geodetic data are computed using triangular dislocation model.')
        
        self.G     = G
        self.G_sar = G_sar        
 

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
        # whether we used 3D data
        if gdim == 3:
            ndata = nsta*3
        else:
            ndata = nsta*2

        G         = np.zeros((ndata,3*nelem))
        element   = element-1
        rec       = np.column_stack((xy, np.zeros(len(xy))))
        for i in range(len(element)):
            # !!! please make sure this is right !!!
            src = node[element[i]]
            src[:,2] = -src[:,2]

            # strike slip
            if ss == 1:
                disp = TDdispHS(rec, src, [1,0,0], nu)
                if gdim == 3:
                    G[:, 3*i] = disp.flatten()
                elif gdim == 2:
                    G[:, 3*i] = disp[:,[0,1]].flatten()
            # dip slip
            if ds == 1:
                disp = TDdispHS(rec, src, [0,1,0], nu)
                if gdim == 3:
                    G[:, 3*i+1] = disp.flatten()
                elif gdim == 2:
                    G[:, 3*i+1] = disp[:,[0,1]].flatten()

#               import matplotlib.pyplot as plt
#               from matplotlib.tri import Triangulation
#               tri = Triangulation(node[:,0], node[:,1], element)
#               plt.triplot(tri)
#               plt.scatter(xy[:,0], xy[:,1])
#               plt.quiver(xy[:,0], xy[:,1], disp[:,0], disp[:,1])
#               plt.show()
            # open slip
            if op == 1:
                disp = TDdispHS(rec, src, [0,0,1], nu)
                if gdim == 3:
                    G[:, 3*i+2] = disp.flatten()
                elif gdim == 2:
                    G[:, 3*i+2] = disp[:,[0,1]].flatten()
    
        return G
    

    def MakeGSAR_Tridisloc(self, unit, node, element, xy, mu, nu, ss, ds, op):
        '''
        Make Green function for InSAR data and using triangular dislocation
        Written by Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley
        Nov 17. 2016

        2024-05-09: Fix a bug nelement--> nsta
        
        Input:
            unit    = array, n*3 vector of LOS direction
            node    = array, n*3 vertices of mesh
            element = array, n*3 indices of triangular
            xy      = array, n*2 location of SAR data
            mu      = float, 3e10
            nu      = float, Possion ratio 0.25
            ss      = int, strike slip
            ds      = int, dip slip
            op      = int, open slip
        
        Output:
            G       = array, green function
        
        '''
        
        # import libs
        from numpy import sin, cos, deg2rad
    
        # number of SAR data
        nsta  = len(xy)
        # number of elements
        nelem = len(element)
    
        # Descending track unit vector
        if len(unit) == 0:
            lookangle = 23.0
            track = 13.9
            unit  = np.array([-cos(deg2rad(track))*sin(deg2rad(lookangle)),
                               sin(deg2rad(track))*sin(deg2rad(lookangle)),
                               -cos(deg2rad(lookangle))])
            unit  = -1.0*unit
            
        # init Green function G
        G         = np.zeros((nsta, 3*nelem))
        element   = element-1
        rec       = np.column_stack((xy, np.zeros(len(xy))))
    
        # for each element
        for i in range(nelem):
            range_1 = np.zeros((nsta,1))
            range_2 = np.zeros((nsta,1))
            range_3 = np.zeros((nsta,1))
    
            src      = node[element[i]]
            src[:,2] = -src[:,2]
    
            # strike slip
            if ss == 1:
                disp = TDdispHS(rec, src, [1,0,0], nu)
                for j in range(nsta):
                    range_1[j] = unit[j].dot(disp[j])
                G[:,3*i] = range_1.T
            # dip slip
            if ds == 1:
                disp = TDdispHS(rec, src, [0,1,0], nu)
                for j in range(nsta):
                    range_2[j] = unit[j].dot(disp[j])
                    G[:,3*i+1] = range_2.T
            # open slip
            if op == 1:
                disp = TDdispHS(rec, src, [0,0,1], nu)
                for j in range(nsta):
                    range_3[j] = unit[j].dot(disp[j])
                    G[:,3*i+2] = range_3.T
    
        return G

