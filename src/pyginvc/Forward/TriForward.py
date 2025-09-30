#!/usr/bin/env python
# Written by Zhao Bin, at UC Berkeley. April 22 2016
# MOD by Zhao Bin, at UC Berkeley. Dec 4, 2016. Adding greenfile as input parameter for run_forward
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
# Mod by Zhao Bin, Dec. 7, 2018. We use HDF5 to store green functions and slotions

from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Triangle import Triangle
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

class TriForward(object):
    '''
    Output is a class representing forward surface displacement from slip model
    '''
    
    def __init__(self, dict_fault, dict_data, dict_green, dict_weight):
        '''
        Constructor.
        
        Input:
            flt         = an instance of class Fault
            data        = an instance of class GeoData
            green       = an instance of class GreenFunction
            dict_weight = a dict
        '''
        
        #
        # GeoData
        #
        gpsfile = dict_data['gpsfile']
        sarfile = dict_data['sarfile']
        levfile = dict_data['levfile']
        gfiletype = dict_data['gfiletype']
        self.data = GeoData(gpsfile, sarfile, levfile, gfiletype)
    
        #
        # Fault
        #
        vertexfile  = dict_fault['vertexfile']
        elementfile = dict_fault['elementfile']
        self.flt  = Triangle(vertexfile, elementfile, origin=[])

        #
        # GreenFunction
        #
        if dict_green['grnmethod'] == 'meade':
            from pyginvc.Greens.Meade import Meade
            self.green = Meade(self.flt, self.data, dict_green)
        elif dict_green['grnmethod'] == 'nikkhoo':
            from pyginvc.Greens.Nikkhoo import Nikkhoo
            self.green = Nikkhoo(self.flt, self.data, dict_green)
        elif dict_green['grnmethod'] == 'poly3d':
            from pyginvc.Greens.TriPoly3D import TriPoly3D
            self.green = TriPoly3D(self.flt, self.data, dict_green)
        
        
        #
        # Weight
        self.wsar  = np.array(dict_weight['wsar'])
        #
        
        #
        # Do Forward
        #
        self.run_forward()



    def run_forward(self):
        '''
        Run forward.
        '''
    
        # import libs
        from numpy import hstack, vstack
        import datetime
        
        #
        # Data
        #
        len_gps  = len(self.data.d_gps)
        len_lev  = len(self.data.d_lev)
        len_sar  = len(self.data.d_sar)
        len_geod = len_gps+len_lev
        len_all  = len_geod+len_sar
        d_gps    = self.data.d_gps.reshape(len_gps)
        d_lev    = self.data.d_lev.reshape(len_lev)
        d_sar    = self.data.d_sar.reshape(len_sar)
        D        = hstack((d_gps, d_lev, d_sar))
    
        # get fault slip
        slip     = self.flt.slip.flatten()
    
        # compute the Moment
        shearmodulus= self.green.modulus
        [Mo, Mw]  = self.moment_tri_element(self.flt.vertex_enu, self.flt.element, slip, shearmodulus)
    
        # print the status
        now       = datetime.datetime.now()
        now       = now.strftime('%y%m%d:%H%M:%S')
        logging.info('Geodetic Moment Magnitude M0 = %E' %(Mo))
        logging.info('Geodetic Moment Magnitude Mw = %f' %(Mw))
    
        dhat = np.empty(0) 
        # compute the modeled observation
        if len_geod > 0:
            if len_sar > 0:
                dhat  = vstack((self.green.G, self.green.G_sar[:,0:len(slip)])).dot(slip)
            else:
                dhat  = self.green.G.dot(slip)
        else:
            if len_sar > 0:
                dhat  = self.green.G_sar[:,0:len(slip)].dot(slip) 
        r         = D - dhat
    
        # if we have SAR data
        if len_sar > 0:
            r_sar               = r[len_geod:len_all]
            wrss_sar            = self.wsar*r_sar.dot(self.wsar*r_sar)
            # print the status
            now = datetime.datetime.now()
            now = now.strftime('%y%m%d:%H%M:%S')
            logging.info('SAR Weighted Residual Sum of Squares (WRSS) = %f' %(wrss_sar))
            logging.info('SAR WRSS / (N) = %f ' %(wrss_sar/(len_sar)))
            
        # if we have GPS data
        if len_gps > 0:
            r_gps        = self.data.W.dot(D[0:len_gps] - dhat[0:len_gps])
            # print the status
            now = datetime.datetime.now()
            now = now.strftime('%y%m%d:%H%M:%S')
            logging.info('GPS Weighted Residual Sum of Squares (WRSS) = %f' %(r_gps.dot(r_gps)))
            logging.info('GPS WRSS / (N) = %f ' %(r_gps.dot(r_gps)/(len_gps)))
    
        # if we have LEV data
        if len_lev > 0:
            r_lev        = self.data.W.dot(D[len_gps:len_geod] - dhat[len_gps:len_geod])
            # print the status
            now = datetime.datetime.now()
            now = now.strftime('%y%m%d:%H%M:%S')
            logging.info('LEV Weighted Residual Sum of Squares (WRSS) = %f' %(r_lev.dot(r_lev)))
            logging.info('LEV WESS / (N) = %f ' %(r_lev.dot(r_lev)/(len_lev)))
    
        self.dhat = dhat
        self.r    = r
        self.slip = slip.reshape(len(slip),1)
        self.misfit = np.empty(0)
        # output the modeled observation
#       out = Output(self.flt, self.data, self.green)
#     out.WriteData(dhat, r)

        
    @staticmethod
    def moment_tri_element(vertex, element, slip, shearmodulus):
        '''
        Calculate moment and magnitude for triangular element
        Written by Zhao Bin, Institute of Seismology, CEA. May 2017
        Mod by Zhao Bin, Apr. 24, 2019. Now 6.067 is used.
    
        Input:
            vertex   : coordinate for vertex in ENU
            element  : index of vertices for each triangular
            slip     : list/array slip value for all elememts
            shearmodulus:
        Output:
            Mo       : total moment
            Mw       : magnitude
        '''
        import numpy as np
        
        nelem   = len(element)
        Mo      = np.zeros(nelem)
        
        for i in range(nelem):
            id1 = element[i,0]-1
            id2 = element[i,1]-1
            id3 = element[i,2]-1
            
            len1 = np.sqrt((vertex[id1,0]-vertex[id2,0])**2 + 
                           (vertex[id1,1]-vertex[id2,1])**2 +
                           (vertex[id1,2]-vertex[id2,2])**2)
            len2 = np.sqrt((vertex[id1,0]-vertex[id3,0])**2 + 
                           (vertex[id1,1]-vertex[id3,1])**2 +
                           (vertex[id1,2]-vertex[id3,2])**2)
            len3 = np.sqrt((vertex[id3,0]-vertex[id2,0])**2 + 
                           (vertex[id3,1]-vertex[id2,1])**2 +
                           (vertex[id3,2]-vertex[id2,2])**2)
            s     = 0.5*(len1+len2+len3)
            area  = np.sqrt(s*(s-len1)*(s-len2)*(s-len3))
            Mo[i] = 1e6 * area * np.abs(np.sqrt(slip[3*i]**2 + slip[3*i+1]**2)) * shearmodulus
        Mo_total  = np.sum(Mo)/1000.0
        Mw_total  = 2.0/3.0*np.log10(Mo_total) - 6.067
        
        return Mo_total, Mw_total        
