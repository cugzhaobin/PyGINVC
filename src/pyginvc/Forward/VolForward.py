#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA. May 15, 2021

from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Volcano import Volcano
from pyginvc.Greens.Mogi import Mogi
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

class VolForward(object):
    '''
    Output is a class representing forward surface displacement from volcano Mogi model
    '''
    
    def __init__(self, dict_volcano, dict_data, dict_green, dict_weight):
        '''
        Constructor.
        
        Input:
            dict_volcano     = an dict contains parameters about volcano
            dict_data        = an dict contains parameters about data
            dict_green       = an dict contains parameters about Greens Function
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
        # Volcano
        #
        volfile  = dict_volcano['volcanofile']
        self.vol = Volcano(volfile, origin=[])

        #
        # GreenFunction
        #
        self.green = Mogi(self.vol, self.data, dict_green)

        self.wsar  = dict_weight['wsar']
        
        
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
    
        # get volumn change of volcano
        volumn   = self.vol.volumn
    
        dhat = np.empty(0) 
        # compute the modeled observation
        if len_geod > 0:
            if len_sar > 0:
                dhat  = vstack((self.green.G, self.green.G_sar)).dot(volumn)
            else:
                dhat  = self.green.G.dot(volumn)
        else:
            if len_sar > 0:
                dhat  = self.green.G_sar.dot(volumn) 
        dhat = dhat.flatten()
        r    = D - dhat
    
        # if we have SAR data
        if len_sar > 0:
            r_sar               = r[len_geod:len_all]
            wrss_sar            = self.wsar*r_sar.dot(self.wsar*r_sar)
            # print the status
            logging.info('SAR Weighted Residual Sum of Squares (WRSS) = %f' %(wrss_sar))
            logging.info('SAR WRSS / (N) = %f ' %(wrss_sar/(len_sar)))
            
        # if we have GPS data
        if len_gps > 0:
            r_gps        = self.data.W.dot(D[0:len_gps] - dhat[0:len_gps])
            # print the status
            logging.info('GPS Weighted Residual Sum of Squares (WRSS) = %f' %(r_gps.dot(r_gps)))
            logging.info('GPS WRSS / (N) = %f ' %(r_gps.dot(r_gps)/(len_gps)))
    
        # if we have LEV data
        if len_lev > 0:
            r_lev        = self.data.W.dot(D[len_gps:len_geod] - dhat[len_gps:len_geod])
            # print the status
            logging.info('LEV Weighted Residual Sum of Squares (WRSS) = %f' %(r_lev.dot(r_lev)))
            logging.info('LEV WESS / (N) = %f ' %(r_lev.dot(r_lev)/(len_lev)))
    
        self.dhat = dhat
        self.r    = r

        return
