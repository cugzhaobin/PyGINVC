#!/usr/bin/env python2
# Written by Zhao Bin, at UC Berkeley. April 22 2016
# MOD by Zhao Bin, at UC Berkeley. Dec 4, 2016. Adding greenfile as input parameter for run_forward
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
# Mod by Zhao Bin, Dec. 7, 2018. We use HDF5 to store green functions and slotions
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Fault import Fault
from pyginvc.Greens.Okada import Okada
from pyginvc.Inversion.GeoInversion import GeoInversion
from pyginvc.Export.Output import Output
import logging
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")


class Forward(object):
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
        faultfile = dict_fault['faultfile']
        nsegs     = dict_fault['nsegs']
        ndeps     = dict_fault['ndeps']
        doSubFault= dict_fault['doSubFault']
        origin    = dict_fault['origin']
        self.flt  = Fault(faultfile, nsegs, ndeps, doSubFault, origin=origin)

        #
        # GreenFunction
        #
        self.green = Okada(self.flt, self.data, dict_green)
        
        #
        # Weight
        self.wsar      = dict_weight['wsar']
        self.data.wsar = dict_weight['wsar']
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
        slip      = self.flt.dis_geom_grid[:,7:10].flatten()
    
        # compute the Moment
        shearmodel= 3E10
        [Mo, Mw]  = GeoInversion.Moment(self.flt.dis_geom_grid, slip, shearmodel) 
    
        # print the status
        logging.info('Geodetic Moment Magnitude M0 = %E' %(Mo))
        logging.info('Geodetic Moment Magnitude Mw = %f' %(Mw))
    
    
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
            wrss_sar            = (self.wsar*r_sar).dot(self.wsar*r_sar)
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
            logging.info('LEV WRSS / (N) = %f ' %(r_lev.dot(r_lev)/(len_lev)))
    
        # output the modeled observation
        out = Output(self.flt, self.data, self.green)
        out.WriteData(dhat=dhat, r=r)
#       out.WriteData()
