#!/usr/bin/env python
from numpy import zeros
import numpy as np
import logging
from pyginvc.GeoData.GPSData import GPSData
from pyginvc.GeoData.SARData import SARData
from pyginvc.GeoData.LEVData import LEVData
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")

class GeoData(GPSData, LEVData, SARData):
    '''
    GeoData representing Geodetic Data used for slip inversion.
    '''
    llh_gps = []
    llh_lev = []
    llh_sar = []
    d_gps   = []
    d_lev   = []
    d_sar   = []
    unit    = []
    W       = []
    ndim    = 3

    def __init__(self, gpsfile, sarfile, levfile, gfiletype):
        '''
        Constructor.

        Load All Data for inversion, including GPS, InSAR and Level data
        Written by Zhao Bin @ UC Berkeley April 2016

        Input:
            gpsfile     = GPS file name
            sarfile     = SAR file name
            levfile     = Level file name
            gfiletype   = GPS file type, GMT2D or GMT3D or IOS3D
        '''
        self.llh_gps = np.array([])
        self.llh_sar = np.array([])
        self.llh_lev = np.array([])
        self.d_gps   = np.array([])
        self.d_sar   = np.array([])
        self.d_lev   = np.array([])
        self.unit    = np.array([])
        self.LoadAllData(gpsfile, sarfile, levfile, gfiletype)
        return

    def LoadAllData(self, gpsfile, sarfile, levfile, gfiletype):
        '''
        Load all geodetic data from input files

        Input:
            gpsfile     = GPS file name
            sarfile     = SAR file name
            levfile     = Level file name
            gfiletype   = GPS file type, GMT2D or GMT3D or IOS3D
        Output:

        '''

        # GPS Data
        self.LoadGPSData(gpsfile, gfiletype)        

        # LEV Data
        self.LoadLEVData(levfile)

        # SAR Data
        self.LoadSARData(sarfile)
       
        # concatenate data, covariance, and weight matricies
        # only combine GPS and level data, excluding InSAR data
        up     = np.hstack((self.W_gps, zeros((self.d_gps.size, self.d_lev.size))))
        dw     = np.hstack((zeros((self.d_lev.size, self.d_gps.size)), self.W_lev))
        self.W = np.vstack((up, dw))
        
        # print processing status
        logging.info('Finished load geodetic data for inversion.')
        return 

    def DumpData(self):
        '''
        '''
        import h5py

        with h5py.File('data.h5', 'w') as h5:
            h5['d_gps'] = self.d_gps
            h5['d_lev'] = self.d_lev
            h5['d_sar'] = self.d_sar
            h5['W']     = self.W
            h5['llh_gps'] = self.llh_gps
            h5['llh_sar'] = self.llh_sar
            h5['llh_lev'] = self.llh_lev
            h5['W_gps']   = self.W_gps
            h5['ndim']    = self.ndim
            h5['unit']    = self.unit
        return
