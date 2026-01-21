#!/usr/bin/env python
from numpy import zeros
import numpy as np
import logging
from scipy import linalg
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

    def __init__(self, gpsfiles, sarfiles, levfiles, gfiletype):
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
        self.LoadGeoData(gpsfiles, sarfiles, levfiles, gfiletype)
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

    def LoadGeoData(self, gpsfiles, sarfiles, levfiles, gfiletype):
        
        gpsfiles = self._to_list(gpsfiles)
        sarfiles = self._to_list(sarfiles)
        levfiles = self._to_list(levfiles)

        # 
        gps_data, sar_data, lev_data = [], [], []
        for gpsfile in gpsfiles:
            try:
                gps_data.append(GPSData(gpsfile, gfiletype))
            except Exception as e:
                logging.info(f"Loading GPS file {gpsfile}: {e}")
        for sarfile in sarfiles:
            try:
                sar_data.append(SARData(sarfile))
            except Exception as e:
                logging.info(f"Loading SAR file {sarfile}: {e}")
        for levfile in levfiles:
            try:
                lev_data.append(LEVData(levfile))
            except Exception as e:
                logging.info(f"Loading LEV file {levfile}: {e}")

        # Now combine the data together
        llh_gps = [item.llh_gps for item in gps_data]
        llh_sar = [item.llh_sar for item in sar_data]
        llh_lev = [item.llh_lev for item in lev_data]
        d_gps   = [item.d_gps for item in gps_data]
        d_sar   = [item.d_sar for item in sar_data]
        d_lev   = [item.d_lev for item in lev_data]
        unit    = [item.unit  for item in sar_data]
        w_gps   = [item.W_gps for item in gps_data]
        w_sar   = [item.W_sar for item in sar_data]
        w_lev   = [item.W_lev for item in lev_data]
        cov_gps = [item.cov_gps for item in gps_data]
        n_gps   = [len(item.llh_gps) for item in gps_data]
        n_sar   = [len(item.llh_sar) for item in sar_data]
        sta_gps = [item.station_gps for item in gps_data]
        w       = w_gps+w_lev

        self.llh_gps = np.vstack(llh_gps)
        self.llh_sar = np.vstack(llh_sar)
        self.llh_lev = np.vstack(llh_lev)
        self.d_gps   = np.hstack(d_gps)
        self.cov_gps = linalg.block_diag(*cov_gps)
        self.d_sar   = np.hstack(d_sar)
        self.d_lev   = np.hstack(d_lev)
        self.unit    = np.vstack(unit)
        self.W       = linalg.block_diag(*w)
        self.W_gps   = w_gps
        self.W_sar   = linalg.block_diag(*w_sar)
        self.n_gps   = np.array(n_gps)
        self.n_sar   = np.array(n_sar)
        self.station_gps = np.hstack(sta_gps)

        if gfiletype[-2:] == '2D':
            self.ndim = 2
        if gfiletype[-2:] == '3D':
            self.ndim = 3

        # print processing status
        logging.info('Finished load geodetic data for inversion.')
        
    def _to_list(self, file_input):
        if isinstance(file_input, str):
            return [file_input]
        elif isinstance(file_input, list):
            return file_input
        
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
