#!/usr/bin/env python
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

class GeoData:
    '''
    GeoData representing Geodetic Data used for slip inversion.
    '''

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
        self.gpsfiles = gpsfiles
        self.sarfiles = sarfiles
        self.levfiles = levfiles
        self.gfiletype= gfiletype

        self.llh_gps = np.array([])
        self.llh_sar = np.array([])
        self.llh_lev = np.array([])
        self.d_gps   = np.array([])
        self.d_sar   = np.array([])
        self.d_lev   = np.array([])
        self.unit    = np.array([])
       # self.W       = np.empty(0)
       # self.W_gps   = np.empty(0)
       # self.W_sar   = np.empty(0)
       # self.W_lev   = np.empty(0)

        self.cov_gps = np.empty(0)
        self.cov_sar = np.empty(0)
        self.ndim    = 2
        self.n_gps   = 0
        self.n_sar   = 0
        self.station_gps = np.empty(0, dtype=str)
    
    def load_data(self):
        self.load_geo_data()
        logging.info("All geodetic data loaded")  
    

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
        #up     = np.hstack((self.W_gps, zeros((self.d_gps.size, self.d_lev.size))))
        #dw     = np.hstack((zeros((self.d_lev.size, self.d_gps.size)), self.W_lev))
        #self.W = np.vstack((up, dw))
        
        # print processing status
        logging.info('Finished load geodetic data for inversion.')

    def load_geo_data(self):
        
        gpsfiles = self._to_list(self.gpsfiles)
        sarfiles = self._to_list(self.sarfiles)
        levfiles = self._to_list(self.levfiles)
        gfiletype = self.gfiletype

        gps_data, sar_data, lev_data = [], [], []
        for gpsfile in gpsfiles:
            try:
                gps = GPSData(gpsfile, gfiletype)
                gps.load()
                gps_data.append(gps)
            except Exception as e:
                logging.info(f"Loading GPS file {gpsfile}: {e}")
        for sarfile in sarfiles:
            try:
                sar = SARData(sarfile)
                sar.load()
                sar_data.append(sar)
            except Exception as e:
                logging.info(f"Loading SAR file {sarfile}: {e}")
        for levfile in levfiles:
            try:
                lev_data.append(LEVData(levfile))
            except Exception as e:
                logging.info(f"Loading LEV file {levfile}: {e}")

        # Now combine the data together
        if gps_data:
            self.llh_gps = np.vstack([item.llh_gps for item in gps_data])
            self.d_gps   = np.hstack([item.d_gps for item in gps_data])
            self.cov_gps = np.hstack([item.cov_gps for item in gps_data])
            self.n_gps   = np.array([len(item.llh_gps) for item in gps_data])
            self.station_gps = np.hstack([item.station_gps for item in gps_data])

        if sar_data:
            self.llh_sar = np.vstack([d.llh_sar for d in sar_data])
            self.d_sar   = np.hstack([d.d_sar for d in sar_data])
            self.unit    = np.vstack([d.unit  for d in sar_data])
            self.cov_sar = np.hstack([d.cov_sar for d in sar_data])
            self.n_sar   = np.array([len(d.llh_sar) for d in sar_data])

        if lev_data:
            self.llh_lev = np.vstack([d.llh_lev for d in lev_data])
            self.d_lev   = np.hstack([d.d_lev for d in lev_data])
        
        if gfiletype[-2:] == '2D':
            self.ndim = 2
        if gfiletype[-2:] == '3D':
            self.ndim = 3

        # print processing status
        logging.info('Finished load geodetic data for inversion.')
        
    def _to_list(self, file_input):
        if file_input is None:
            return []
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
            #h5['W']     = self.W
            h5['llh_gps'] = self.llh_gps
            h5['llh_sar'] = self.llh_sar
            h5['llh_lev'] = self.llh_lev
            #h5['W_gps']   = self.W_gps
            h5['ndim']    = self.ndim
            h5['unit']    = self.unit
            h5['cov_gps'] = self.cov_gps
            h5['cov_sar'] = self.cov_sar
        return

    def process_data(self, gps_blk_idx, sar_blk_idx):
        if len(gps_blk_idx) == 0:
            logging.info('Block index for GPS is null')
            return
        else:
            obs, llh, cov = [], [], []
            for i in np.unique(gps_blk_idx):
                if i < 0: continue
                idx = (gps_blk_idx==i)
                obs.append(self.d_gps.reshape(-1,self.ndim)[idx])
                llh.append(self.llh_gps[idx])
                cov.append(self.cov_gps.reshape(-1,self.ndim)[idx])

            self.llh_gps  = np.vstack(llh)
            self.d_gps    = np.vstack(obs).flatten()
            self.cov_gps  = np.vstack(cov).flatten()

        if len(sar_blk_idx) == 0:
            logging.info('Block index for SAR is null')
            return
        else:
            llh, los, unit, cov = [], [], [], []
            for i in np.unique(sar_blk_idx):
                if i < 0: continue
                idx = (sar_blk_idx==i)
                llh.append(self.llh_sar[idx])
                los.append(self.los[idx])
                unit.append(self.unit[idx])
                cov.append(self.cov_sar[idx][:,idx])

            self.llh_sar = np.vstack(llh)
            self.d_sar   = np.vstack(los)
            self.unit    = np.vstack(unit)
            self.cov_sar = linalg.block_diag(*cov)
    
    def get_sigma(self):
        if self.cov_gps.ndim == 1:
            return np.sqrt(self.cov_gps).reshape(-1, self.ndim)
        elif self.cov_gps.ndim == 2:
            return np.sqrt(np.diag(self.cov_gps)).reshape(-1, self.ndim)