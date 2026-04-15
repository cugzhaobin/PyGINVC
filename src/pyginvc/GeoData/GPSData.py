#!/usr/bin/env python
import numpy as np
import os, sys, logging

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")

class GPSData(object):
    '''
    GPSData representing GPS Data used for slip inversion.
    '''
    def __init__(self, gpsfile, gfiletype):
        '''
        Constructor.

        Load All Data for inversion, including GPS, InSAR and Level data
        Written by Zhao Bin @ UC Berkeley April 2016

        Input:
            gpsfile     = GPS file name
            gfiletype   = GPS file type, GMT2D or GMT3D or IOS3D
        '''
        self.gpsfile     = gpsfile
        self.gpsfiletype = gfiletype

        self.llh_gps     = np.empty((0,2))
        self.d_gps       = np.array([])
        self.cov_gps     = np.array([])
        self.Icov_gps    = np.array([])
        self.W_gps       = np.array([])
        self.station_gps = np.array([], dtype=object)
        self.ndim        = 2  # Default to 3D data

    def load(self):
        self.load_gps_data()
        self.build_weight()

    def load_gps_data(self):
        '''
        Load GPS data.
        Input:
            gpsfile     = GPS file name
            gfiletype   = GPS file type, GMT2D or GMT3D or IOS3D
        '''
        gpsfile   = self.gpsfile
        gfiletype = self.gpsfiletype
        if not os.path.exists(gpsfile):
            logging.fatal(f"{gpsfile} file does not exit.")

        if gfiletype in ('GMT2D', 'IOS2D'):
            vel, llh_gps, cov_gps, station_gps = self.GetGMTdata(gpsfile)
        
        elif gfiletype == 'GMT3D':
            vel, llh_gps, cov_gps, station_gps = self.GetGMT3Ddata(gpsfile);
        
        elif gfiletype == 'IOS3D':
            vel, llh_gps, cov_gps, station_gps = self.GetIOS3Ddata(gpsfile);
        else:
            logging.fatal(f'Unknown GPS file type: {gfiletype}')

        self.nsta        = len(llh_gps)
        self.ndim        = vel.shape[1]
        self.llh_gps     = llh_gps
        self.d_gps       = vel.flatten()
        self.cov_gps     = cov_gps
        self.station_gps = station_gps
        #self.Icov_gps    = np.linalg.inv(cov_gps)
        #self.W_gps       = np.linalg.cholesky(self.Icov_gps)
        logging.info(f'{self.nsta} GPS stations in {gfiletype} format are loaded.')

    def build_weight(self):
        self.W_gps = 1/self.cov_gps

    @staticmethod
    def GetGMTdata(filename, scale=1.0):
        '''   
         Read coseismic data in GMT format
         Written by Zhao Bin when he was in UC Berkeley, Jan 28 2016
         Mod by Zhao Bin, Mar. 27, 2019.
    
         Input:
            filename  = GPS file name in GMT format, the unit is mm
         Output:
            vel       = vector of velocities/offsets: east, north
            cov       = covariance matrix of velocities
            llh       = latitute longitute height: lat, lon
            station   = station name
        
        '''
        # check the file is exit
        if not os.path.isfile(filename):
            logging.critical('{} does not exist! Please edit YAML file!'.format(filename))
    
        try:
            # load the data
            data = np.genfromtxt(filename, comments='#')
            if data.ndim == 1:
                data = data.reshape((1, -1))
    
            # load stations
            if data.shape[1]>=7:
                station = np.genfromtxt(filename, usecols=7, comments='#', dtype=None, encoding='utf-8')
            else:
                station = np.array([])

            # extract data and rescale the data, unit is mm
            llh      = data[:,[1,0]]
            vel      = data[:,2:4]
            esig     = data[:,4]*scale
            nsig     = data[:,5]*scale
            esig[esig == 0] = 1.0
            nsig[nsig == 0] = 1.0
            
            # make cov
            nsta     = len(data)
        #   cov      = np.zeros((2*nsta, 2*nsta), dtype='float32')
            cov      = np.zeros(2*nsta)
    
        #   for i in range(nsta):
        #       k = 2*i
        #       cov[k,k]     = esig[i]**2
        #       cov[k+1,k+1] = nsig[i]**2
            cov[0::2] = esig**2
            cov[1::2] = nsig**2
            return vel, llh, cov, station
        except Exception as e:
            logging.error(f'Error reading GMT 2D data: {str(e)}')
            raise ValueError(f'Invalid GMT 2D data format in {filename}') from e
    
    @staticmethod
    def GetGMT3Ddata(filename, scale=1.0):
        '''   
         Read coseismic data in GMT format
         Written by Zhao Bin @  UC Berkeley, April 15 2016
    
         Input:
            filename  = GPS file name in GMT format
                        Lon Lat Ve    Vn    Vu    Se    Sn    Su  Site
                        deg deg mm    mm    mm    mm    mm    mm
         Output:
            vel       = vector of velocities
            cov       = covariance matrix of velocities
            llh       = latitute longitute height
            station   = station name
        
        '''    
        # check the file is exit
        if not os.path.isfile(filename):
            logging.critical('{} does not exist! Please edit YAML file!'.format(filename))
            sys.exit()
    
        try:
            # load the data
            data = np.genfromtxt(filename, comments='#')  
            if data.ndim == 1: data = data.reshape((1, -1))
        
            # load stations  
            station = np.genfromtxt(filename, usecols=8, comments='#', dtype=None, encoding='utf-8') if data.shape[1] == 9 else np.array([])
        
            # print read info
            nsta = len(data)
        
            # extract data and rescale the data
            llh      = data[:,[1,0]]
            vel      = data[:,2:5]
            esig     = data[:,5]*scale
            nsig     = data[:,6]*scale
            usig     = data[:,7]*scale
            esig[esig == 0] = 1.0
            nsig[nsig == 0] = 1.0
            usig[usig == 0] = 1.0
        
            # make cov
            cov = np.zeros((3*nsta, 3*nsta))
        
            for i in range(nsta):
                k = 3*i
                cov[k,k]     = esig[i]**2
                cov[k+1,k+1] = nsig[i]**2
                cov[k+2,k+2] = usig[i]**2

            return vel, llh, cov, station
        except Exception as e:
            logging.error(f'Error reading GMT 3D data: {str(e)}')
            # print processing status
           

    @staticmethod
    def GetIOS3Ddata(filename, scale=1.0):
        '''   
         Read coseismic data in extended GMT format defined by Institute of Seismology
         Written by Zhao Bin, Mar. 27, 2019.
    
         Input:
            filename  = GPS file name in GMT format
                        Lon Lat Ve    Vn    Se    Sn    Cne   Site  Vu   Su
                        deg deg mm    mm    mm    mm    mm          mm   mm
         Output:
            vel       = vector of velocities
            cov       = covariance matrix of velocities
            llh       = latitute longitute height
            station   = station name
        
        '''
    
        # check the file is exit
        if not os.path.isfile(filename):
            logging.critical('{} does not exist! Please edit YAML file!'.format(filename))
            sys.exit()

        try:
            # load the data
            data = np.genfromtxt(filename, usecols=[0,1,2,3,8,4,5,9], comments='#')  
            if data.ndim == 1: data = data.reshape((1, -1))
        
            # load stations
            station = np.genfromtxt(filename, usecols=7, comments='#', dtype=None, encoding='utf-8')
        
            # print read info
            nsta = len(data)
            if nsta < 2:
                logging.fatal('Number of GNSS station < 2!')
                sys.exit()

            # extract data and rescale the data
            llh      = data[:,[1,0]]
            vel      = data[:,2:5]
            esig     = data[:,5]*scale
            nsig     = data[:,6]*scale
            usig     = data[:,7]*scale
            esig[esig == 0] = 1.0
            nsig[nsig == 0] = 1.0
            usig[usig == 0] = 1.0
            
            # make cov
            cov = np.zeros([3*nsta, 3*nsta])
        
            for i in range(nsta):
                k = 3*i
    
                cov[k,k]     = esig[i]**2
                cov[k+1,k+1] = nsig[i]**2
                cov[k+2,k+2] = usig[i]**2

            # return the result
            return vel, llh, cov, station
        except Exception as e:
            logging.error(f'Error reading IOS3D data: {str(e)}')
            raise ValueError(f'Invalid IOS3D data format in {filename}') from e
            
