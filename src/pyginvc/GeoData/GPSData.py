#!/usr/bin/env python
from numpy import array, genfromtxt, zeros
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
    llh_gps = []
    d_gps   = []
    ndim    = 3

    def __init__(self, gpsfile, gfiletype):
        '''
        Constructor.

        Load All Data for inversion, including GPS, InSAR and Level data
        Written by Zhao Bin @ UC Berkeley April 2016

        Input:
            gpsfile     = GPS file name
            gfiletype   = GPS file type, GMT2D or GMT3D or IOS3D
        '''
        self.llh_gps     = np.array([])
        self.d_gps       = np.array([])
        self.ndim        = 3  # Default to 3D data
        self.cov_gps     = np.array([])
        self.Icov_gps    = np.array([])
        self.W_gps       = np.array([])
        self.station_gps = np.array([])
        
        self.LoadGPSData(gpsfile, gfiletype)
        return

    def LoadGPSData(self, gpsfile, gfiletype):
        '''
        Load GPS data.
        Input:
            gpsfile     = GPS file name
            gfiletype   = GPS file type, GMT2D or GMT3D or IOS3D
        '''

        if gpsfile != '':
            if gfiletype == 'GMT2D' or gfiletype == 'IOS2D':
                self.ndim = 2
                [vel, llh_gps, cov_gps, station_gps] = self.GetGMTdata(gpsfile)
    
    	        # print processing status
                logging.info('{} GPS stations in {} format are loaded.'.format(len(llh_gps),gfiletype))
            
            if gfiletype == 'GMT3D':
                self.ndim = 3
                [vel, llh_gps, cov_gps, station_gps] = self.GetGMT3Ddata(gpsfile);
    
    	        # print processing status
                logging.info('{} GPS stations in {} format are loaded.'.format(len(llh_gps),gfiletype))
            
            if gfiletype == 'IOS3D':
                self.ndim = 3
                [vel, llh_gps, cov_gps, station_gps] = self.GetIOS3Ddata(gpsfile);
    
    	        # print processing status
                logging.info('{} GPS stations in {} format are loaded.'.format(len(llh_gps),gfiletype))
            
            self.llh_gps     = llh_gps
            self.d_gps       = vel.flatten()
            self.cov_gps     = cov_gps
            self.Icov_gps    = np.linalg.inv(cov_gps)
            self.W_gps       = np.linalg.cholesky(self.Icov_gps)
            self.station_gps = station_gps
        else:
            self.llh_gps     = zeros((0,2))
            self.d_gps       = zeros(0)
            self.cov_gps     = zeros((0,0))
            self.Icov_gps    = zeros((0,0))
            self.W_gps       = zeros((0,0))
            self.station_gps = zeros(0)



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
            sys.exit()
    
        try:
            # load the data
            data = genfromtxt(filename, comments='#')
            if data.ndim == 1:
                data = data.reshape((1, len(data)))
    
            # load stations
            station = np.genfromtxt(filename, usecols=7, comments='#', dtype='S') if data.shape[1]==10 else array([])

            nsta = len(data)
            nvec = nsta
    
            # init some paramters
            llh = zeros((nsta,2))
            vel = zeros((nsta,2))
    
            # extract data and rescale the data, unit is mm
            llh      = data[:,[1,0]]
            vel      = data[:,2:4]*scale
            esig     = data[:,4]*scale
            nsig     = data[:,5]*scale
            esig[esig == 0] = 1.0
            nsig[nsig == 0] = 1.0
            
            # make cov
            cov      = zeros([2*nvec, 2*nvec], dtype='float32')
    
            for i in range(nsta):
                k = 2*i
                cov[k,k]     = esig[i]**2
                cov[k+1,k+1] = nsig[i]**2
                cov[k,k+1]   = 0.0
                cov[k+1,k]   = cov[k,k+1]
    
            # print the status
            logging.info('{} GPS stations are read in.'.format(nsta))
    
            # return the result
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
            data = genfromtxt(filename, comments='#')  
            if data.ndim == 1: data = data.reshape((1, len(data)))
        
            # load stations  
            station = np.genfromtxt(filename, usecols=8, comments='#', dtype='S') if data.shape[1] == 9 else array([])
        
            # print read info
            nsta = len(data)
        
            # init some variables
            vel = zeros((nsta,3))
            llh = zeros((nsta,2))
        
            # extract data and rescale the data
            llh      = data[:,[1,0]]
            vel      = data[:,2:5]*scale
            esig     = data[:,5]*scale
            nsig     = data[:,6]*scale
            usig     = data[:,7]*scale
            esig[esig == 0] = 1.0
            nsig[nsig == 0] = 1.0
            usig[usig == 0] = 1.0
        
            # make cov
            cov = zeros([3*nsta, 3*nsta])
        
            for i in range(nsta):
                k = 3*i
        
                cov[k,k]     = esig[i]**2
                cov[k+1,k+1] = nsig[i]**2
                cov[k+2,k+2] = usig[i]**2
                cov[k,k+1] = 0.0
                cov[k,k+2] = 0.0
                cov[k+1,k] = 0.0
                cov[k+2,k] = 0.0
                logging.info(f'{nsta} GPS stations are read in GMT 3D format.')
                # return the result
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
    
        from numpy import array, genfromtxt, zeros
        import os, sys
    
        # check the file is exit
        if not os.path.isfile(filename):
            logging.critical('{} does not exist! Please edit YAML file!'.format(filename))
            sys.exit()

        try:
            # load the data
            data = genfromtxt(filename, usecols=[0,1,2,3,8,4,5,9], comments='#')  
            if data.ndim == 1: data = data.reshape((1, len(data)))
        
            # load stations
            station = np.genfromtxt(filename, usecols=7, comments='#', dtype='S') if data.shape[1]==10 else array([])
        
            # print read info
            nsta = len(data)
            if nsta < 2:
                logging.fatal('Number of GNSS station < 2!')
                sys.exit()
        
            # init some variables
            llh = zeros((nsta,2))
            vel = zeros((nsta,3))
            
            # extract data and rescale the data
            llh      = data[:,[1,0]]
            vel      = data[:,2:5]*scale
            esig     = data[:,5]*scale
            nsig     = data[:,6]*scale
            usig     = data[:,7]*scale
            esig[esig == 0] = 1.0
            nsig[nsig == 0] = 1.0
            usig[usig == 0] = 1.0
            
            # make cov
            cov = zeros([3*nsta, 3*nsta])
        
            for i in range(nsta):
                k = 3*i
    
                cov[k,k]     = esig[i]**2
                cov[k+1,k+1] = nsig[i]**2
                cov[k+2,k+2] = usig[i]**2
                cov[k,k+1] = 0.0
                cov[k,k+2] = 0.0
                cov[k+1,k] = 0.0
                cov[k+2,k] = 0.0
        
            # print processing status
            logging.info('%d GPS stations are read in.' %(nsta))
         
            # return the result
            return vel, llh, cov, station
        except Exception as e:
            logging.error(f'Error reading IOS3D data: {str(e)}')
            raise ValueError(f'Invalid IOS3D data format in {filename}') from e
            
    def plotGPSVector(self, showfig=False):
        import matplotlib.pyplot as plt
        
        if len(self.llh) == 0:
            return
        
        dispE = self.d_gps[0::self.ndim]
        dispN = self.d_gps[1::self.ndim]
        plt.quiver(self.llh_gps[:,0], self.llh_gps[:,1], dispE, dispN)
        if showfig:
            plt.show()
