#!/usr/bin/env python
import numpy as np
import os, logging
from scipy.spatial.distance import cdist
from pyproj import Proj

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")

class SARData(object):
    '''
    SARData representing InSAR data used for slip inversion.
    '''

    def __init__(self, sarfile, scale=1.0, sigma=None, lamda=None):
        '''
        Constructor.

        Written by Zhao Bin @ UC Berkeley April 2016
        Mod by Zhao Bin, adding parameters sig and lam for construct covariance

        Input:
            sarfile     = SAR file name
            scale       = scale LOS value
            sigma       = 
            lamda       =
        '''
        self.sarfile = sarfile
        self.scale   = scale
        self.sigma   = sigma
        self.lamda   = lamda

        self.llh_sar = np.empty((0,2))
        self.d_sar   = np.empty(0)
        self.unit    = np.empty((0,3))
        self.cov_sar = np.empty(0)

        #self.LoadSARData(sarfile, scale=1.0)
        #self.MakeCVM(sigma, lamda)
        #return

    def load(self):
        self.load_sar_data(self.scale)
        self.build_spatial_cov()

    def set_wsar(self, wsar):
        '''
        '''
        self.wsar = wsar

    def load_sar_data(self, scale=1.0):
        '''
        Load InSAR data from input files

        Input:
            sarfile     = SAR file name
            scale       = scale LOS value. e.g. from meter to mm
        Output:

        '''

        #
        # Read SAR Data
        #
        filename = self.sarfile
        if not filename:
            # Empty sarfile is a placeholder for "this track has no data".
            return
        if not os.path.exists(filename):
            raise FileNotFoundError(
                'SAR file %s does not exist.' % filename)
        [llh_sar, d_sar, unit] = self.GetSARdata(filename)
        self.llh_sar  = llh_sar
        self.d_sar    = d_sar.flatten()*scale
        self.unit     = unit
    
        # print processing status
        logging.info('{} InSAR points are loaded.'.format(len(self.llh_sar)))

    def build_spatial_cov(self):
        '''
        Construct variance-covariance matrix
        Written by Zhao Bin, Jul. 31, 2019.
        
        Input:
            sigma     = Array/ scalar representing amplitude
            lamda     = Array/ scalar representing scale length
        Output:
            cov_sar = covariance of InSAR data
        '''
        sigma  = self.sigma
        lamda  = self.lamda
        nlos   = len(self.llh_sar)
        if nlos == 0:
            # No InSAR points: average distance / projection below are
            # meaningless (nan) and would raise inside pyproj.
            self.cov_sar = np.empty(0)
        elif sigma is None or lamda is None:
            self.cov_sar = np.ones(nlos)
        else:
            origin  = np.mean(self.llh_sar, axis=0)
            P       = Proj(ellps='WGS84', proj='tmerc', lon_0=origin[1])
            en      = P(self.llh_sar[:, 1], self.llh_sar[:, 0])
            en_org  = P(origin[1], origin[0])
            localxy = np.column_stack([(en[0] - en_org[0]) / 1000.0,
                                       (en[1] - en_org[1]) / 1000.0])
            dist    = cdist(localxy, localxy, metric="euclidean")
            cov_sar = sigma**2 * np.exp(-dist**2/(2*lamda**2))
            self.cov_sar = cov_sar
    
    @staticmethod
    def GetSARdata(filename):
        '''
        Read coseismic InSAR LOS data
        Written by Zhao Bin when his is in UC Berkeley, Jan 28 2016
    
        Input:
            filename  = InSAR LOS file name, the unit is mm
                      - Lon Lat LOS E N U
        Output:
            llh       = latitute longitute height
            ran       = vector of rnage change
            unit      = unit vector
            cov_sar   = covariance matrix of LOS
        '''    
        # load the data
        data = np.genfromtxt(filename, usecols=[0,1,2,3,4,5], comments='#')  
    
        # print read info
        nsta = len(data)
    
        # init some variables
        ran  = np.zeros(nsta)
        llh  = np.zeros((nsta, 2))
        unit = np.zeros((nsta, 3))
    
        # extract data
        llh    = data[:,[1,0]]
        ran    = data[:,2]
        unit   = data[:,3:6]
    
        # print processing status
        logging.info('{} SAR points are read in.'.format(nsta))
        return llh, ran, unit