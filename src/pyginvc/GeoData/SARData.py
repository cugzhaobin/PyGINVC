#!/usr/bin/env python
import numpy as np
import os, sys, logging
from pyginvc.libs import geotools as gt

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
        self.llh_sar = np.array([])
        self.d_sar   = np.array([])
        self.unit    = np.array([])
        self.cov_sar = np.array([])

        self.LoadSARData(sarfile, scale=1.0)
        self.MakeCVM(sigma, lamda)
        return

    def set_wsar(self, wsar):
        '''
        '''
        self.wsar = wsar

    def LoadSARData(self, sarfile, scale=1.0):
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
        if sarfile != '':
            [llh_sar, d_sar, unit] = self.GetSARdata(sarfile)
            self.llh_sar  = llh_sar
            self.d_sar    = d_sar.flatten()*scale
            self.unit     = unit
        else:
            self.d_sar    = np.array([])
            self.llh_sar  = np.array([])
            self.unit     = np.array([])
    
        	# print processing status
        logging.info('{} InSAR points are loaded.'.format(len(self.llh_sar)))
        
        return 

    def MakeCVM(self, sigma, lamda):
        '''
        Construct variance-covariance matrix
        Written by Zhao Bin, Jul. 31, 2019.
        
        Input:
            sigma     = Array/ scalar representing amplitude
            lamda     = Array/ scalar representing scale length
        Output:
            cov_sar = covariance of InSAR data
        '''

        nlos      = len(self.llh_sar)
        if sigma == 'None' or sigma == '':
            self.cov_sar = np.eye(nlos)
        else:
            if lamda == 'None' or sigma == '':
                self.cov_sar = np.eye(nlos)*sigma**2
            else:
                distances = np.zeros((nlos, nlos))
                for i in range(nlos):
                    for j in range(i+1,nlos):
                        en             = gt.llh2localxy(self.llh_sar[j,[1,0]], self.llh_sar[i,[1,0]])
                        distances[i,j] = np.sqrt(sum(en**2))
                        distances[j,i] = distances[i,j]
                self.cov_sar   = sigma*(np.exp(-(distances**2)/(2*lamda**2)))

        Icov_sar         = np.linalg.inv(self.cov_sar)
        self.W_sar       = np.linalg.cholesky(Icov_sar)

        return

    
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

        # check the file is exit
        if os.path.isfile(filename) is False:
            logging.debug('InSAR file is not exist! Please edit YAML file!')
            sys.exit()
    
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
