#!/usr/bin/env python
import numpy as np
import os
import logging
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")

class LEVData(object):
    '''
    LEVData representing Leveling Data used for slip inversion.
    '''
    llh_lev = []
    d_lev   = []

    def __init__(self, levfile):
        '''
        Constructor.

        Written by Zhao Bin @ UC Berkeley April 2016

        Input:
            levfile     = Level file name
        '''
        self.LoadLEVData(levfile)
        return


    def LoadLEVData(self, levfile):
        '''
        Load leveling data from input files

        Input:
            levfile     = Level file name
        Output:

        '''
        from numpy import zeros
        #
        # Read Level Data
        #
        if levfile != '':
            logging.info('Loading level data now')
    #       [d_lev, llh_lev, Len, sigma_lev, Icov_lev, W_lev] = GetLEVdata(levfile)
            self.llh_lev     = zeros((0,0))
            self.cov_lev     = np.eye(0)
        else:
            self.llh_lev     = zeros((0,2))
            self.d_lev       = zeros(0)
            self.Icov_lev    = zeros((0,0))
            self.W_lev       = zeros((0,0))
    
    	    # print processing status
            logging.info('%d leveling stations are loaded.' %(len(self.llh_lev)))


    @staticmethod
    def GetLEVdata(filename):
        '''   
         Read coseismic data in GMT format
         Written by Zhao Bin @ UC Berkeley, April 15 2016
    
         Input:
            filename  = LEV file name, Data Format is
                        lon, lat, height, length of line, elevation change, sigma
                        
         Output:
            d         = vector of elevation changes or rates
            llh       = longitute latitute height
            Len       = distance from reference benchmark to benchmark
            sigma     = vector of standard deviations: sqrt(diag(cov))
            Icov      = inverse of covariance matrix
            W         = weight matrix       
        '''
    
        # check the file is exit
        if os.path.isfile(filename) is False:
            logging.critical(filename+' is not exist! Please input another file!')
#           sys.exit()
    
        # load the data
        data = np.genfromtxt(filename, usecols=[0,1,2,3,4,5], comments='#')  
    
        # print read info
        nsta = len(data)
    
        # init some variables
        llh = np.zeros((nsta,3))
        Len = np.zeros((nsta,1))
    
        # extract data
        llh = data[0:nsta, 1:3]
        Len = data[1:nsta, 3]
        d   = data[0:nsta, 4]
        s   = data[0:nsta, 5]
        
        # Form weight matrix
        W   = np.zeros((len(d)-1, len(d)-1))
        W[0,0] = np.sqrt(1/s[0]**2/Len[0])
    
        # Allows for multiple line segments
        for i in range(1, len(d)-1):
            if Len[i] > Len[i-1]:
                W[i,i]    = np.sqrt(1/(s[i]**2*Len[i]-s[i-1]**2*Len[i-1]))
                W[i,i-1]  = -1*W[i,i]
            else:
                W[i,i]    = np.sqrt(1/(s[i]**2*2))
                
                
        # Form inverse of the covariance matrix
        Icov   = np.tranpose(W)*W
        rcov0  = np.invert(Icov)
        rcov1  = np.array([])
        
        for i in range(0, len(d)-1):
            rcov1[i,0] = 0.0
            rcov1[0,i] = 0.0
            
            for j in range(0, len(d)-1):
                rcov1[i+1,j+1] = rcov1[i,j]
                
                
        # Use root-mean-square error in absolute vertical heights of 1.0mm/yr
        sigv  =  0.000158
        
        for i in range(0, len(d)):
            rcov1[i,0] = rcov1[i,i] + sigv**2
            
        Icov  = np.invert(rcov1)
        sigma = np.sqrt(np.diag(rcov1))
        Len   = np.hstack([0, Len])        
        W     = np.linglg.cholesky(np.mat(Icov))
    
    
        # print the status
        logging.info('%d Leveling stations are read in.' %(llh))
        
        return d, llh, Len, sigma, Icov, W
