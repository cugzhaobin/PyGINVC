#!/usr/bin/env python
# Written by Zhao Bin, May 15, 2021

import numpy as np
import sys, os, logging
import pyginvc.libs.geotools as gt
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class Volcano(object):
    '''
    Triangle is a class representing volcano
    '''

    def __init__(self, volfile, origin=[]):
    
        '''
        Load volcano location and volumn changes
        Written by Zhao Bin, Institute of Seismology, CEA. May 15, 2021
        
        Input:
            volfile : file containins volcano location and volumn changes
        
        '''
        # check the file is exit
        isfile = os.path.isfile(volfile)
        if (isfile == False):
           logging.warning ('{} is not exist! Please input another file!'.format(vertexfile))
           sys.exit()
     
        # read in location of volcanos
        location = np.genfromtxt(volfile, comments='#', usecols=[0,1,2])
        if location.ndim == 1: location = location.reshape((1,len(location)))
        
        # read in volumn changes
        volumn  = np.genfromtxt(volfile, comments='#', usecols = [3])
        
        # calculate the origin
        origin = np.mean(location, axis=0)[0:2]
        
        # convert geo coordinate into local coordinate
        location_enu = np.zeros((len(location),3))
        for i in range(len(location)):
            en = gt.llh2localxy(location[i,0:2], origin)
            location_enu[i,:] = np.array([en[0], en[1], location[i,2]])
        
        # return results
        self.location_llh = location
        self.location_enu = location_enu
        self.origin     = origin
        self.nvol       = len(location)
        self.volumn     = volumn

