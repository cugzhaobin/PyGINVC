#!/usr/bin/env python
# Written by Zhao Bin, May 15, 2021

import numpy as np
import datetime
import logging
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class Mogi(object):
    '''
    Okada is a class representing rectangualr dislocation
    '''    
    
    def __init__(self, vol, data, dict_green):
        '''
        Constructor.
        
        Input:
            vol        = an instance of class Volcano
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        '''
        
        import h5py, os
        self.vol = vol
        greenfile = dict_green['greenfile']
        if greenfile == "":
            self.GenGreens(vol, data, dict_green)
        elif greenfile == "SAVE":
            self.GenGreens(vol, data, dict_green)
            with h5py.File('greenfunc.h5', 'w') as h5:
                h5.create_dataset('G', data = self.G, compression='gzip')
                h5.create_dataset('G_sar', data = self.G_sar, compression='gzip')
        elif os.path.isfile(greenfile):
            with h5py.File(greenfile, 'r') as h5:
                self.G     = h5['G'][()]
                self.G_sar = h5['G_sar'][()]
        else:
            self.GenGreens(vol, data, dict_green)
        return


    def GenGreens(self, vol, data, dict_green):
        '''
        Generate Green's Function.

        Input:
            vol        = an instance of class Volcano
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        Output:
            self.G     = Green's function for GPS
            self.G_sar = Green's function for InSAR
        '''
        import numpy as np
        from pyginvc.libs import geotools as gt

        llh_gps       = data.llh_gps
        llh_lev       = data.llh_lev
        llh_sar       = data.llh_sar
        ndim          = data.ndim
        unit          = data.unit
        origin        = vol.origin
        nu            = dict_green['nu']
        if 'verbose' not in dict_green.keys():
            verbose   = True
        else:
            verbose   = dict_green['verbose']
    
        # convert the LLH to local coordinate
        xy_gps        = np.zeros((len(llh_gps),2))
        xy_lev        = np.zeros((len(llh_lev),2))
        xy_sar        = np.zeros((len(llh_sar),2))
        if len(llh_gps) > 0:
            for i in range(len(llh_gps)):
                xy_gps[i,:] = gt.llh2localxy(llh_gps[i], origin)
            
        if len(llh_lev) > 0:
            for i in range(len(llh_lev)):
                xy_lev[i,:] = gt.llh2localxy(llh_lev[i], origin)
            
        if len(llh_sar) > 0:
            for i in range(len(llh_sar)):
                xy_sar[i,:] = gt.llh2localxy(llh_sar[i], origin)
    
        # if we have GPS data
        G_dis = np.array([])
        if len(xy_gps) > 0 and len(vol.location_llh)>0:
            G_dis = self.MakeGGPS(xy_gps, nu, ndim)
    
            # print the status
            if verbose:
                logging.info('Green function for %d GPS stations are computed.' %(len(xy_gps)))
        G = G_dis
            
        # if we have level data
        G_dis = np.array([])    
        if len(xy_lev) > 0:
            G_dis = self.MakeGLEV(xy_lev, nu)
            G = np.hstack((G, G_dis))
    
            # print the status
            now = datetime.datetime.now()
            now = now.strftime('%y%m%d:%H%M:%S')
            if verbose:
                logging.info('Green function for %d Level stations are computed.' %(len(xy_lev)))
        
        # if we have SAR data
        G_sar   = np.array([])
        if len(xy_sar) > 0:
            G_sar = self.MakeGSAR(unit, xy_sar, nu)
    
            # print the status
            if verbose:
                logging.info('Green function for %d InSAR points are computed.' %(len(xy_sar)))
        
        self.G     = G
        self.G_sar = G_sar

        # print the status
        if verbose:
            logging.info('Green function for geodetic data are computed using Okada model.')
        
        return
        
    def MakeGGPS(self, xy, nu, gdim):
        '''
        Make green's function for GPS data
        Written by Zhao Bin when his is in UC Berkeley, Jan 28 2016
        Subroutine to make design matrix for slip
    
        Input:
            xy        = xy coordinates of GPS stations
            nu        = Poisson's ratio
            gdim      = 2/3
        Output:
            G         = Green's Function for GPS
     
        '''
        nvol = self.vol.nvol
        nsta = xy.shape[0]
    
        # whether we used 3D data
        if gdim == 3:
            ndata = nsta*3
        else:
            ndata = nsta*2
    
        # init for G
        G = np.zeros((ndata, nvol))
    
        # for each subfaults
        for i in range(nvol):
            volgeom = self.vol.location_enu[i]
            u       = self.relative_disp(nu, volgeom, xy)
            if gdim == 3:
                G[:,i] = u.flatten()
            elif gdim == 2:
                G[:,i] = u[:,[0,1]].flatten()

        return G
        
        
    def MakeGSAR(self, unit, xy, nu):
        '''
        Make Green's Function for InSAR data
        This treats range change data as absolute with respect to a reference
        level of zero deformation.  In practice the reference level is uncertain
        and the inversion should allow for a best-fitting reference level.
        Hence the appearance of 3*nf+1 parameters below.
    
        Input:
            unit         -  unit vector of InSAR data
            xy           -  numpy array, xy coordinates of stations
            nu           -  Poisson's ratio
    
        Output:
            G            -  numpy array, desigin matrix
        
        '''
        from numpy import sin, cos, deg2rad
        nvol = self.vol.nvol
        nsta = xy.shape[0]
    
        # Descending track unit vector
        if len(unit) == 0:
            lookangle = 23
            track     = 13.9
            unit      = np.array([-cos(deg2rad(track))*sin(deg2rad(lookangle)), 
                                  sin(deg2rad(track))*sin(deg2rad(lookangle)), 
                                 -cos(deg2rad(lookangle))])
            unit      = -1.0*unit
        G = np.zeros((nsta, nvol))
    
        los = np.zeros(nsta)
        # for each volcano
        for i in range(nvol):
            volgeom = self.vol.location_enu[i]
            u       = self.relative_disp(nu, volgeom, xy)
            for j in range(nsta):
                los[j] = unit[j,:].dot(u[j,:])
            G[:,i] = los
        return G

    
    def relative_disp(self, nu, volgeom, xy):
        '''
        Compute relative green function following DISMODEL. In fact we do use all data 
        
        Input:
            nu       =  Poisson's ratio
            volgeom  =  numpy array with 3 column and 1 record
            xy       =  numpy array with 2 column, east and north in kilometer
    
        Output:
            u        =  green function of all stations in xy for a record of subfault
        '''
    
        e = xy[:,0] - volgeom[0]
        n = xy[:,1] - volgeom[1]
        r = np.sqrt(e**2, n**2)
        d = volgeom[2]

        prefactor = (1-nu)*1e3/(np.pi*d**2)
        uz        = prefactor*(1+r**2/d**2)**(-1.5)
        ur        = (prefactor/d)*r*(1+r**2/d**2)**(-1.5)

        # convert to radial displacement to east north
        theta     = np.arctan2(n, e)
        ue        = ur*np.cos(theta)
        un        = ur*np.sin(theta)
        u         = np.column_stack((ue, un, uz))
        return u
