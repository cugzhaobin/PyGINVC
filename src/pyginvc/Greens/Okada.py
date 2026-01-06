#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 7 2016
# Revised by Zhao Bin, Apr. 26, 2019. to Object Oriented Code.

import logging
import numpy as np
from numpy import cos, sin, deg2rad, hstack
from pyginvc.libs import geotools as gt
from pyginvc.Greens import okada85 as okada
from pyginvc.Greens.BaseGreen import BaseGreen
#from okada_wrapper import dc3dwrapper
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class Okada(BaseGreen):
    '''
    Okada is a class representing rectangualr dislocation
    '''    
    
    def __init__(self, flt, data, dict_green):
        '''
        Constructor.
        
        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        '''
        super(Okada, self).__init__(flt, data, dict_green)
        return


    def GenGreens(self, flt, data, dict_green):
        '''
        Generate Green's Function.

        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        Output:
            self.G     = Green's function for GPS
            self.G_sar = Green's function for InSAR
        '''

        llh_gps       = data.llh_gps
        llh_lev       = data.llh_lev
        llh_sar       = data.llh_sar
        ndim          = data.ndim
        unit          = data.unit
        origin        = flt.origin
        dis_geom_grid = flt.dis_geom_grid
        
        greentype     = dict_green['greentype']
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
        if len(xy_gps) > 0 and len(dis_geom_grid)>0:
            G_dis = self.MakeGGPS(dis_geom_grid, xy_gps, nu, greentype[0], greentype[1], greentype[2], ndim)
    
            # print the status
            if verbose:
                logging.info('Green function for %d GPS stations are computed.' %(len(xy_gps)))
        G = G_dis
            
        # if we have level data
        G_dis = np.array([])    
        if len(xy_lev) > 0:
            G_dis = self.MakeGLEV(dis_geom_grid, xy_lev, nu, greentype[0], greentype[1], greentype[2])
            G = np.hstack((G, G_dis))
    
            if verbose:
                logging.info('Green function for %d Level stations are computed.' %(len(xy_lev)))
        
        # if we have SAR data
        G_sar   = np.array([])
        if len(xy_sar) > 0:
            G_sar = self.MakeGSAR(unit, dis_geom_grid, xy_sar, nu, greentype[0], greentype[1], greentype[2])
    
            # print the status
            if verbose:
                logging.info('Green function for %d InSAR points are computed.' %(len(xy_sar)))
        
        self.G     = G
        self.G_sar = G_sar
        if 'gps_ramp' in dict_green.keys():
            self.G_gps_ramp = self.MakeGGPSRamp(xy_gps, ndim)
        if 'sar_ramp' in dict_greens.keys():
            self.G_sar_ramp = self.MakeGSARRamp(xy_sar, ndim)

        # print the status
        if verbose:
            logging.info('Green function for geodetic data are computed using Okada model.')
        
        return
        
    def MakeGGPS(self, dis_geom, xy, nu, ss, ds, op, gdim):
        '''
        Make green's function for GPS data
        Written by Zhao Bin when his is in UC Berkeley, Jan 28 2016
        Subroutine to make design matrix for slip
    
        Input:
            dis_geom  = Dislocation Geometry
            xy        = xy coordinates of GPS stations
            nu        = Poisson's ratio
            ss        = 0/1
            ds        = 0/1
            op        = 0/1
            gdim      = 2/3
        Output:
            G         = Green's Function for GPS
     
        '''
        if dis_geom.ndim == 1: dis_geom = dis_geom.reshape(len(dis_geom),10)
        if xy.ndim == 1: xy = xy.reshape(len(xy),2)
        nf   = dis_geom.shape[0]
        nsta = xy.shape[0]
    
        # whether we used 3D data
        if gdim == 3:
            ndata = nsta*3
        else:
            ndata = nsta*2
    
        # init for G
        G = np.zeros((ndata,3*nf))
    
        # for each subfaults
        for i in range(nf):
            # if we want to constrain strike slip 
            if bool(ss) is True:
                source         = np.hstack((dis_geom[i,0:7], [1,0,0]))
                u_rel_1        = self.relative_disp(nu, source, xy) 
                if gdim == 3:
                    G[:,3*i]   = u_rel_1.T.flatten()
                else:
                    u_rel_1    = u_rel_1[0:2,:]
                    G[:,3*i]   = u_rel_1.T.flatten()
         
            # if we want to constrain dip slip 
            if bool(ds) is True:
                source         = np.hstack((dis_geom[i,0:7], [0,1,0]))
                u_rel_2        = self.relative_disp(nu, source.flatten(), xy) 
                if gdim == 3:
                    G[:,3*i+1] = u_rel_2.T.flatten()
                else:
                    u_rel_2    = u_rel_2[0:2,:]
                    G[:,3*i+1] = u_rel_2.T.flatten()
    
            # if we want to constrain openning slip 
            if bool(op) is True:
                source         = np.hstack((dis_geom[i,0:7], [0,0,1]))
                u_rel_3        = self.relative_disp(nu, source.flatten(), xy) 
                if gdim == 3:
                    G[:,3*i+2] = u_rel_3.T.flatten()
                else:
                    u_rel_3    = u_rel_3[0:2,:]
                    G[:,3*i+2] = u_rel_3.T.flatten()
    
        return G
        
        
    def MakeGSAR(self, unit, dis_geom, xy, nu, ss, ds, op):
        '''
        Make Green's Function for InSAR data
        This treats range change data as absolute with respect to a reference
        level of zero deformation.  In practice the reference level is uncertain
        and the inversion should allow for a best-fitting reference level.
        Hence the appearance of 3*nf+1 parameters below.
    
        Input:
            unit         -  unit vector of InSAR data
            dis_geom     -  numpy array, Dislocation Geometry
            xy           -  numpy array, xy coordinates of stations
            nu           -  Poisson's ratio
            ss           -  if we will estimate strike   slip, ss=1, else ss=0
            ds           -  if we will estimate dip      slip, ds=1, else ds=0
            op           -  if we will estimate openning slip, op=1, else op=0
    
        Output:
            G            -  numpy array, desigin matrix
        
        '''
    
    
        # 
        if dis_geom.ndim == 1: 
            dis_geom = dis_geom.reshape(len(dis_geom),10)
        if xy.ndim == 1: 
            xy = xy.reshape(len(xy),2)
            
        nf     = dis_geom.shape[0]
        ndata  = xy.shape[0]

        # Descending track unit vector
        if len(unit) == 0:
            lookangle = 23
            track     = 13.9
            unit      = np.array([-cos(deg2rad(track))*sin(deg2rad(lookangle)), 
                                  sin(deg2rad(track))*sin(deg2rad(lookangle)), 
                                 -cos(deg2rad(lookangle))])
            unit      = -1.0*unit
        G = np.zeros((ndata, 3*nf+3))
    
        # for each subfaults
        for i in range(nf):   
            range_1 = np.zeros((ndata,1))
            range_2 = np.zeros((ndata,1))
            range_3 = np.zeros((ndata,1))
            # SS motion.
            if ss != 0:
                u        = self.disloc(nu, hstack((dis_geom[i,0:7],[1,0,0])), xy)
                for j in range(ndata):
                    range_1[j,:] = unit[j,:].dot(u[j,:])
                G[:,3*i] = range_1.T
    
    	    # DS motion.
            if ds != 0: 
                u          = self.disloc(nu, hstack((dis_geom[i,0:7],[0,1,0])), xy)
                for j in range(ndata):
                    range_2[j,:] = unit[j,:].dot(u[j,:])
                G[:,3*i+1] = range_2.T
    
    	    # Opening
            if op != 0: 
                u          = self.disloc(nu, hstack((dis_geom[i,0:7],[0,0,1])), xy)
                for j in range(ndata):
                    range_3[j,:] = unit[j,:].dot(u[j,:])
                G[:,3*i+2] = range_3.T
    
        # The following is for adding a constant and slope to the SAR data
#       for j in range(ndata):
#           G[j,3*nf+0] = 1.0
#           G[j,3*nf+1] = xy[i,0]
#           G[j,3*nf+2] = xy[i,1]
     
        return G

    @staticmethod
    def disloc(nu, disgeom, x):
        '''
        Subroutine to make design matrix for slip

        Input:
            nu         = Poisson's ratio
    
            Dislocation Geometry in a Global E-N-Up coordinate system
            disgeom[0] = len = fault length in strike direction (km)
            disgeom[1] = wid = fault width in dip direction (km)
            disgeom[2] = dep = depth  of lower edge of fault (km)
            disgeom[3] = dip = dip angle (degrees)
            disgeom[4] = strik =  strike, clockwise from N (degrees)
            disgeom[5] = delE  = East offset of midpoint from origin (km)
            disgeom[6] = delN  = North offset of midpoint from origin (km)
            disgeom[7] = ss  =  strike slip motion (m)
            disgeom[8] = ds  =  dip slip motion (m)
            disgeom[9] = op  =  opening  motion (m)
    
            Vector of receiver coordinates at which displacement are computed
            x  = vector [East (km), North (km)]
        Output: Displacements in the Global coordinate system
    
        '''
    
        alp     = 1 - 2*nu
        dip     = disgeom[3]
        if dip == 90 or dip == 270:
           cs   = 0.0
           sn   = 1.0
        else:
           dipr = np.deg2rad(dip)
           cs   = np.cos(dipr)
           sn   = np.sin(dipr)

        strkr   = np.deg2rad(90 - disgeom[4])
    
        # get the size of xy
        x       = np.mat(x)
        [m,n]   = x.shape
        x       = np.array(x)
    
        delE    = disgeom[5]
        delN    = disgeom[6]
    
        # init u
        u       = np.zeros((m,3))
        ut      = np.zeros(9)
    
        for k in range(m):
            # transform the global receiver coordinates into local coordinates
            zlocal = (complex(x[k,0],x[k,1]) - complex(delE, delN)) * np.exp(complex(0,-strkr)) + 0.5*disgeom[0]
            xrec   = np.array([zlocal.real, zlocal.imag])
     
            # compute green function
#           [success, ut, grad_u] = dc3dwrapper(alp, [xrec[0],xrec[1],0], disgeom[2], dip, 
#                                               [0, disgeom[0]], [0, disgeom[1]], [disgeom[7], disgeom[8], disgeom[9]])
            okada.okadaio(alp, [xrec[0],xrec[1]], [disgeom[2], 0, disgeom[0], 0, disgeom[1], sn, cs, disgeom[7], disgeom[8], disgeom[9]], ut)
            
            # transformaiton back to global coordinates
            zglob  = complex(ut[0], ut[1])*np.exp(complex(0, strkr))
            u[k,0] = zglob.real
            u[k,1] = zglob.imag
            u[k,2] = ut[2]
    
    
        return u
    
    
    def relative_disp(self, nu, dis_geom, xy):
        '''
        Compute relative green function following DISMODEL. In fact we do use all data 
        
        Input:
            nu       =  Poisson's ratio
            dis_geom =  numpy array with 10 column and 1 record
            xy       =  numpy array with 2 column, east and north in kilometer
    
        Output:
            u_rel    =  green function of all stations in xy for a record of subfault
        '''
    
        # get the number of stations in xy
        xy         = np.mat(xy)
        [nsta, m]  = xy.shape
        xy         = np.array(xy)
      
        # init variables
        u          = np.zeros((nsta, 3))
        for i in range(nsta):
            u[i,:] = self.disloc(nu, dis_geom, np.array([xy[i,0], xy[i,1]]))
    
        # compute  displacements realtive to LAST station
        u_rel      = np.zeros((3,nsta))
        u_rel[0,:] = u[:,0].T
        u_rel[1,:] = u[:,1].T
        u_rel[2,:] = u[:,2].T
        u_rel      = u.T
    
        return u_rel
