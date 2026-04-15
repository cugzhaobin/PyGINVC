#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 7 2016
# Revised by Zhao Bin, Apr. 26, 2019. to Object Oriented Code.

import logging
import numpy as np
from scipy import linalg
from numpy import cos, sin, deg2rad
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

    def build_greens(self):
        if not self.load_greens():
            self.generate_greens()
            self.save_greens()
        self.rotate_greens()
        logging.info("Building Green's functions finished.")

    def generate_greens(self):
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
        data          = self.data
        flt           = self.flt
        dict_green    = self.dict_green

        llh_gps       = data.llh_gps
        llh_lev       = data.llh_lev
        llh_sar       = data.llh_sar
        ndim          = data.ndim
        unit          = data.unit
        origin        = flt.origin
        dis_geom_grid = flt.dis_geom_grid
        
        ss, ds, op    = dict_green['greentype']
        nu            = dict_green['nu']
        verbose       = dict_green.get('verbose', True)
        
        # convert the LLH to local coordinate
        xy_gps        = self._convert_to_local(llh_gps, origin)
        xy_sar        = self._convert_to_local(llh_sar, origin)
        xy_lev        = self._convert_to_local(llh_lev, origin)

        """
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
        """
    
        # if we have GPS data
        if xy_gps.size > 0 and dis_geom_grid.size > 0:
            self.G = self.make_gps_green(dis_geom_grid, xy_gps, nu, ss, ds, op, ndim)
            if verbose:
                logging.info(f"Green's functions for {len(xy_gps)} GPS stations are computed.")
            
        # if we have level data
        if xy_lev.size > 0:
            G_lev  = self.make_gps_green(dis_geom_grid, xy_lev, nu, ss, ds, op)
            self.G = np.vstack((self.G, G_lev))
            if verbose:
                logging.info(f"Green's functions for {len(xy_lev)} level stations are computed.")
        
        # if we have SAR data
        if xy_sar.size > 0:
            self.G_sar = self.make_sar_green(unit, dis_geom_grid, xy_sar, nu, ss, ds, op)
            if verbose:
                logging.info(f"Green's functions for {len(xy_sar)} InSAR points are computed.")
        
        if dict_green.get('gps_ramp', False):
            self.G_gps_ramp = self.make_green_gps_ramp(xy_gps, ndim, method=1)
        else:
            self.G_gps_ramp = np.empty((self.G.shape[0],0))
        if dict_green.get('sar_ramp', False):
            sizes = data.n_sar
            ramp  = []
            for i in range(len(sizes)):
                start = sum(sizes[:i])
                end   = sum(sizes[:i+1])
                ramp.append(self.make_green_sar_ramp(xy_sar[start:end]))
            self.G_sar_ramp = linalg.block_diag(*ramp)
        else:
            self.G_sar_ramp = np.empty((self.G_sar.shape[0],0))

        logging.info("Green's function are computed using Okada model.")
        
    def make_gps_green(self, dis_geom, xy, nu, ss, ds, op, gdim):
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
        if dis_geom.ndim == 1: 
            dis_geom = dis_geom.reshape(-1,10)
        if xy.ndim == 1: 
            xy = xy.reshape(-1,2)

        nf    = self.flt.nf
        nobs  = len(xy) * self.data.ndim    
        G     = np.zeros((nobs, 3*nf))
    
        # for each subfaults
        for i in range(nf):
            # if we want to constrain strike slip 
            if bool(ss) is True:
                source      = np.hstack((dis_geom[i,0:7], [1,0,0]))
                u           = self.relative_disp(nu, source, xy)
                G[:, 3*i+0] = u[:gdim].T.ravel()
         
            # if we want to constrain dip slip 
            if bool(ds) is True:
                source      = np.hstack((dis_geom[i,0:7], [0,1,0]))
                u           = self.relative_disp(nu, source, xy) 
                G[:, 3*i+1] = u[:gdim].T.ravel()
    
            # if we want to constrain openning slip 
            if bool(op) is True:
                source      = np.hstack((dis_geom[i,0:7], [0,0,1]))
                u           = self.relative_disp(nu, source, xy)
                G[:, 3*i+2] = u[:gdim].T.ravel()

        return G
        
        
    def make_sar_green(self, unit, dis_geom, xy, nu, ss, ds, op):
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
            xy = xy.reshape(-1,2)
            
        nf     = self.flt.nf
        ndata  = self.data.n_sar

        # Descending track unit vector
        if len(unit) == 0:
            lookangle = 23
            track     = 13.9
            unit      = np.array([-cos(deg2rad(track))*sin(deg2rad(lookangle)), 
                                  sin(deg2rad(track))*sin(deg2rad(lookangle)), 
                                 -cos(deg2rad(lookangle))])
            unit      = -1.0 * unit
        G = np.zeros((ndata, 3*nf))
    
        # for each subfaults
        for i in range(nf):   
            # SS motion.
            if bool(ss):
                u          = self.disloc(nu, np.hstack((dis_geom[i,0:7],[1,0,0])), xy)
                G[:,3*i+0] = (unit * u).sum(axis=1)
    
    	    # DS motion.
            if bool(ds): 
                u          = self.disloc(nu, np.hstack((dis_geom[i,0:7],[0,1,0])), xy)
                G[:,3*i+1] = (unit * u).sum(axis=1)
    
    	    # Opening
            if bool(op): 
                u          = self.disloc(nu, np.hstack((dis_geom[i,0:7],[0,0,1])), xy)
                G[:,3*i+2] = (unit * u).sum(axis=1)
 
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
        nsta = xy.shape[0]
        u    = np.zeros((nsta, 3))
        for i in range(nsta):
            u[i,:] = self.disloc(nu, dis_geom, np.array([xy[i,0], xy[i,1]]))
    
        # compute  displacements realtive to LAST station
        u_rel      = np.zeros((3,nsta))
        u_rel[0,:] = u[:,0].T
        u_rel[1,:] = u[:,1].T
        u_rel[2,:] = u[:,2].T
        u_rel      = u.T
    
        return u_rel
