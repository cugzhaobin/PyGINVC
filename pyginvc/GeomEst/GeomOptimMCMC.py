#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Determine the optimal fault geometry using Bayesian inversion methodology.

Created on July. 10, 2019.

@author: Zhao Bin, Institute of Seismology, CEA.
"""

from pyginvc.Geometry.Fault import Fault
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Greens.Okada import Okada
import numpy as np


class GeomOptimization:
    '''
    '''
    
    def __init__(self, dict_fault, dict_data, dict_green, dict_weight, dict_param):
        '''
        Construction the class.
        Input:
            dict_fault = dictionary type contains fault parameters
            dict_data  = dictionary type contains data parameters
            dict_green = dictionary type contains Green's Function parameters
            dict_weight= dictionary type contains weight parameters
            dict_param = dictionary type contains other parameters
        '''

        #
        # GeoData
        #
        gpsfile         = dict_data['gpsfile']
        sarfile         = dict_data['sarfile']
        levfile         = dict_data['levfile']
        gfiletype       = dict_data['gfiletype']
        self.data       = GeoData(gpsfile, sarfile, levfile, gfiletype)

        # 
        # Fault
        #
        faultfile       = dict_fault['faultfile']
        flt             = Fault(faultfile, 1, 1, False, origin=[])
        fltboundfile    = dict_fault['faultbound']
        x               = flt.geom2x(flt.geom_grid, flt.origin)
        xlb, xub        = flt.GetFaultBounds(fltboundfile, 1)
        self.x          = x
        self.xlb        = xlb
        self.xub        = xub
        self.flt        = flt

        #
        # Green function parameter
        #
        self.dict_green = dict_green

        #
        # data weighting
        #
        self.wgps       = dict_weight['wgps']
        self.wsar       = dict_weight['wsar']
        sig             = dict_weight['sig']
        lam             = dict_weight['lam']
        self.data.MakeVCM(sig=sig, lam=lam)

        #
        # bayesian
        #
        self.nsteps     = dict_param['nsteps']
        self.nburn      = dict_param['nburn']
        self.ss         = dict_param['ss']
        self.ds         = dict_param['ds']

    def getdata(self):
        '''
        '''
        return np.hstack((self.d_gps, self.d_lev, self.d_sar))


    def forward(self, ileng, iwid, idep, idip, istrk, ilat, ilon, iss, ids):
        '''
        '''
        
        x = np.array([[ileng, iwid, idep, idip, istrk, ilat, ilon]])
        [geom_short, dis_geom_short] = Fault.x2geom(x, self.flt.origin)
        self.flt.geom_grid = geom_short
        self.flt.dis_geom_grid = dis_geom_short
        self.green = Okada(self.flt, self.data, self.dict_green)
        
        m_gps, m_sar = [], []
        slip = [float(iss), float(ids), 0.0]
        if len(self.data.d_gps) > 0:
            m_gps = self.green.G.dot(slip)
        if len(self.data.d_sar) > 0:
            m_sar = self.green.G_sar[:,[0,1,2]].dot(slip)
            
        return np.hstack((m_gps, m_sar))


    def bayesian_inversion(self, method='emcee'):
        '''
        Run Bayesian Inversion.
        
        Input:
            method = emcee
        '''
        
        def log_prior(theta, leng, wid, dep, dip, strk, lat, lon, ss, ds):
            '''
            Calculate prior probability. In this case, uniform distribution is Implementated.
            
            Input:
                theta = variable
                leng  = min and max values of fault length
                wid   = min and max values of fault width
                dep   = min and max values of top depth of fault
                dip   = min and max values of fault dip
                strk  = min and max values of fault strike
                lat   = min and max values of fault location lat
                lon   = min and max values of fault location lon
            Output:
                ln of prior probability.
            '''
            ileng, iwid, idep, idip, istrk, ilat, ilon, iss, ids = theta

            if min(leng) < ileng < max(leng) and \
               min(strk) < istrk < max(strk) and \
               min(wid) < iwid < max(wid) and \
               min(dep) < idep < max(dep) and \
               min(dip) < idip < max(dip) and \
               min(lon) < ilon < max(lon) and \
               min(lat) < ilat < max(lat) and \
               min(ss) < iss < max(ss) and \
               min(ds) < ids < max(ds):
                return 0.0
            else:
                return -np.inf

        def log_likelihood(theta, leng, wid, dep, dip, strk, lat, lon, ss, ds):
            '''
            Calculate likelihood.
            Mod by Zhao Bin, Jul. 30, 2019. Consider the co-variance matrix of InSAR data.
            
            Input:
                theta = variable
                leng  = min and max values of fault length
                wid   = min and max values of fault width
                dep   = min and max values of top depth of fault
                dip   = min and max values of fault dip
                strk  = min and max values of fault strike
                lat   = min and max values of fault location lat
                lon   = min and max values of fault location lon
            Output:
                ln of likelihood.
            '''
            
            ileng, iwid, idep, idip, istrk, ilat, ilon, iss, ids = theta
            obs = np.hstack((self.data.d_gps, self.data.d_lev, self.data.d_sar))
            mod = self.forward(ileng, iwid, idep, idip, istrk, ilat, ilon, iss, ids)
            res = (obs-mod).reshape(len(obs),1)
            cov = self.data.cov_sar
            icov= np.linalg.inv(cov)
            
            print('wrms = %15.4f mm' %(np.sqrt(np.sum(res**2)/len(obs))))
            print('chi2/(N) = %15.4f mm^2' %(res.T.dot(icov).dot(res)[0,0]/len(obs)))
            
            return -0.5*res.T.dot(icov).dot(res)

        def log_posterior(theta, leng, wid, dep, dip, strk, lat, lon, ss, ds):
            '''
            Calculate posterior probability based on prior distribution and likelihood distribution.
            
            Input:
                theta = variable
                leng  = min and max values of fault length
                wid   = min and max values of fault width
                dep   = min and max values of top depth of fault
                dip   = min and max values of fault dip
                strk  = min and max values of fault strike
                lat   = min and max values of fault location lat
                lon   = min and max values of fault location lon
            Output:
                ln of posterior probability.
            '''
            ileng, iwid, idep, idip, istrk, ilat, ilon, iss, ids = theta
            
            if min(leng) < ileng < max(leng) and \
               min(strk) < istrk < max(strk) and \
               min(wid) < iwid < max(wid) and \
               min(dep) < idep < max(dep) and \
               min(dip) < idip < max(dip) and \
               min(lon) < ilon < max(lon) and \
               min(lat) < ilat < max(lat) and \
               min(ss) < iss < max(ss) and \
               min(ds) < ids < max(ds):
                return log_prior(theta, leng, wid, dep, dip, strk, lat, lon, ss, ds)+\
                        log_likelihood(theta, leng, wid, dep, dip, strk, lat, lon, ss, ds)
            else:
                return -np.inf

        leng = [self.xlb[0,0], self.xub[0,0]]
        wid  = [self.xlb[0,1], self.xub[0,1]]
        dep  = [self.xlb[0,2], self.xub[0,2]]
        dip  = [self.xlb[0,3], self.xub[0,3]]
        strk = [self.xlb[0,4], self.xub[0,4]]
        lat  = [self.xlb[0,5], self.xub[0,5]]
        lon  = [self.xlb[0,6], self.xub[0,6]]
        ss   = self.ss
        ds   = self.ds
        if method == 'emcee':
            import emcee, corner
            ndim     = 9
            nwalkers = 20
            nburn    = self.nburn
            nsteps   = self.nsteps
        
            starting_guess = np.random.random((nwalkers, ndim))
            starting_guess[:,0] = np.random.uniform(min(leng), max(leng), nwalkers)
            starting_guess[:,1] = np.random.uniform(min(wid), max(wid), nwalkers)
            starting_guess[:,2] = np.random.uniform(min(dep), max(dep), nwalkers)
            starting_guess[:,3] = np.random.uniform(min(dip), max(dip), nwalkers)
            starting_guess[:,4] = np.random.uniform(min(strk), max(strk), nwalkers)
            starting_guess[:,5] = np.random.uniform(min(lat), max(lat), nwalkers)
            starting_guess[:,6] = np.random.uniform(min(lon), max(lon), nwalkers)
            starting_guess[:,7] = np.random.uniform(min(ss), max(ss), nwalkers)
            starting_guess[:,8] = np.random.uniform(min(ds), max(ds), nwalkers)

            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, threads=1,
                    args=[leng, wid, dep, dip, strk, lat, lon, ss, ds])
            sampler.run_mcmc(starting_guess, nsteps)
            emcee_trace = sampler.chain[:,nburn:,:].reshape(-1, ndim)
            fig = corner.corner(emcee_trace)
            fig.savefig("posterior_distribution.png")
            np.save('mcmc_trace', emcee_trace)


    @staticmethod
    def trace2geom(tracefile, faultfile):
        '''
        '''
        import os
        flt = Fault(faultfile, 1, 1, False, origin=[])
        if os.path.isfile(tracefile):
            trace = np.load(tracefile)
            leng  = np.quantile(trace[:,0], [0.16, 0.5, 0.84])
            wid   = np.quantile(trace[:,1], [0.16, 0.5, 0.84])
            dep   = np.quantile(trace[:,2], [0.16, 0.5, 0.84])
            dip   = np.quantile(trace[:,3], [0.16, 0.5, 0.84])
            strk  = np.quantile(trace[:,4], [0.16, 0.5, 0.84])
            lat   = np.quantile(trace[:,5], [0.16, 0.5, 0.84])
            lon   = np.quantile(trace[:,6], [0.16, 0.5, 0.84])
            x     = np.array([leng, wid, dep, dip, strk, lat, lon]).T
            [geom, dis_geom] = Fault.x2geom(x, flt.origin)
            for i in range(len(geom)):
                print("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" %(geom[i,0], geom[i,1], geom[i,2], geom[i,3], geom[i,4],
                    geom[i,5], geom[i,6], 1, 1, 0))
        else:
            print(' The input trace file does not exist!')
