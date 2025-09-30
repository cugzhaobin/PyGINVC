#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Determine the optimal fault geometry using Bayesian inversion methodology or other methods.

Created on July. 10, 2019.

@author: Zhao Bin, Institute of Seismology, CEA.
"""

from scipy.linalg import block_diag
import numpy as np
from pyginvc.Geometry.Fault import Fault
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Greens.Okada import Okada

class GeomOptimization:
    '''
    GeomOptimization is a class representing determination of fault geometry from geodetic data.
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
        self.dict_weight= dict_weight

        #
        # data weighting
        #
        self.wgps       = dict_weight['wgps']
        self.wsar       = dict_weight['wsar']
        sigma           = dict_weight['sigma']
        lamda           = dict_weight['lamda']
        self.data.MakeCVM(sigma=sigma, lamda=lamda)

        #
        # bayesian
        #
        self.nsteps     = dict_param['nsteps']
        self.nburn      = dict_param['nburn']
        self.ss         = dict_param['ss']
        self.ds         = dict_param['ds']
        self.mpi        = bool(dict_param['mpi'])

    def getdata(self):
        '''
        Get geodetic data
        '''
        return np.hstack((self.d_gps, self.d_lev, self.d_sar))


    def forward(self, ileng, iwid, idep, idip, istrk, ilat, ilon, iss, ids):
        '''
        Calcualte the prediction from input fault geometry and slips.
        
        Input:
            ileng  = length of fault patch
            iwid   = width of fault patch
            idep   = top depth of fault patch
            idip   = dip angle of fault patch
            istrk  = strike angle of fault patch
            ilat   = latitutde of the middle point of top edge of fault plane
            ilon   = longitude of the middle point of top edge of fault plane
            iss    = strike-slip
            ids    = dip-slip
            
        Output:
            predicted surface displacement 
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
            cov = block_diag(self.data.cov_gps, self.data.cov_lev, self.data.cov_sar)
            icov= np.linalg.inv(cov)
            
            print('rms      = %15.4f mm' %(np.sqrt(np.sum(res**2)/len(obs))))
            print('chi2/(N) = %15.4f' %(res.T.dot(icov).dot(res)[0,0]/len(obs)))
            
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
                ss    = min and max values of strike-slip
                ds    = min and max values of dip-slip
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
            import emcee, corner, sys
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

            if self.mpi is True:
                try:
                    from schwimmbad import MPIPool
                except ImportError:
                    pass
                with MPIPool() as pool:
                    if not pool.is_master():
                        pool.wait()
                        sys.exit(0)
                    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, threads=1,
                            args=[leng, wid, dep, dip, strk, lat, lon, ss, ds])
                    sampler.run_mcmc(starting_guess, nsteps)
            else:
                sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, threads=1,
                        args=[leng, wid, dep, dip, strk, lat, lon, ss, ds])
                sampler.run_mcmc(starting_guess, nsteps)
            emcee_trace = sampler.chain[:,nburn:,:].reshape(-1, ndim)
            fig = corner.corner(emcee_trace)
            fig.savefig("posterior_distribution.png")
            np.save('mcmc_trace', emcee_trace)


    def trace2geom(tracefile, faultfile):
        '''
        Determine the fault geometry and uncertainties from Bayesian inversion.
        
        Input:
            tracefile  = the default file name is mcmc_trace.
            faultfile  = file name of faultgeom format   
            
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


    def run_geomest(self, method='differential_evolution'):
        '''
        Determine the fault geometry using functions of scipy.optimize package.

        Input:
            method = minimize/fmin/fmin_ncg/fmin_slsqp/differential_evolution
	    '''
        from  scipy import optimize


        xlb    = self.xlb
        xub    = self.xub
        option = (self.flt, self.data, self.dict_green, self.dict_weight)
        bnds   = [(xlb[0,0],xub[0,0]), (xlb[0,1],xub[0,1]), (xlb[0,2],xub[0,2]), \
                  (xlb[0,3],xub[0,3]), (xlb[0,4],xub[0,4]), (xlb[0,5],xub[0,5]), (xlb[0,6],xub[0,6])]

        x0     = self.x
        if (method == 'minimize'):
            res  = optimize.minimize(self.disl_wr, x0, args=option, method='TNC', bounds=bnds)
            x    = res.x
        elif (method == 'fmin'):
	        res  = optimize.fmin(self.disl_wr, x0, args=option)
	        x    = res
        elif (method == 'fmin_ncg'):
	        res  = optimize.fmin(self.disl_wr, x0, args=option)
	        x    = res
        elif (method == 'fmin_slsqp'):
	        res  = optimize.fmin_slsqp(self.disl_wr, x0, args=option, bounds=bnds)
	        x    = res.x
        elif (method == 'differential_evolution'):
	        res  = optimize.differential_evolution(self.disl_wr, bnds, args=option, strategy='rand1bin', disp=True, polish=True)
	        x    = res.x
	
        print('**************************************')
        print('             length  width     depth     dip    strike   lat       lon')
        print(" Intial: %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" %(x0[0,0], x0[0,1], x0[0,2], x0[0,3], x0[0,4], x0[0,5], x0[0,6]))
        print(" Lbound: %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" %(xlb[0,0], xlb[0,1], xlb[0,2], xlb[0,3], xlb[0,4], xlb[0,5], xlb[0,6]))
        print(" Ubound: %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" %(xub[0,0], xub[0,1], xub[0,2], xub[0,3], xub[0,4], xub[0,5], xub[0,6]))
        print('Optimal: %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' %(x[0], x[1], x[2], x[3], x[4], x[5], x[6]))
        x = np.array(x).reshape((1,7))
        geom, dis_geom = Fault.x2geom(x, self.flt.origin)
        print('%f %f %f %f %f %f %f %f %f %f' %(geom[0][0], geom[0][1], geom[0][2], geom[0][3], geom[0][4], geom[0][5], geom[0][6], 0.0, 0.0, 0.0))


    @staticmethod
    def disl_wr(x, flt, data, dict_green, dict_weight):
        '''
        Object function
        
        Input:
            x           = 
            flt         = an instance of class Fault
            data        = an instance of class GeoData
            dict_green  = a dict for calculating Green's function
            dict_weight = a dict
        '''
        
        import datetime        
        import numpy as np
        from numpy import vstack, hstack, zeros

    	
        nf                           = 1
        x                            = np.array(x).reshape((1,7))
        [geom_short, dis_geom_short] = Fault.x2geom(x, flt.origin)
        flt.geom_grid                = geom_short
        flt.dis_geom_grid            = dis_geom_short
        green                        = Okada(flt, data, dict_green)
        G                            = green.G
        G_sar                        = green.G_sar
    	        
        len_gps                      = len(data.d_gps)
        len_sar                      = len(data.d_sar)
        len_lev                      = len(data.d_lev)
        len_geod                     = len_gps + len_lev
        len_all                      = len_geod + len_sar
        D                            = hstack((data.d_gps, data.d_lev));
        W                            = data.W
        wsar                         = dict_weight['wsar']
        d_sar                        = data.d_sar
    		
    	# get d2I and d2R
        if len_geod > 0:
    	    if len_sar > 0:
    		    d2I            = hstack((W.dot(D), wsar*d_sar))
    		    d2R            = hstack((W.dot(D), wsar*d_sar))
    	    else:
    		    d2I            = W.dot(D)
    		    d2R            = W.dot(D)
        else:
    	    if len_sar > 0:
    		    d2I            = wsar*d_sar
    		    d2R            = wsar*d_sar
    	    else:
    		    print(' WARNING: final error, no data to do inversion')
    		    return
    		
    	# init parameters
        r              = zeros(len_all)
        G2I            = zeros((len(D)+nf, nf))
        slip           = zeros(3*nf)
    		    
    	        
        if len(d_sar) > 0:
    	    sar_switch = 1
    	    slip           = zeros(3*nf+3)
    #	    sig_slip       = zeros(3*nf+3)
    #	    slipb          = zeros(3*nf+3)
        else:
    	    sar_switch = 0
    		       
    	# if we use InSAR data to do inversion
        if len(d_sar) > 0:
    		# and if we use GPS and Level data to do inversion, combing all Green functions
    	    if len_geod > 0:
    		    G2I    = vstack((hstack((W.dot(G), zeros((len_geod,3)))), wsar*G_sar))
    	    else:
    		    G2I    = wsar*G_sar
        else:
            G2I        = W.dot(G)
    		
    	# invert the matrix G2I                    
        Ginv           = np.linalg.pinv(G2I)
        # compute the slip
        slip[:]        = Ginv.dot(d2I)
    #    cov_slip       = Ginv.dot(Ginv.T)
    	        
        if len(d_sar) > 0:
    	    if len(G) == 0:
    		    dhat = G_sar.dot(slip)
    	    else:
    		    dhat = vstack((hstack((G, zeros((len_geod,3)))), G_sar)).dot(slip)
        else:
            dhat = G.dot(slip)
    		        
        if len_geod > 0:
            r[0:len_geod] = d2R[0:len_geod] - W.dot(dhat[0:len_geod])
            r_gps = r[0:len_geod]
    		    
        r[len_geod:len_all]=d2R[len_geod:len_all] - wsar*dhat[len_geod:len_all]
        r_sar = r[len_geod:len_all]
    	
        now = datetime.datetime.now()
        now = now.strftime('%y%m%d:%H%M:%S')
        if len_gps > 0:
    	    n_gps = len_gps
    	    print(' STATUS  :'+now+' Inverse/disl_wr: GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(r_gps.dot(r_gps)))
    	    print(' STATUS  :'+now+' Inverse/disl_wr: GPS: WRSS/(N) %f' %(r_gps.dot(r_gps)/n_gps))
        if len_sar > 0:
    	    n_sar = len_sar
    	    print(' STATUS  :'+now+' Inverse/disl_wr: SAR Weighted Residual Sum of Squares (WRSS) %10.3f' %(r_sar.dot(r_sar)))
    	    print(' STATUS  :'+now+' Inverse/disl_wr: SAR: WRSS/(N) %f' %(r_sar.dot(r_sar)/n_sar))
    	        
        misfit = r.dot(r)
        return misfit
