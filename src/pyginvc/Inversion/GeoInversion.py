#!/usr/bin/env python
# Written by Zhao Bin, when he was at UC Berkeley. Feb 14 2016
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
# Mod by Zhao Bin, Dec. 7, 2018. We use HDF5 to store solutions

import logging
import numpy as np
from scipy import optimize
from numpy import array, genfromtxt, remainder, zeros
from numpy import vstack, hstack, arange, sqrt, diag, sum
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")


class GeoInversion(object):
    '''
    GeoInversion represents a class of slip inversion using geodetic data
    '''
    
    def __init__(self, flt, data, green, lap, dict_weight, dict_bound):
        '''
        Constructor.
        
        Input:
            flt         = an instance of class Fault
            data        = an instance of class GeoData
            green       = an instance of class GreenFunction
            lap         = an instance of class Lapacian
            dict_weight = a dict
            dict_bound  = a dict
        '''
        
        self.flt         = flt
        self.data        = data
        self.green       = green
        self.lap         = lap
        self.dict_weight = dict_weight
        self.dict_bound  = dict_bound
       
        self.inversion()
        return
 
            
    def inversion(self):
        '''
        Invert for slip model
        Mod by Zhao Bin, Mar 14, 2020. Rescale W with wgps
        '''

        d_gps         = self.data.d_gps
        d_lev         = self.data.d_lev
        d_sar         = self.data.d_sar
        wgps          = self.dict_weight['wgps']
        W             = wgps*self.data.W
        wsar          = self.dict_weight['wsar']
        G             = self.green.G
        G_sar         = self.green.G_sar
        G_lap         = self.lap.G_lap
        dis_geom_grid = self.flt.dis_geom_grid
        nsegs         = self.flt.nsegs
        ndeps         = self.flt.ndeps
        smoothfactor  = self.dict_weight['smoothfactor']
        dict_bound    = self.dict_bound
        
        # set SAR weights in data object
        self.data.set_wsar(wsar)

        # count number of smoothing factors
        scount         = 0
        smo_fact_start, smo_fact_end,  smo_fact_step = smoothfactor
        smo_facts      = arange(smo_fact_start, smo_fact_end, smo_fact_step)
        nsmooth        = len(smo_facts) if len(smo_facts)>0 else 1

        # Data dimensions
        len_sar        = len(d_sar)
        len_geod       = len(d_gps) + len(d_lev)
        len_all        = len_geod   + len(d_sar)
        nf             = nsegs*ndeps
        d_lap          = zeros(3*nf)
        D              = hstack((d_gps, d_lev))
    
        # get d2I and d2R
        if len_geod > 0:
            if len_sar > 0:
                d2I    = hstack((W.dot(D), wsar*d_sar, d_lap))
                d2R    = hstack((W.dot(D), wsar*d_sar))
            else:
                d2I    = hstack((W.dot(D), d_lap))
                d2R    = W.dot(D)
        else:
            if len_sar > 0:
                d2I    = hstack((wsar*d_sar, d_lap))
                d2R    = wsar*d_sar
            else:
                logging.warning('Final error, no data to do inversion')
                return
    
        
        
        # init parameters
        slip           = zeros((nsmooth, 3*nf))
        sig_slip       = zeros((nsmooth, 3*nf))
        slipb          = zeros((nsmooth, 3*nf))
        r              = zeros(len_all)
        misfit         = zeros(nsmooth)
        misfit_sar     = zeros(nsmooth)
        misfit_gps     = zeros(nsmooth)
        smoothness     = zeros(nsmooth)
        ruff           = zeros(nsmooth)
        moment         = zeros((nsmooth,2))
        ss_slip        = zeros((nsmooth, nf))
        ss_sig_slip    = zeros((nsmooth, nf))
        dip_slip       = zeros((nsmooth, nf))
        dip_sig_slip   = zeros((nsmooth, nf))
        openning       = zeros((nsmooth, nf))
        op_sig_slip    = zeros((nsmooth, nf))
        #G2I            = zeros((len(D)+nf, nf))
        sar_switch     = 0
        
        if len(d_sar) > 0:
            slip           = zeros((nsmooth, 3*nf+3))
            sig_slip       = zeros((nsmooth, 3*nf+3))
            slipb          = zeros((nsmooth, 3*nf+3))
            sar_switch     = 1

        if 'sar_plane_range' in dict_bound.keys() and len(G_sar)>3:
            splane_range   = dict_bound['sar_plane_range']
            if len(splane_range) == 0:
                G_sar[:,-3:] = 0.0
    
        if dict_bound['bound_switch']:
            [bu, bl]       = self.set_bounds(nsegs, ndeps, sar_switch, **dict_bound)  
       
            # print the status
            logging.info('Boundaries for constrained linear is constructed.')
        
        # for each smooth factor   
        for smo_fact in smo_facts:
            
            # print the status
            logging.info('Inverting for No. %d with smooth factor %f' %(scount+1, smo_fact))
            
            # scale the G_lap
            G_laps     = smo_fact*G_lap
                        
            # if we use InSAR data to do inversion
            if len_sar > 0:
                # MOD By Zhao Bin, Jan 7 2017
                if len(wsar) == len_sar:
                    wsar = wsar.reshape((len_sar,1))
    
                # and if we use GPS and Level data to do inversion, combing all Green functions
                if len_geod > 0:
                    G2I    = vstack((hstack((W.dot(G), zeros((len_geod,3)))), wsar*G_sar, hstack((G_laps, zeros((nf*3,3))))))
                else:
                    G2I    = vstack((wsar*G_sar, hstack((G_laps, zeros((nf*3,3))))))
            else:
                G2I        = vstack((W.dot(G), G_laps))
    
            # invert the matrix G2I                    
            Ginv           = np.linalg.pinv(G2I)
            # compute the slip
            slip[scount,:] = Ginv.dot(d2I)
            cov_slip       = Ginv.dot(Ginv.T)
            sig_slip[scount,:] = sqrt(diag(cov_slip))
    
            # print the status
            logging.info('Unconstrained linear least square finished.')
    
            # if constrained linear least square method is used
            if dict_bound['bound_switch']:
                res = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls', lsmr_tol='auto')
                slipb[scount,:] = res.x
                slip[scount,:]  = slipb[scount,:]       
    
                # print the status
                logging.info('Constrained linear least square finished.')
            dhat = None 
                
            if len(d_sar) > 0:
                if len(G) == 0:
                    dhat = G_sar.dot(slip[scount,:])
                else:
                    dhat = vstack((hstack((G, zeros((len_geod,3)))), G_sar)).dot(slip[scount,:])
    
            else:
                dhat = G.dot(slip[scount,:])
            
            if len_geod > 0:
                r[0:len_geod] = d2R[0:len_geod] - W.dot(dhat[0:len_geod])
                r_gps         = r[0:len_geod]
            
            # By Zhao Bin Jan 7 2017
            if len(wsar) == len_sar:  wsar = wsar.flatten()
            r[len_geod:len_all] = d2R[len_geod:len_all]-wsar*dhat[len_geod:len_all]
            r_sar               = r[len_geod:len_all]
    #       rnk[scount] = rank(G2I)
    
            if len(d_gps) > 0:
                n_gps = len(d_gps)
                logging.info('GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(r_gps.dot(r_gps)))
                logging.info('GPS: WRSS/(N) %f' %(r_gps.dot(r_gps)/n_gps))
                misfit_gps[scount] = r_gps.dot(r_gps)
                
                
            if len(d_sar) > 0:
                n_sar = len(d_sar)
                logging.info('SAR Weighted Residual Sum of Squares (WRSS) %10.3f' %(r_sar.dot(r_sar)))
                logging.info('SAR: WRSS/(N) %f' %(r_sar.dot(r_sar)/n_sar))
                misfit_sar[scount] = r_sar.dot(r_sar)
            
            misfit[scount] = r.dot(r)
            
            smoothness[scount] = sum( (smo_fact*slip[scount,0:3*nf].dot(G_lap)) * (smo_fact*slip[scount,0:3*nf].dot(G_lap)) )
            ruff[scount]       = sum( (slip[scount,0:3*nf].dot(G_lap)) * (slip[scount,0:3*nf].dot(G_lap)) )
            
            
            # print slip
            for k in range(nf):
                ss_slip[scount, k]      = slip[scount, 3*k]
                ss_sig_slip[scount, k]  = sig_slip[scount, 3*k]
                dip_slip[scount, k]     = slip[scount, 3*k+1]
                dip_sig_slip[scount, k] = sig_slip[scount, 3*k+1]
                openning[scount, k]     = slip[scount, 3*k+2]
                op_sig_slip[scount, k]  = sig_slip[scount, 3*k+2]
    
    
            # compute the moment
            shearmodulus= self.green.modulus
            [Mo, Mw]  = self.Moment(dis_geom_grid, slip[scount,:], shearmodulus)
            moment[scount] = np.array([Mo, Mw])            
            
            # print the status
            logging.info('Geodetic Moment Magnitude M0 = %E' %(Mo))
            logging.info('Geodetic Moment Magnitude Mw = %f' %(Mw))
    
            scount = scount+1
    
        self.ruff       = ruff
        self.misfit_gps = misfit_gps
        self.misfit_sar = misfit_sar
        self.misfit     = misfit
        self.slip       = slip
        self.sig_slip   = sig_slip
        self.r          = r
        self.dhat       = dhat
        self.moment     = moment
        self.smo_facts  = smo_facts
        self.smoothness = smoothness
        self.wsar       = wsar
        self.data.wsar  = wsar
        
        return
        

            
    @staticmethod 
    def Moment(dis_geom_grid, slip, shearmodulus):
        '''
        Estimate moment and moment magnitude
        Mod by Zhao Bin, Apr. 24, 2019. Now 6.067 is used.
    
        Input:
            dis_geom_grid = fault patches
            slip          = estimated fault slip
            shearmodulus  = strength of crust    
            
        Output:
            Mo            = total moment in Nm
            Mw            = moment magnitude
        
        '''    
        
        nf = len(dis_geom_grid)
        Mo = np.zeros(nf)
        
        for i in range(nf):
            Mo[i] = (1e3*dis_geom_grid[i,0] *
                     1e3*dis_geom_grid[i,1] *
                     np.abs(np.sqrt(slip[i*3]**2+slip[3*i+1]**2)) * shearmodulus)
        
        Mo_total = np.sum(Mo)/1e3
        Mw_total = 2.0/3.0*np.log10(Mo_total) - 6.067
        
        return Mo_total, Mw_total
                
    
    @staticmethod
    def set_bounds(nsegs, ndeps, sar_switch, **varargin):
        '''
        Set bounds for estimated parameters - slip[nsegs*ndeps*3, 3]+sar_plane
        Written by Zhao Bin @ UC Berkeley April 18, 2016
        MOD by Zhao Bin, Dec 22, 2016. read lower and upper slip bound from files
        
        Input:
            nsegs       = number of patches along strike direction
            ndeps       = number of patches along dip direction
            sar_switch  = 1 if InSAR data, 0 if no InSAR data
            varargin    = dictionary, containing all parameters
        Output:
            bu          = upper bound of parameter
            bl          = lower bound of parameter
        '''            
    
        nelems    = 3*nsegs*ndeps
        if sar_switch == True:
            nelems = nelems + 3
            
        # Initialize variables
        ss_range_c = 0
        ss_range   = array([0,0])
        ds_range_c = 0
        ds_range   = array([0,0])
        op_range_c = 0
        op_range   = array([0,0])
        surf_slip  = zeros(3*nsegs)
        surf_slip_c  = 0
        s_ss_range   = 0
        s_ss_range_c = 0
        s_ds_range   = 0
        s_ds_range_c = 0
        s_op_range_c = 0
        s_op_range   = 0
        bot_slip     = zeros(3*nsegs)
        bot_slip_c   = 0
        b_ss_range_c = 0 
        b_ds_range_c = 0 
        b_op_range_c = 0
        bu           = zeros(nelems)
        bl           = zeros(nelems)
        slip_lb      = ''
        slip_ub      = ''
                
        for (k,v) in varargin.items():
            if k == 'sar_plane_range':
                splane_range = v
            elif k == 'ss_range':
                ss_range = v
    #           ss_range_c = 1
            elif k == 'ds_range':
                ds_range = v
    #           ds_range_c = 1
            elif k == 'op_range':
                op_range = v
    #           op_range_c = 1
            elif k == 'surf_slip':
                surf_slip = v
                surf_slip_c = 1
            elif k == 's_ss_range':
                s_ss_range = v
    #           s_ss_range_c = 1
            elif k == 's_ds_range':
                s_ds_range = v
    #           s_ds_range_c = 1
            elif k == 's_op_range':
                s_op_range = v
    #           s_op_range_c = 1
            elif k == 'bot_slip':
                bot_slip = v
                bot_slip_c = 1
            elif k == 'b_ss_range':
                b_ss_range = v
    #           b_ss_range_c = 1
            elif k == 'b_ds_range':
                b_ds_range = v
    #           b_ds_range_c = 1
            elif k == 'b_op_range':
                b_op_range = v
    #           b_op_range_c = 1
            elif k == 'slip_lb':
                slip_lb = v
            elif k == 'slip_ub':
                slip_ub = v
            
        for i in range(0, nelems):
            # for strike slip component
            if remainder(i,3) == 0:
                if i <= nsegs*3:
                    if surf_slip_c == 1:
                         bu[i] = surf_slip[i]+s_ss_range
                         bl[i] = surf_slip[i]-s_ss_range
                         if bu[i] > ss_range[1]:
                             bu[i] = ss_range[1]
                         if bl[i] < ss_range[1]:
                             bl[i] = ss_range[1]
                    else:
                         bu[i] = ss_range[1]
                         bl[i] = ss_range[0]
                elif i>nsegs*(ndeps-1)*3 & i<=nsegs*ndeps*3:
                     if (bot_slip_c == 1):
                         bu[i] = bot_slip[i-nsegs*(ndeps-1)*3] + b_ss_range
                         bl[i] = bot_slip[i-nsegs*(ndeps-1)*3] - b_ss_range
                         if bu[i] > ss_range[1]:
                             bu[i] = ss_range[1]
                         if bl[i] < ss_range[0]:
                             bl[i] = ss_range[1]
                     else:
                         bu[i] = ss_range[1]
                         bl[i] = ss_range[0]
                else:
                     bu[i] = ss_range[1]
                     bl[i] = ss_range[0]
                
            # for dip slip component
            if remainder(i,3) == 1:
                 if i <= nsegs*3:
                     if surf_slip_c == 1:
                         bu[i] = surf_slip[i]+s_ds_range
                         bl[i] = surf_slip[i]-s_ds_range
                         if bu[i] > ds_range[1]:
                             bu[i] = ds_range[1]
                         if bl[i] < ds_range[1]:
                             bl[i] = ds_range[1]
                     else:
                         bu[i] = ds_range[1]
                         bl[i] = ds_range[0]
                 elif i>nsegs*(ndeps-1)*3 & i<=nsegs*ndeps*3:
                     if (bot_slip_c == 1):
                         bu[i] = bot_slip[i-nsegs*(ndeps-1)*3] + b_ds_range
                         bl[i] = bot_slip[i-nsegs*(ndeps-1)*3] - b_ds_range
                         if bu[i] > ds_range[1]:
                             bu[i] = ds_range[1]
                         if bl[i] < ds_range[0]:
                             bl[i] = ds_range[1]
                     else:
                         bu[i] = ds_range[1]
                         bl[i] = ds_range[0]
                 else:
                     bu[i] = ds_range[1]
                     bl[i] = ds_range[0]           
            # openning component uses op_range, s_op_range, b_op_range   
            if remainder(i,3) == 2:
                 if i <= nsegs*3:
                     if surf_slip_c == 1:
                         bu[i] = surf_slip[i]+s_op_range
                         bl[i] = surf_slip[i]-s_op_range
                         if bu[i] > op_range[1]:
                             bu[i] = op_range[1]
                         if bl[i] < op_range[1]:
                             bl[i] = op_range[1]
                     else:
                         bu[i] = op_range[1]
                         bl[i] = op_range[0]
                 elif i>nsegs*(ndeps-1)*3 & i<=nsegs*ndeps*3:
                     if (bot_slip_c == 1):
                         bu[i] = bot_slip[i-nsegs*(ndeps-1)*3] + b_op_range
                         bl[i] = bot_slip[i-nsegs*(ndeps-1)*3] - b_op_range
                         if bu[i] > op_range[1]:
                             bu[i] = op_range[1]
                         if bl[i] < op_range[0]:
                             bl[i] = op_range[1]
                     else:
                         bu[i] = op_range[1]
                         bl[i] = op_range[0]
                 else:
                     bu[i] = op_range[1]
                     bl[i] = op_range[0]                
    
            if i>nsegs*ndeps*3:
                if len(splane_range) == 0:
                    logging.warning('A default range of [-1e3, 1e3] for SAR plane is assigned.')
                    splane_range = [-1e3, 1e3]
                bu[i] = splane_range[1]
                bl[i] = splane_range[0]            
                
        # added by zhao bin Dec 22, 2016. 
        # if slip bound files are exist, then overwrite the bound array bu and bl
        if slip_lb != '' and slip_ub != '':
            nelems       = 3*nsegs*ndeps
            bl[0:nelems] = genfromtxt(slip_lb, usecols = [7,8,9]).reshape(nelems)
            bu[0:nelems] = genfromtxt(slip_ub, usecols = [7,8,9]).reshape(nelems)
            logging.info('The lower bound and upper bound of slip are assgined from {} and {}.'.format(slip_lb, slip_ub))
    
        return bu, bl 

        
if __name__ == '__main__':
    print('inverse')
