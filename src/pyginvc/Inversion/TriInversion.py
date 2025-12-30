#!/usr/bin/env python
# Written by Zhao Bin, when he was at UC Berkeley. Feb 14 2016
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
# Mod by Zhao Bin, Dec. 7, 2018. We use HDF5 to store solutions

import numpy as np
import logging, sys
from scipy import optimize
from numpy import vstack, hstack, zeros, arange, sqrt, diag, sum, genfromtxt

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")


class TriInversion(object):
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
        smoothfactor  = self.dict_weight['smoothfactor']
        dict_bound    = self.dict_bound
        self.data.set_wsar(wsar)

        # count number of smoothing factors
        scount         = 0
        
        # smooth factor
        smo_fact_start = smoothfactor[0]
        smo_fact_end   = smoothfactor[1]
        smo_fact_step  = smoothfactor[2]
    
        len_sar        = len(d_sar)
        len_geod       = len(d_gps) + len(d_lev)
        len_all        = len_geod   + len(d_sar)
        nf             = len(self.flt.element)
        d_lap          = zeros(3*nf)
        D              = hstack((d_gps, d_lev));
    
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
    
        # get the number of smoothing
        smo_facts      = arange(smo_fact_start, smo_fact_end, smo_fact_step)
        if len(smo_facts) == 0:
            logging.critical('Please check the input smoothing parameters')
            sys.exit()
        else:
            nsmooth    = len(smo_facts)
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
        G2I            = zeros((len(D)+nf, nf))
        sar_switch     = False
        dhat           = np.empty(0)
        
        if len(d_sar) > 0 and sar_switch == True:
            slip           = zeros((nsmooth, 3*nf+3))
            sig_slip       = zeros((nsmooth, 3*nf+3))
            slipb          = zeros((nsmooth, 3*nf+3))

        if 'sar_plane_range' in dict_bound.keys() and len(G_sar)>3:
            splane_range   = dict_bound['sar_plane_range']
            if len(splane_range) == 0:
                G_sar[:,-3:] = 0.0
    
        if dict_bound['bound_switch']:
            [bu, bl]       = self.set_bounds(len(self.flt.element), sar_switch, **dict_bound)  
       
            # print the status
            logging.info('Boundaries for constrained linear is constructed.')
        
        # for each smooth factor   
        for smo_fact in arange(smo_fact_start, smo_fact_end, smo_fact_step):
            
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
#                   G2I    = vstack((hstack((W.dot(G), zeros((len_geod,3)))), wsar*G_sar, hstack((G_laps, zeros((nf*3,3))))))
                    G2I    = vstack((W.dot(G), wsar*G_sar, G_laps))
                else:
#                   G2I    = vstack((wsar*G_sar, hstack((G_laps, zeros((nf*3,3))))))
                    G2I    = vstack((wsar*G_sar, G_laps))
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
                res = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls')
                slipb[scount,:] = res.x
                slip[scount,:] = slipb[scount,:]
                [m,n] = G2I.shape           
    
                # print the status
                logging.info('Constrained linear least square finished.')
    
            if len(d_sar) > 0:
                if len(G) == 0:
                    dhat = G_sar.dot(slip[scount,:])
                else:
#                   dhat = vstack((hstack((G, zeros((len_geod,3)))), G_sar)).dot(slip[scount,:])
                    dhat = vstack((G, G_sar)).dot(slip[scount,:])
    
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
            [Mo, Mw]  = self.moment_tri_element(self.flt.vertex_enu, self.flt.element, slip[scount,:], shearmodulus)
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
    def moment_tri_element(vertex, element, slip, shearmodulus):
        '''
        Calculate moment and magnitude for triangular element
        Written by Zhao Bin, Institute of Seismology, CEA. May 2017
        Mod by Zhao Bin, Apr. 24, 2019. Now 6.067 is used.
    
        Input:
            vertex   : coordinate for vertex in ENU
            element  : index of vertices for each triangular
            slip     : list/array slip value for all elememts
            shearmodulus:
        Output:
            Mo       : total moment
            Mw       : magnitude
        '''
        
        nelem   = len(element)
        Mo      = np.zeros(nelem)
        
        for i in range(nelem):
            id1 = element[i,0]-1
            id2 = element[i,1]-1
            id3 = element[i,2]-1
            
            len1 = np.sqrt((vertex[id1,0]-vertex[id2,0])**2 + 
                           (vertex[id1,1]-vertex[id2,1])**2 +
                           (vertex[id1,2]-vertex[id2,2])**2)
            len2 = np.sqrt((vertex[id1,0]-vertex[id3,0])**2 + 
                           (vertex[id1,1]-vertex[id3,1])**2 +
                           (vertex[id1,2]-vertex[id3,2])**2)
            len3 = np.sqrt((vertex[id3,0]-vertex[id2,0])**2 + 
                           (vertex[id3,1]-vertex[id2,1])**2 +
                           (vertex[id3,2]-vertex[id2,2])**2)
            s     = 0.5*(len1+len2+len3)
            area  = np.sqrt(s*(s-len1)*(s-len2)*(s-len3))
            Mo[i] = 1e6 * area * np.abs(np.sqrt(slip[3*i]**2 + slip[3*i+1]**2)) * shearmodulus
        Mo_total  = np.sum(Mo)/1000.0
        Mw_total  = 2.0/3.0*np.log10(Mo_total) - 6.067
        
        return Mo_total, Mw_total        
        
    
    @staticmethod
    def set_bounds(nelements, sar_switch, **varargin):
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
    
        slip_lb, slip_ub = '', ''
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

        nelems    = 3*nelements
        if sar_switch == True:
            nelems = nelems + 3

        bl = np.zeros((nelements,3))
        bu = np.zeros((nelements,3))

        bl[:,0] = ss_range[0]
        bu[:,0] = ss_range[1]
        bl[:,1] = ds_range[0]
        bu[:,1] = ds_range[1]
        bl[:,2] = op_range[0]
        bu[:,2] = op_range[1]

        bl = bl.flatten()
        bu = bu.flatten()
            
        # if slip bound files are exist, then overwrite the bound array bu and bl
        if slip_lb != '' and slip_ub != '':
            nelems    = 3*nelements
            bl[0:nelems] = genfromtxt(slip_lb, usecols = [3,4,5]).reshape(nelems)
            bu[0:nelems] = genfromtxt(slip_ub, usecols = [3,4,5]).reshape(nelems)
    
        return bu, bl 

        
if __name__ == '__main__':
    print('inverse')
