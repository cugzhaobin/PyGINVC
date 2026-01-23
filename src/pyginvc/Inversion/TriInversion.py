#!/usr/bin/env python
# Written by Zhao Bin, when he was at UC Berkeley. Feb 14 2016
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
# Mod by Zhao Bin, Dec. 7, 2018. We use HDF5 to store solutions

import numpy as np
import logging
from scipy.linalg import block_diag
from scipy import optimize

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
        G_gps_ramp    = self.green.G_gps_ramp
        G_sar_ramp    = self.green.G_sar_ramp
        smoothfactor  = self.dict_weight['smoothfactor']
        dict_bound    = self.dict_bound
        n_ramp        = G_gps_ramp.shape[1]+G_sar_ramp.shape[1]
        # self.data.set_wsar(wsar)

        n_sar = self.data.n_sar
        WSAR  = np.zeros_like(self.data.W_sar)
        assert len(n_sar) == len(wsar), "Number of wsar should be the same as sar data"
        for i in range(len(n_sar)):
            idx0 = sum(n_sar[:i])
            idx1 = sum(n_sar[:i+1])
            WSAR[idx0:idx1,idx0:idx1] = wsar[i]*self.data.W_sar[idx0:idx1,idx0:idx1]
    
        len_geod       = len(d_gps) + len(d_lev)
        len_all        = len_geod   + len(d_sar)
        nf             = self.flt.nf
        d_lap          = np.zeros(3*nf)
        D              = np.hstack((d_gps, d_lev))
    
        if len_all < 1:
            logging.fatal('No data to do inversion')
            return
        # get d2I and d2R
        # if len_geod > 0:
        #     if len_sar > 0:
        #         d2I    = hstack((W.dot(D), wsar*d_sar, d_lap))
        #         d2R    = hstack((W.dot(D), wsar*d_sar))
        #     else:
        #         d2I    = hstack((W.dot(D), d_lap))
        #         d2R    = W.dot(D)
        # else:
        #     if len_sar > 0:
        #         d2I    = hstack((wsar*d_sar, d_lap))
        #         d2R    = wsar*d_sar
        #     else:
        #         logging.warning('Final error, no data to do inversion')
        #         return
        d2I = np.hstack((W@D, WSAR@d_sar, d_lap))
        d2R = np.hstack((W@D, WSAR@d_sar))
    
        # get the number of smoothing
        # smo_facts      = arange(smo_fact_start, smo_fact_end, smo_fact_step)
        # if len(smo_facts) == 0:
        #     logging.critical('Please check the input smoothing parameters')
        #     sys.exit()
        # else:
        #     nsmooth    = len(smo_facts)
            
        # smooth factor
        smo_fact_start, smo_fact_end,  smo_fact_step = smoothfactor
        smo_facts      = np.arange(smo_fact_start, smo_fact_end, smo_fact_step)
        nsmooth        = len(smo_facts) if len(smo_facts)>0 else 1
        
        # init parameters
        slip           = np.zeros((nsmooth, 3*nf))
        sig_slip       = np.zeros((nsmooth, 3*nf))
        # slipb          = zeros((nsmooth, 3*nf))
        r              = np.zeros(len_all)
        misfit         = np.zeros(nsmooth)
        misfit_sar     = np.zeros(nsmooth)
        misfit_gps     = np.zeros(nsmooth)
        smoothness     = np.zeros(nsmooth)
        ruff           = np.zeros(nsmooth)
        moment         = np.zeros((nsmooth,2))
        ramp           = np.zeros((nsmooth, n_ramp))
        # ss_slip        = zeros((nsmooth, nf))
        # ss_sig_slip    = zeros((nsmooth, nf))
        # dip_slip       = zeros((nsmooth, nf))
        # dip_sig_slip   = zeros((nsmooth, nf))
        # openning       = zeros((nsmooth, nf))
        # op_sig_slip    = zeros((nsmooth, nf))
        # G2I            = zeros((len(D)+nf, nf))   
        sar_switch     = False
        
        # if len(d_sar) > 0 and sar_switch == True:
        #     slip           = zeros((nsmooth, 3*nf+3))
        #     sig_slip       = zeros((nsmooth, 3*nf+3))
        #     slipb          = zeros((nsmooth, 3*nf+3))

        # if 'sar_plane_range' in dict_bound.keys() and len(G_sar)>3:
        #     splane_range   = dict_bound['sar_plane_range']
        #     if len(splane_range) == 0:
        #         G_sar[:,-3:] = 0.0
    
        if dict_bound['bound_switch']:
            [bu, bl] = self.set_bounds(len(self.flt.element), sar_switch, **dict_bound)
            if n_ramp>0:
                bl = np.hstack([bl, np.full(n_ramp, -np.inf)])
                bu = np.hstack([bu, np.full(n_ramp,  np.inf)])
       
            # print the status
            logging.info('Boundaries for constrained linear is constructed.')
        
        # for each smooth factor   
        for i, smo_fact in enumerate(smo_facts):
            
            # print the status
            logging.info('Inverting for No. %d with smooth factor %f' %(i+1, smo_fact))
            
            # scale the G_lap
            G_laps     = smo_fact*G_lap
                        
            # if we use InSAR data to do inversion
#             if len_sar > 0:
#                 # MOD By Zhao Bin, Jan 7 2017
#                 if len(wsar) == len_sar:
#                     wsar = wsar.reshape((len_sar,1))
    
#                 # and if we use GPS and Level data to do inversion, combing all Green functions
#                 if len_geod > 0:
# #                   G2I    = vstack((hstack((W.dot(G), zeros((len_geod,3)))), wsar*G_sar, hstack((G_laps, zeros((nf*3,3))))))
#                     G2I    = vstack((W.dot(G), wsar*G_sar, G_laps))
#                 else:
# #                   G2I    = vstack((wsar*G_sar, hstack((G_laps, zeros((nf*3,3))))))
#                     G2I    = vstack((wsar*G_sar, G_laps))
#             else:
#                 G2I        = vstack((W.dot(G), G_laps))

            if len(G) == 0:
                G = G.reshape(-1,3*nf)
            if len(G_sar) == 0:
                G_sar = G_sar.reshape(-1,3*nf)
            WG     = W@G
            WG_sar = WSAR @ G_sar
            if len(G) == 0:
                WG = WG.reshape(-1, nf*3)
            if len(G_sar) == 0:
                WG_sar = WG_sar.reshape(-1, nf*3)
            G2I = np.vstack((WG, WG_sar, G_laps))
            Gramp = np.vstack([ block_diag(G_gps_ramp, G_sar_ramp),
                               np.zeros((G_laps.shape[0], n_ramp)) ])
            G2I = np.column_stack((G2I, Gramp))
    
            # invert the matrix G2I                    
            Ginv        = np.linalg.pinv(G2I)
            slip[i]     = Ginv.dot(d2I)[0:3*nf]
            cov_slip    = Ginv.dot(Ginv.T)[0:3*nf,0:3*nf]
            sig_slip[i] = np.sqrt(np.diag(cov_slip))

            logging.info('Unconstrained linear least square finished.')
    
            # if constrained linear least square method is used
            if dict_bound['bound_switch']:
                res = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls')
                slip[i] = res.x[0:3*nf]
                ramp[i] = res.x[3*nf:] 
    
                logging.info('Constrained linear least square finished.')
    
#             if len(d_sar) > 0:
#                 if len(G) == 0:
#                     dhat = G_sar.dot(slip[scount,:])
#                 else:
# #                   dhat = vstack((hstack((G, zeros((len_geod,3)))), G_sar)).dot(slip[scount,:])
#                     dhat = vstack((G, G_sar)).dot(slip[scount,:])
    
#             else:
#                 dhat = G.dot(slip[scount,:])

            dhat = np.column_stack((
                np.vstack((G, G_sar)), block_diag(G_gps_ramp, G_sar_ramp)
                )) @ res.x
            
            if len_geod > 0:
                r[0:len_geod] = d2R[0:len_geod] - W.dot(dhat[0:len_geod])
                r_gps         = r[0:len_geod]

            r[len_geod:len_all] = d2R[len_geod:len_all] - WSAR @ dhat[len_geod:len_all]
            r_sar               = r[len_geod:len_all]
    #       rnk[scount] = rank(G2I)
    
    
            if len(d_gps) > 0:
                n_gps = len(d_gps)
                logging.info('GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(r_gps.dot(r_gps)))
                logging.info('GPS: WRSS/(N) %f' %(r_gps.dot(r_gps)/n_gps))
                misfit_gps[i] = r_gps.dot(r_gps)
                
                
            if len(d_sar) > 0:
                n_sar = len(d_sar)
                logging.info('SAR Weighted Residual Sum of Squares (WRSS) %10.3f' %(r_sar.dot(r_sar)))
                logging.info('SAR: WRSS/(N) %f' %(r_sar.dot(r_sar)/n_sar))
                misfit_sar[i] = r_sar.dot(r_sar)
            
            misfit[i] = r.dot(r)
            
            smoothness[i] = sum( (smo_fact*slip[i].dot(G_lap)) * (smo_fact*slip[i].dot(G_lap)) )
            ruff[i]       = sum( (slip[i].dot(G_lap)) * (slip[i].dot(G_lap)) )
            
            
            # print slip
            # for k in range(nf):
            #     ss_slip[scount, k]      = slip[scount, 3*k]
            #     ss_sig_slip[scount, k]  = sig_slip[scount, 3*k]
            #     dip_slip[scount, k]     = slip[scount, 3*k+1]
            #     dip_sig_slip[scount, k] = sig_slip[scount, 3*k+1]
            #     openning[scount, k]     = slip[scount, 3*k+2]
            #     op_sig_slip[scount, k]  = sig_slip[scount, 3*k+2]
    
    
            # compute the moment
            shearmodulus = self.green.modulus
            [Mo, Mw]     = self.moment_tri_element(self.flt.vertex_enu, self.flt.element, slip[i], shearmodulus)
            moment[i]    = np.array([Mo, Mw])            
            
            # print the status
            logging.info('Geodetic Moment Magnitude M0 = %E' %(Mo))
            logging.info('Geodetic Moment Magnitude Mw = %f' %(Mw))
    
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
        self.WSAR       = WSAR
        self.ramp       = ramp
        
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
        
        idx1    = element[:,0]-1
        idx2    = element[:,1]-1
        idx3    = element[:,2]-1
        v1      = vertex[idx1]
        v2      = vertex[idx2]
        v3      = vertex[idx3]

        vec1    = v2 - v1
        vec2    = v3 - v1
        cross   = np.cross(vec1, vec2)
        area    = 0.5 * np.linalg.norm(cross, axis=1)
        slip    = slip.reshape(-1,3)
        slip_mag= np.linalg.norm(slip[:,0:2], axis=1)
        Mo      = 1e6 * area * slip_mag * shearmodulus
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
        # if sar_switch == True:
        #     nelems = nelems + 3

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
            bl[0:nelems] = np.genfromtxt(slip_lb, usecols = [3,4,5]).reshape(nelems)
            bu[0:nelems] = np.genfromtxt(slip_ub, usecols = [3,4,5]).reshape(nelems)
    
        return bu, bl 

        
if __name__ == '__main__':
    print('inverse')
