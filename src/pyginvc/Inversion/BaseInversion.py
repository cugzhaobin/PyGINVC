#!/usr/bin/env python
import logging
import numpy as np
from scipy import optimize
from pyginvc.GeoData.GeoData import GeoData

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S"
)

class BaseInversion:
    def __init__(self, data: GeoData, dict_weight: dict, dict_bound: dict):
        """
        data: 
        """
        self.data        = data
        self.dict_weight = dict_weight
        self.dict_bound  = dict_bound

        self.slip        = None
        self.sig_slip    = None
        self.r           = None
        self.dhat        = None
        self.misfit      = None
        self.misfit_gps  = None
        self.misfit_sar  = None
        self.smoothness  = None
        self.ruff        = None
        self.moment      = None
        self.smo_facts   = None
        self.WSAR        = None

    def get_smoothing_factors(self):
        smf1, smf2, step = self.dict_weight['smoothfactor']
        smo_facts = np.arange(smf1, smf2, step)
        return smo_facts if len(smo_facts) > 0 else np.array([smf1])

    def get_weight(self):
        """
        Get weight of GPS and SAR data. relative weight * orignal weight
        """

        # number of data from different SAR datasets (different orbit or satellitte)
        n_sar = self.data.n_sar

        # relative weight between GPS and SAR
        wsar  = self.dict_weight['wsar']
        if self.data.cov_sar.ndim == 1:
            cov_sar = np.diag(self.data.cov_sar)
        elif self.sar.cov_sar.ndim == 2:
            cov_sar = self.data.cov_sar
    
        WSAR  = np.zeros_like(cov_sar)
        if WSAR.size > 0:
            for i in range(len(n_sar)):
                idx0 = sum(n_sar[:i])
                idx1 = sum(n_sar[:i+1])
                WSAR[idx0:idx1,idx0:idx1] = wsar[i]*self.data.W_sar[idx0:idx1,idx0:idx1]
        
        cov_gps = self.data.cov_gps
        WGPS    = np.diag(1/np.sqrt(cov_gps))
        return WGPS, WSAR

    def assemble_data_vector(self, nfault: int):
        """
        Assemble data vector.
        nfault : number of fault
        """
        data      = self.data
        d_gps     = data.d_gps
        d_lev     = data.d_lev
        d_sar     = data.d_sar
        d_geod    = np.hstack([x for x in (d_gps, d_lev) if x.size>0])
        
        len_geod  = len(d_geod)
        len_sar   = len(d_sar)
        
        # weight
        W, WSAR   = self.get_weight()
        d_lap     = np.zeros(3 * nfault)
        
        di_block, dr_block = [], []
        if len_geod > 0:
            wd = W @ d_geod
            di_block.append(wd.flatten())
            dr_block.append(wd.flatten())

        if len_sar > 0:
            ws = WSAR @ d_sar
            di_block.append(ws)
            dr_block.append(ws)

        if d_lap.size > 0:
            di_block.append(d_lap)

        d2I = np.hstack(di_block)
        d2R = np.hstack(dr_block)

        return d2I, d2R

    def solve_least_squares(self, G2I, d2I, bl=None, bu=None):
        if self.dict_bound.get('bound_switch', False) and bl is not None:
            res = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls', lsmr_tol='auto')
            logging.info("Constrained LSQ finished.")
            return res.x
        else:
            x, *_ = np.linalg.lstsq(G2I, d2I, rcond=None)
            logging.info("Unconstrained LSQ finished.")
            return x

    def compute_residuals(self, d2R, G2R, x, W, WSAR, len_geod, len_sar):
        r = np.zeros(len_geod + len_sar)
        dhat = G2R @ x

        mis_gps, mis_sar = 0.0, 0.0

        if len_geod > 0:
            r[:len_geod] = d2R[:len_geod] - W @ dhat[:len_geod]
            mis_gps = r[:len_geod].T @ r[:len_geod]
            logging.info(f"GPS WRSS: {mis_gps:.2f}")

        if len_sar > 0:
            r[len_geod:] = d2R[len_geod:] - WSAR @ dhat[len_geod:]
            mis_sar = r[len_geod:].T @ r[len_geod:]
            logging.info(f"SAR WRSS: {mis_sar:.2f}")

        misfit = mis_gps + mis_sar
        return r, dhat, mis_gps, mis_sar, misfit

    @staticmethod
    def set_slip_bounds(nsegs, ndeps, sar_switch, **varargin):
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
        # if sar_switch == True:
        #     nelems = nelems + 3
            
        # Initialize variables
        # ss_range_c = 0
        ss_range   = np.array([0,0])
        # ds_range_c = 0
        ds_range   = np.array([0,0])
        # op_range_c = 0
        op_range   = np.array([0,0])
        surf_slip  = np.zeros(3*nsegs)
        surf_slip_c  = 0
        s_ss_range   = 0
        # s_ss_range_c = 0
        s_ds_range   = 0
        # s_ds_range_c = 0
        # s_op_range_c = 0
        s_op_range   = 0
        bot_slip     = np.zeros(3*nsegs)
        bot_slip_c   = 0
        # b_ss_range_c = 0 
        # b_ds_range_c = 0 
        # b_op_range_c = 0
        bu           = np.zeros(nelems)
        bl           = np.zeros(nelems)
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
            
        for i in range(nelems):
            # for strike slip component
            if np.remainder(i,3) == 0:
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
            if np.remainder(i,3) == 1:
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
            if np.remainder(i,3) == 2:
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
            bl[0:nelems] = np.genfromtxt(slip_lb, usecols = [7,8,9]).reshape(nelems)
            bu[0:nelems] = np.genfromtxt(slip_ub, usecols = [7,8,9]).reshape(nelems)
            logging.info('The lower bound and upper bound of slip are assgined from {} and {}.'.format(slip_lb, slip_ub))
    
        return bu, bl