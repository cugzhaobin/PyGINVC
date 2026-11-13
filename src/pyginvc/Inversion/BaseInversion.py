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
        """
        Parse ``smoothfactor`` into an array of smoothing factors.

        Accepted forms:
            scalar              -> [value]
            [a]                 -> [a]
            [a, b]              -> explicit factors [a, b]
            [a, b, step]        -> a, a+step, ..., b (end point included when
                                   it lies on the grid)

        A three-element list keeps the historical "start, stop, step" meaning,
        but the stop value is no longer silently dropped.
        """
        raw = self.dict_weight.get('smoothfactor', [0.1, 0.5, 0.1])
        if raw is None:
            raise ValueError('smoothfactor must not be empty.')
        if np.isscalar(raw):
            values = [float(raw)]
        else:
            values = [float(v) for v in np.asarray(raw).ravel()]
        if not values:
            raise ValueError('smoothfactor must contain at least one value.')

        if len(values) == 3:
            start, stop, step = values
            if not np.isfinite([start, stop, step]).all():
                raise ValueError(
                    'smoothfactor must be finite, got %r' % (values,))
            if step == 0:
                raise ValueError('smoothfactor step must be non-zero.')
            if (stop - start) * step < 0:
                raise ValueError(
                    'smoothfactor step (%r) points away from stop (%r) '
                    'starting at start (%r).' % (step, stop, start))
            n = int(np.floor((stop - start) / step + 1e-9)) + 1
            smo_facts = start + step * np.arange(max(n, 1))
        else:
            smo_facts = np.asarray(values, dtype=float)

        if not np.isfinite(smo_facts).all() or (smo_facts < 0).any():
            raise ValueError(
                'smoothfactor values must be finite and non-negative, got %r'
                % (smo_facts.tolist(),))
        return smo_facts

    def get_weight(self):
        """
        Get weight of GPS and SAR data. relative weight * original weight

        Weights are built from the stored covariance arrays.  Older versions
        read ``data.W_sar``, but that attribute is no longer created by
        :class:`GeoData`, which made any InSAR inversion fail with
        ``AttributeError``.
        """

        # InSAR: per-observation variance (the diagonal when a full
        # covariance matrix was supplied) and per-file relative weights.
        cov_sar = np.asarray(getattr(self.data, 'cov_sar', np.array([])))
        n_sar   = np.asarray(getattr(self.data, 'n_sar', 0),
                             dtype=int).reshape(-1)
        var_sar = cov_sar if cov_sar.ndim == 1 else np.diag(cov_sar)
        var_sar = np.asarray(var_sar, dtype=float).ravel()

        if var_sar.size == 0:
            WSAR = np.zeros((0, 0))
        else:
            if n_sar.sum() != var_sar.size:
                raise ValueError(
                    'cov_sar has %d entries but n_sar sums to %d.'
                    % (var_sar.size, n_sar.sum()))
            w_sar = 1.0 / np.sqrt(np.where(var_sar > 0, var_sar, 1.0))

            wsar = np.asarray(self.dict_weight.get('wsar', 1.0),
                              dtype=float).ravel()
            if wsar.size == 1:
                wsar = np.full(n_sar.size, wsar[0])
            if wsar.size != n_sar.size:
                raise ValueError(
                    'wsar has %d values but %d InSAR file(s) are present.'
                    % (wsar.size, n_sar.size))

            offset = 0
            for weight, n in zip(wsar, n_sar):
                w_sar[offset:offset + n] *= weight
                offset += n
            WSAR = np.diag(w_sar)

        # GPS: same diagonal-variance treatment as before.
        cov_gps = np.asarray(getattr(self.data, 'cov_gps', np.array([])))
        var_gps = cov_gps if cov_gps.ndim == 1 else np.diag(cov_gps)
        var_gps = np.asarray(var_gps, dtype=float).ravel()
        if var_gps.size == 0:
            WGPS = np.zeros((0, 0))
        else:
            WGPS = np.diag(1.0 / np.sqrt(np.where(var_gps > 0, var_gps, 1.0)))
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

        # The weighted observation block is the same for the inversion and
        # the residual vector; build it once instead of twice.  Only d2I gets
        # the zero block matching the smoothing rows; d2R stays at the data
        # length because callers use it directly against G2R @ x.
        parts = []
        if len_geod > 0:
            parts.append(np.asarray(W @ d_geod).ravel())
        if len_sar > 0:
            parts.append(np.asarray(WSAR @ d_sar).ravel())
        d2R = np.hstack(parts) if parts else np.array([])

        d2I = (np.hstack([d2R, d_lap]) if d_lap.size > 0
               else d2R.copy())
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

    def compute_residuals(self, d2R, G2R, x, W, WSAR, len_geod, len_sar,
                          dhat=None):
        r = np.zeros(len_geod + len_sar)
        # Callers that already formed the prediction can pass it in, avoiding
        # a duplicate (and potentially large) matrix-vector product.
        if dhat is None:
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
                if i < nsegs*3:
                    if surf_slip_c == 1:
                         bu[i] = surf_slip[i]+s_ss_range
                         bl[i] = surf_slip[i]-s_ss_range
                         if bu[i] > ss_range[1]:
                             bu[i] = ss_range[1]
                         if bl[i] < ss_range[0]:
                             bl[i] = ss_range[0]
                    else:
                         bu[i] = ss_range[1]
                         bl[i] = ss_range[0]
                elif i>nsegs*(ndeps-1)*3 and i<=nsegs*ndeps*3:
                     if (bot_slip_c == 1):
                         bu[i] = bot_slip[i-nsegs*(ndeps-1)*3] + b_ss_range
                         bl[i] = bot_slip[i-nsegs*(ndeps-1)*3] - b_ss_range
                         if bu[i] > ss_range[1]:
                             bu[i] = ss_range[1]
                         if bl[i] < ss_range[0]:
                             bl[i] = ss_range[0]
                     else:
                         bu[i] = ss_range[1]
                         bl[i] = ss_range[0]
                else:
                     bu[i] = ss_range[1]
                     bl[i] = ss_range[0]
                
            # for dip slip component
            if np.remainder(i,3) == 1:
                 if i < nsegs*3:
                     if surf_slip_c == 1:
                         bu[i] = surf_slip[i]+s_ds_range
                         bl[i] = surf_slip[i]-s_ds_range
                         if bu[i] > ds_range[1]:
                             bu[i] = ds_range[1]
                         if bl[i] < ds_range[0]:
                             bl[i] = ds_range[0]
                     else:
                         bu[i] = ds_range[1]
                         bl[i] = ds_range[0]
                 elif i>nsegs*(ndeps-1)*3 and i<=nsegs*ndeps*3:
                     if (bot_slip_c == 1):
                         bu[i] = bot_slip[i-nsegs*(ndeps-1)*3] + b_ds_range
                         bl[i] = bot_slip[i-nsegs*(ndeps-1)*3] - b_ds_range
                         if bu[i] > ds_range[1]:
                             bu[i] = ds_range[1]
                         if bl[i] < ds_range[0]:
                             bl[i] = ds_range[0]
                     else:
                         bu[i] = ds_range[1]
                         bl[i] = ds_range[0]
                 else:
                     bu[i] = ds_range[1]
                     bl[i] = ds_range[0]           
            # openning component uses op_range, s_op_range, b_op_range   
            if np.remainder(i,3) == 2:
                 if i < nsegs*3:
                     if surf_slip_c == 1:
                         bu[i] = surf_slip[i]+s_op_range
                         bl[i] = surf_slip[i]-s_op_range
                         if bu[i] > op_range[1]:
                             bu[i] = op_range[1]
                         if bl[i] < op_range[0]:
                             bl[i] = op_range[0]
                     else:
                         bu[i] = op_range[1]
                         bl[i] = op_range[0]
                 elif i>nsegs*(ndeps-1)*3 and i<=nsegs*ndeps*3:
                     if (bot_slip_c == 1):
                         bu[i] = bot_slip[i-nsegs*(ndeps-1)*3] + b_op_range
                         bl[i] = bot_slip[i-nsegs*(ndeps-1)*3] - b_op_range
                         if bu[i] > op_range[1]:
                             bu[i] = op_range[1]
                         if bl[i] < op_range[0]:
                             bl[i] = op_range[0]
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
        has_lb = bool(str(slip_lb or '').strip())
        has_ub = bool(str(slip_ub or '').strip())
        if has_lb != has_ub:
            raise ValueError(
                'slip_lb and slip_ub must be given together; got %r and %r.'
                % (slip_lb, slip_ub))
        if has_lb and has_ub:
            nelems = 3*nsegs*ndeps
            lb = np.genfromtxt(slip_lb, usecols=[7, 8, 9]).reshape(-1)
            ub = np.genfromtxt(slip_ub, usecols=[7, 8, 9]).reshape(-1)
            if lb.size != nelems or ub.size != nelems:
                raise ValueError(
                    'slip_lb/slip_ub must contain %d slip values (%d patches '
                    'x 3 components); got %d and %d.'
                    % (nelems, nsegs*ndeps, lb.size, ub.size))
            bl[0:nelems] = lb
            bu[0:nelems] = ub
            logging.info('The lower bound and upper bound of slip are assgined from {} and {}.'.format(slip_lb, slip_ub))
    
        return bu, bl
