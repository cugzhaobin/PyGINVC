#!/usr/bin/env python
# Written by Zhao Bin, when he was at UC Berkeley. Feb 14 2016
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
# Mod by Zhao Bin, Dec. 7, 2018. We use HDF5 to store solutions

import logging
import numpy as np
from scipy import optimize
from scipy import linalg
from pyginvc.Geometry.Patch import Fault
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Inversion.BaseInversion import BaseInversion


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%M-%Y %H:%M:%S")


class GeoInversion(BaseInversion):
    '''
    GeoInversion represents a class of slip inversion using geodetic data
    '''
    
    def __init__(self, flt: Fault, data: GeoData, green, lap, dict_weight, dict_bound):
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

    def run_inversion_clean(self):
        """Invert for slip distribution"""

        # unpack data
        d_gps         = self.data.d_gps
        d_lev         = self.data.d_lev
        d_sar         = self.data.d_sar
        G_lap         = self.lap.G_lap
        dict_bound    = self.dict_bound

        nf            = self.flt.nf
        nsegs         = self.flt.nsegs
        ndeps         = self.flt.ndeps

        d2I, d2R      = self.assemble_data_vector(nf)

        # data weight
        W, WSAR       = self.get_weight()
        
        # dimensions
        len_geod       = len(d_gps) + len(d_lev)
        len_sar        = len(d_sar)
        len_all        = len_geod + len_sar
        n_ramp         = self.green.G_gps_ramp.shape[1]+self.green.G_sar_ramp.shape[1]

        # smoothing factors
        smo_facts      = self.get_smoothing_factors()
        nsmooth        = len(smo_facts)

        # init parameters
        slip           = np.zeros((nsmooth, 3*nf))
        sig_slip       = np.zeros((nsmooth, 3*nf))
        r              = np.zeros(len_all)
        misfit         = np.zeros(nsmooth)
        misfit_sar     = np.zeros(nsmooth)
        misfit_gps     = np.zeros(nsmooth)
        smoothness     = np.zeros(nsmooth)
        ruff           = np.zeros(nsmooth)
        moment         = np.zeros((nsmooth,2))
        ramp           = np.zeros((nsmooth, n_ramp))
        sar_switch     = 0

        if dict_bound['bound_switch']:
            bu, bl = self.set_slip_bounds(nsegs, ndeps, sar_switch, **dict_bound)

            if n_ramp>0:
                bl = np.hstack([bl, np.full(n_ramp, -np.inf)])
                bu = np.hstack([bu, np.full(n_ramp,  np.inf)])

            logging.info('Boundaries for constrained linear is constructed.')

        # inversion loop
        shearmodulus = self.green.modulus

        for i, smo_fact in enumerate(smo_facts):
            logging.info(f'Inversion {i+1} | smoothing factor = {smo_fact}')
            
            G2I, G2R = self.assemble_design_matrix(smo_fact)
            # if constrained linear least square method is used
            if dict_bound['bound_switch']:
                res     = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls', lsmr_tol='auto')
                x       = res.x
                logging.info('Constrained linear least square finished.')
            else:
                x, *_   = np.linalg.lstsq(G2I, d2I, rcond=None)
                logging.info('Unconstrained linear least square finished.')

            slip[i] = x[:3*nf]
            ramp[i] = x[3*nf:]
            
            dhat    = G2R @ x

            if len_geod > 0:
                r[0:len_geod] = d2R[0:len_geod] - W @ dhat[0:len_geod]
                r_gps         = r[0:len_geod]
                misfit_gps[i] = r_gps.dot(r_gps)
                logging.info('GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_gps[i]))
                logging.info('GPS: WRSS/(N) %f' %(misfit_gps[i]/len(d_gps)))

            if len_sar > 0:
                r[len_geod:]  = d2R[len_geod:] - WSAR @ dhat[len_geod:]
                r_sar         = r[len_geod:]
                misfit_sar[i] = r_sar.dot(r_sar)
                logging.info('SAR Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_sar[i]))
                logging.info('SAR: WRSS/(N) %f' %(misfit_sar[i]/len_sar))

            misfit[i]     = r.dot(r)
            # ||smo*L@s||^2, matching the penalty rows smo*G_lap in G2I.
            lap_slip      = G_lap @ slip[i]
            ruff[i]       = float(lap_slip.dot(lap_slip))
            smoothness[i] = smo_fact ** 2 * ruff[i]

            # compute the moment
            Mo, Mw     = self.flt.moment(slip[i].reshape(-1,3), shear_modulus=shearmodulus)
            moment[i]  = np.array([Mo, Mw])

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
            
    def assemble_design_matrix(self, smo_fact):
        """
        Assemble desgin matrix for inversion.
        """
        G              = self.green.G
        G_sar          = self.green.G_sar
        G_gps_ramp     = self.green.G_gps_ramp
        G_sar_ramp     = self.green.G_sar_ramp
        n_ramp         = G_gps_ramp.shape[1]+G_sar_ramp.shape[1]

        # scale the G_lap
        G_laps       = smo_fact*self.lap.G_lap
        G_blocks1    = []
        G_blocks2    = []
        ramp_blocks1 = []
        ramp_blocks2 = []

        W, WSAR      = self.get_weight()
        if  G.size > 0:
            G_blocks1.append(W @ G)
            G_blocks2.append(G)
            ramp_blocks1.append(W @ G_gps_ramp)
            ramp_blocks2.append(G_gps_ramp)

        if G_sar.size > 0:
            G_blocks1.append(WSAR @ G_sar)
            G_blocks2.append(G_sar)
            ramp_blocks1.append(WSAR @ G_sar_ramp)
            ramp_blocks2.append(G_sar_ramp)

        G_blocks1.append(G_laps)
        G2I   = np.vstack(G_blocks1)
        Gramp = np.vstack([linalg.block_diag(*ramp_blocks1),
                            np.zeros((G_laps.shape[0], n_ramp))])
        G2I   = np.column_stack((G2I, Gramp))

        G2R   = np.column_stack((np.vstack(G_blocks2),
                                    linalg.block_diag(*ramp_blocks2)))
        return G2I, G2R

        
if __name__ == '__main__':
    print('inverse')
