import logging
import numpy as np
from scipy import optimize
from pyginvc.Inversion.BaseInversion import BaseInversion

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class BlockBackslipInversion(BaseInversion):
    def __init__(self, data, fault, green, lap, dict_weight, dict_bound):
        super().__init__(data, dict_weight, dict_bound)
        self.fault       = fault
        self.green       = green
        self.lap         = lap
            
    def assemble_design_matrix(self, nblock, G_gps_block, smo_fact):
        """
        Assemble desgin matrix for inversion.
        """

        nparam     = nblock*3 + self.fault.nf*3
        G_gps_part, G_sar_part = np.empty(0), np.empty(0)
        if self.green.G.size > 0:
            G_gps_part = np.column_stack([G_gps_block, -self.green.G])
        
        # SAR part
        if self.green.G_sar.size >0:
            G_sar_part = np.zeros((self.green.G_sar.shape[0], nparam))
            G_sar_part[:,3*nblock:] = -self.green.G_sar
        
        G_lap_part = np.zeros((self.lap.G_lap.shape[0], nparam))
        G_lap_part[:,3*nblock:] = self.lap.G_lap

        W, WSAR = self.get_weight()
        
        if G_sar_part.size > 0:
            G2I = np.vstack([W @ G_gps_part, WSAR @ G_sar_part, smo_fact*G_lap_part])
            G2R = np.vstack([W @ G_gps_part, WSAR @ G_sar_part])
            G   = np.vstack([G_gps_part, G_sar_part])
        else:
            G2I = np.vstack([W @ G_gps_part, smo_fact*G_lap_part])
            G2R = np.vstack([W @ G_gps_part])
            G   = np.vstack([G_gps_part])

        return G2I, G2R, G

    def run_inversion(self, nblock, G_gps_block):

        d2I, d2R = self.assemble_data_vector(self.fault.nf)
        
        # unpack data
        dict_bound    = self.dict_bound

        nf            = self.fault.nf
        nsegs         = self.fault.nsegs
        ndeps         = self.fault.ndeps

        # dimensions
        len_gps        = len(self.data.d_gps)
        len_geod       = len_gps + len(self.data.d_lev)
        len_sar        = len(self.data.d_sar)
        len_all        = len_geod + len_sar

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
        euler          = np.zeros((nsmooth,3*nblock))

        nparam         = nblock*3 + nf*3
        bu             = np.full(nparam, np.inf)
        bl             = np.full(nparam,-np.inf)
        slp_bu, slp_bl = self.set_slip_bounds(nsegs, ndeps, 0, **dict_bound)
        bl[3*nblock:]  = slp_bl
        bu[3*nblock:]  = slp_bu
        
        W, WSAR        = self.get_weight()

        # inversion loop
        for i, smo_fact in enumerate(smo_facts):
            logging.info(f'Inversion {i+1} | smoothing factor = {smo_fact}')

            G2I, G2R, G = self.assemble_design_matrix(nblock, G_gps_block, smo_fact)

            # if constrained linear least square method is used
            if dict_bound['bound_switch']:
                res     = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls', lsmr_tol='auto')
                x       = res.x
                logging.info('Constrained linear least square finished.')
            else:
                x, *_   = np.linalg.lstsq(G2I, d2I, rcond=None)
                logging.info('Unconstrained linear least square finished.')

            slip[i]  = x[3*nblock:]
            euler[i] = x[:3*nblock]
            
            
            dhat  = G @ x
            if len_geod > 0:
                r[0:len_geod] = d2R[:len_geod] -  W @ dhat[:len_geod]
                r_gps         = r[:len_geod]
                misfit_gps[i] = r_gps.dot(r_gps)
                logging.info('GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_gps[i]))
                logging.info('GPS: WRSS/(N) %f' %(misfit_gps[i]/len_gps))

            if len_sar > 0:
                r[len_geod:]  = d2R[len_geod:] -  WSAR @ dhat[len_geod:]
                r_sar         = r[len_geod:]
                misfit_sar[i] = r_sar.dot(r_sar)
                logging.info('SAR Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_sar[i]))
                logging.info('SAR: WRSS/(N) %f' %(misfit_sar[i]/len_sar))

            misfit[i]     = misfit_gps[i] + misfit_sar[i]
            lap_slip      = self.lap.G_lap @ slip[i]
            smoothness[i] = np.sum((smo_fact*lap_slip)**2)
            ruff[i]       = np.sum(lap_slip**2)
            shearmodulus  = self.green.modulus
            [Mo, Mw]      = self.fault.moment(slip[i].reshape(-1,3), shear_modulus=shearmodulus)
            moment[i]     = np.array([Mo, Mw])
            gps_mod_slip      = G2R[:len_gps, 3*nblock:] @ slip[i]
            gps_mod_rotation  = G2R[:len_gps, :3*nblock] @ euler[i]
            if len_sar>0:
                sar_mod_slip      = G2R[len_geod:, 3*nblock] @ slip[i]
                sar_mod_rotation  = G2R[len_geod:,:3*nblock] @ x[:3*nblock]
            else:
                sar_mod_slip, sar_mod_rotation = np.array([]), np.array([])
        
        self.ruff       = ruff
        self.misfit_gps = misfit_gps
        self.misfit_sar = misfit_sar
        self.misfit     = misfit
        self.slip       = slip
        self.sig_slip   = sig_slip
        self.r          = r
        self.dhat       = dhat
        self.smo_facts  = smo_facts
        self.smoothness = smoothness
        self.WSAR       = WSAR
        self.euler      = euler
        self.x          = x
        self.gps_mod_slip     = gps_mod_slip
        self.gps_mod_rotation = gps_mod_rotation
        self.sar_mod_slip     = sar_mod_slip
        self.sar_mod_rotation = sar_mod_rotation