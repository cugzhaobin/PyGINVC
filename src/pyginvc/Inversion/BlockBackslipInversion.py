
import logging
import numpy as np
from scipy import optimize
from scipy import linalg

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class BlockBackslipInversion:
    def __init__(self, data, fault, green, lap, dict_weight, dict_bound):
        self.fault       = fault
        self.data        = data
        self.green       = green
        self.lap         = lap
        self.dict_weight = dict_weight
        self.dict_bound  = dict_bound
        
        if fault is None or lap is None or green is None:
            self.is_invert_backslip = False
            logging.info('Fault, Green function or Laplacian is not provided. Backslip inversion will be skipped.')
        else:
            self.is_invert_backslip = True


    def assemble_data_vector(self, nfault):
        data = self.data
        # weighted data vector
        d_gps         = data.d_gps
        d_lev         = data.d_lev
        d_sar         = data.d_sar
        d_geod        = np.hstack([x for x in (d_gps, d_lev) if x.size>0])
        
        len_geod      = len(d_geod)
        len_sar       = len(d_sar)
        wgps          = self.dict_weight['wgps']
        wsar          = self.dict_weight['wsar']
        
        # weight
        WSAR          = np.zeros_like(data.W_sar)
        W             = data.W
        
        
        d_lap         = np.zeros(3 * nfault)
        
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

    def assemble_desigin_matrix(self, nblock, G_gps_block, G_sar_block, smo_fact):
        """
        Assemble degin matrix for inversion.
        """
        if self.is_invert_backslip is False:
            nparam = nblock*3
        else:
            nparam = nblock*3 + self.fault.nf*3
            
        if self.lap is None:
            G_lap_part = np.empty((0,nparam))
        else:
            G_lap_part = np.zeros((self.lap.G_lap.shape[0], nparam))
        
        G_gps_part = np.zeros((G_gps_block.shape[0], nparam))
        G_gps_part[:, 0:G_gps_block.shape[1]] = G_gps_block
        
        # GPS part
        if self.is_invert_backslip:
            G_gps_part[:, 3*nblock:3*nblock+self.green.G.shape[0]] = -1*self.green.G
        
        # SAR part
        if self.green.G_sar > 0 and self.green is not None:
            G_sar_part = np.zeros((self.green.G_sar.shape[0], nparam))
            G_sar_part[:, :G_sar_block.shape[1]] = G_sar_block
            G_sar_part[nblock*3+self.green.G.shape[1]:nblock*3+self.green.G.shape[1]+self.green.G_sar.shape[1], 
                    nblock*3+self.green.G.shape[1]:nblock*3+self.green.G.shape[1]+self.green.G_sar.shape[1]] = -1*self.green.G_sar
        
        WSAR  = np.zeros_like(self.data.W_sar)
        W     = self.data.W
        
        if G_sar_block is not None:
            G2I = np.vstack([W @ G_gps_part, WSAR @ G_sar_part, smo_fact*G_lap_part])
            G2R = np.vstack([W @ G_gps_part, WSAR @ G_sar_part])
            G   = np.vstack([G_gps_part, G_sar_part])
        else:
            G2I = np.vstack([W @ G_gps_part, smo_fact*G_lap_part])
            G2R = np.vstack([W @ G_gps_part])
            G   = np.vstack([G_gps_part])

        return G2I, G2R, G

    def inversion(self, nblock, G_gps_block, G_sar_block):
        if self.fault is None:
            d2I, d2R = self.assemble_data_vector(0)
        else:
            d2I, d2R = self.assemble_data_vector(self.fault.nf)
        
        # unpack data
        smoothfactor  = self.dict_weight['smoothfactor']
        dict_bound    = self.dict_bound

        if self.fault is None:
            nf    = 0
            nsegs = 0
            ndeps = 0
            dis_geom_grid = np.empty(0)
        else:
            nf            = self.fault.nf
            nsegs         = self.fault.nsegs
            ndeps         = self.fault.ndeps
            dis_geom_grid = self.fault.dis_geom_grid

        # dimensions
        len_gps        = len(self.data.d_gps)
        len_geod       = len_gps + len(self.data.d_lev)
        len_sar        = len(self.data.d_sar)
        len_all        = len_geod + len_sar

        # smoothing factors
        smf1, smf2, step = smoothfactor
        smo_facts        = np.arange(smf1, smf2, step)
        nsmooth          = len(smo_facts) if len(smo_facts)>0 else 1


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
        sar_switch     = 0


        # inversion loop
        G_sar_list   = np.empty((0,0))
        if self.green is not None:
            nparam = nblock*3 + nf*3
        else:
            nparam = nblock*3
        for i, smo_fact in enumerate(smo_facts):
            logging.info(f'Inversion {i+1} | smoothing factor = {smo_fact}')
            bl = np.full(nparam, -np.inf)
            bu = np.full(nparam, np.inf)
            bl[3*nblock::3]   = -20.0
            bl[3*nblock+1::3] = 0.0
            bl[3*nblock+2::3] = 0.0
            bu[3*nblock::3]   = 00.0
            bu[3*nblock+1::3] = 1e-3
            bu[3*nblock+2::3] = 1e-3
            

            G2I, G2R, G = self.assemble_desigin_matrix(nblock, G_gps_block, G_sar_block, smo_fact)

            # if constrained linear least square method is used
            if dict_bound['bound_switch']:
                res     = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls', lsmr_tol='auto')
                x       = res.x
                logging.info('Constrained linear least square finished.')
            else:
                x, *_   = np.linalg.lstsq(G2I, d2I, rcond=None)
                logging.info('Unconstrained linear least square finished.')


            if self.fault is not None:
                slip[i]  = x[3*nblock:]
            euler[i] = x[:3*nblock]
            

            WSAR  = np.zeros_like(self.data.W_sar)
            W     = self.data.W
            dhat  = G @ x
            if len_geod > 0:
                r[0:len_geod] = d2R[0:len_geod] -  W @ dhat[0:len_geod]
                r_gps         = r[0:len_geod]
                misfit_gps[i] = r_gps.dot(r_gps)
                logging.info('GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_gps[i]))
                logging.info('GPS: WRSS/(N) %f' %(misfit_gps[i]/len_gps))

            if len_sar > 0:
                r[len_geod:]  = d2R[len_geod:] -  WSAR @ dhat[len_geod:]
                r_sar         = r[len_geod:]
                misfit_sar[i] = r_sar.dot(r_sar)
                logging.info('SAR Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_sar[i]))
                logging.info('SAR: WRSS/(N) %f' %(misfit_sar[i]/len_sar))

            misfit[i]     = r.dot(r)
        
            misfit[i]     = r.dot(r)
            smoothness[i] = np.sum((self.lap @ slip[i])**2)
            ruff[i]       = np.sum((self.lap @ slip[i])**2)
        
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