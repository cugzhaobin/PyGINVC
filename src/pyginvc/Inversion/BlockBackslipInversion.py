#!/usr/bin/env python

import logging
import numpy as np
from scipy import optimize
from scipy.linalg import block_diag
from pyginvc.Inversion.BaseInversion import BaseInversion
from pyginvc.Inversion.Coupling import Coupling

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class BlockBackslipInversion(BaseInversion):
    """
    This class is designed to invert for Euler poles of rigid blocks along with backslip
    """
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
        misfit         = np.zeros(nsmooth)
        misfit_sar     = np.zeros(nsmooth)
        misfit_gps     = np.zeros(nsmooth)
        smoothness     = np.zeros(nsmooth)
        ruff           = np.zeros(nsmooth)
        moment         = np.zeros((nsmooth,2))
        euler          = np.zeros((nsmooth,3*nblock))

        nparam         = nblock*3 + nf*3

        # Per-smoothing-factor results.  These used to be single arrays that the
        # loop overwrote, so only the LAST smoothing factor's answer survived on
        # the object -- useless for a smoothing search or any per-iteration
        # comparison.  Index with [-1] to get the last one, as before.
        x_all            = np.zeros((nsmooth, nparam))
        r                = np.zeros((nsmooth, len_all))
        dhat_all         = np.zeros((nsmooth, len_all))
        gps_mod_slip     = np.zeros((nsmooth, len_gps))
        gps_mod_rotation = np.zeros((nsmooth, len_gps))
        sar_mod_slip     = np.zeros((nsmooth, len_sar))
        sar_mod_rotation = np.zeros((nsmooth, len_sar))
        status           = np.zeros(nsmooth, dtype=int)
        bu             = np.full(nparam, np.inf)
        bl             = np.full(nparam,-np.inf)
        # Falsy (absent, '' or None) means "leave omega free".  Testing for the
        # key with `in` instead made an empty string reach np.genfromtxt(''),
        # so the only way to turn the constraint off was to delete the key.
        eulerfile = dict_bound.get('euler_range')
        if eulerfile:
            euler_val = np.genfromtxt(eulerfile)
            bl[:3*nblock] = euler_val-0.001*abs(euler_val)
            bu[:3*nblock] = euler_val+0.001*abs(euler_val)
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
                res       = optimize.lsq_linear(G2I, d2I, (bl, bu), method='bvls', lsmr_tol='auto')
                x         = res.x
                status[i] = res.status
                # 1/2 = converged, 3 = the unconstrained solution is optimal
                # (no bound active), which is also a valid optimum.
                if res.status not in (1, 2, 3):
                    logging.warning('BVLS did not converge at smoothing factor %g: '
                                    'status=%d (%s).  status 0 = iteration limit, '
                                    '-1 = no progress; the result below is NOT a '
                                    'converged optimum.',
                                    smo_fact, res.status, getattr(res, 'message', ''))
                else:
                    logging.info('Constrained linear least square finished (status=%d).',
                                 res.status)
            else:
                x, *_     = np.linalg.lstsq(G2I, d2I, rcond=None)
                status[i] = 0
                logging.info('Unconstrained linear least square finished.')

            slip[i]  = x[3*nblock:]
            euler[i] = x[:3*nblock]
            x_all[i] = x

            dhat        = G @ x
            dhat_all[i] = dhat
            if len_geod > 0:
                r[i, 0:len_geod] = d2R[:len_geod] -  W @ dhat[:len_geod]
                r_gps            = r[i, :len_geod]
                misfit_gps[i]    = r_gps.dot(r_gps)
                logging.info('GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_gps[i]))
                logging.info('GPS: WRSS/(N) %f' %(misfit_gps[i]/len_gps))

            if len_sar > 0:
                r[i, len_geod:] = d2R[len_geod:] -  WSAR @ dhat[len_geod:]
                r_sar           = r[i, len_geod:]
                misfit_sar[i]   = r_sar.dot(r_sar)
                logging.info('SAR Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit_sar[i]))
                logging.info('SAR: WRSS/(N) %f' %(misfit_sar[i]/len_sar))

            misfit[i]     = misfit_gps[i] + misfit_sar[i]
            lap_slip      = self.lap.G_lap @ slip[i]
            smoothness[i] = np.sum((smo_fact*lap_slip)**2)
            ruff[i]       = np.sum(lap_slip**2)
            shearmodulus  = self.green.modulus
            [Mo, Mw]      = self.fault.moment(slip[i].reshape(-1,3), shear_modulus=shearmodulus)
            moment[i]     = np.array([Mo, Mw])
            # decompose the model with the unweighted G (G2R is W @ G, so using it
            # here would return the weighted model instead of a displacement in mm/yr)
            gps_mod_slip[i]     = G[:len_gps, 3*nblock:] @ slip[i]
            gps_mod_rotation[i] = G[:len_gps, :3*nblock] @ euler[i]
            if len_sar>0:
                sar_mod_slip[i]     = G[len_geod:, 3*nblock:] @ slip[i]
                sar_mod_rotation[i] = G[len_geod:, :3*nblock] @ euler[i]

        self.ruff       = ruff
        self.misfit_gps = misfit_gps
        self.misfit_sar = misfit_sar
        self.misfit     = misfit
        self.slip       = slip
        self.sig_slip   = sig_slip
        self.r          = r
        self.dhat       = dhat_all
        self.smo_facts  = smo_facts
        self.smoothness = smoothness
        self.WSAR       = WSAR
        self.euler      = euler
        self.x          = x_all
        self.status     = status
        self.gps_mod_slip     = gps_mod_slip
        self.gps_mod_rotation = gps_mod_rotation
        self.sar_mod_slip     = sar_mod_slip
        self.sar_mod_rotation = sar_mod_rotation


class LinkedBlockBackslipInversion(BaseInversion):
    """
    Block-backslip inversion with the slip deficit tied to the long-term rate.

        v_obs = B(x) omega - G(x) s(x)
        s(x)  = chi(x) * s_lt(x; omega)

    with one chi per patch, bounded by ``dict_coupling.chi_range`` (default
    [0, 1]) and smoothed by the same scalar Laplacian that smooths a free slip
    field.  This is the Meade-style parameterisation: the fault's slip
    direction and long-term rate come from the relative block rotation, and
    the geodetic data only decide how much of that motion is currently locked.

    Why the iteration
    -----------------
    ``s_lt`` depends on omega, and omega is fitted jointly with chi, so the
    forward model is bilinear in (omega, chi).  The fixed point is:

        1. omega <- weighted block-only least squares
        2. s_lt(omega) from the blocks flanking the fault at its SURFACE trace
        3. solve [B, -G diag(s_lt)] [omega; chi] with chi in [chi_lo, chi_hi]
           and a Laplacian penalty on chi
        4. repeat until omega and chi stop moving

    The final solve uses the converged omega, so the returned slip, residual
    and chi all correspond to the same s_lt.

    Units match the rest of the code: B is the dimensionless ``buildG`` matrix,
    omega is whatever units ``buildG`` produces, and G has mm/yr per mm/yr of
    slip, so chi is dimensionless and s_lt is in mm/yr.
    """

    def __init__(self, data, blocks, fault, green, lap,
                 dict_weight, dict_bound, dict_coupling):
        super().__init__(data, dict_weight, dict_bound)
        self.blocks        = blocks
        self.fault         = fault
        self.green         = green
        self.lap           = lap
        self.dict_coupling = dict_coupling or {}

        chi_range         = self.dict_coupling.get('chi_range', [0.0, 1.0])
        self.chi_lo       = float(chi_range[0])
        self.chi_hi       = float(chi_range[1])
        if not self.chi_lo < self.chi_hi:
            raise ValueError('dict_coupling.chi_range must be increasing, got %r'
                             % (chi_range,))
        self.max_iter      = max(int(self.dict_coupling.get('max_iter', 6)), 1)
        self.tol           = float(self.dict_coupling.get('tol', 1.0e-6))
        self.require_exact = bool(self.dict_coupling.get('require_exact', False))
        self.min_slt       = float(self.dict_coupling.get('min_slt', 1.0))

        self.coupling = Coupling(blocks, fault, green=green,
                                 dict_coupling=self.dict_coupling)
        self.coupling.local_basis()

    # ------------------------------------------------------------- matrices
    def _fault_columns(self, Gslip):
        """
        Gslip @ s = -Gslip @ (s_lt * chi), so the chi coefficient is

            G_chi[:, i] = -sum_k Gslip[:, 3*i+k] * s_lt_vec[i, k].

        Written with slices instead of a per-patch loop: Gslip[:, k::3] is
        exactly the k-th component of every patch.
        """
        Gf = np.zeros((Gslip.shape[0], self.fault.nf))
        sv = self.coupling.s_lt_vec
        for k in range(3):
            Gf -= Gslip[:, k::3] * sv[:, k]
        return Gf

    def _chi_laplacian(self):
        """
        Scalar Laplacian on chi, taken from the (3*nf, 3*nf) slip Laplacian.
        The three slip components have identical blocks, so the ss block
        ``G_lap[::3, ::3]`` is the one to use.  Rows/columns of patches with
        undefined s_lt are zeroed: their chi is not a free parameter and must
        not smooth against its neighbours.
        """
        L = np.array(self.lap.G_lap[::3, ::3], dtype=float)
        bad = ~self.coupling.mask
        if bad.any():
            L[bad, :] = 0.0
            L[:, bad] = 0.0
        return L

    def _initial_euler(self, G_gps_block):
        """Weighted block-only least squares, the starting omega."""
        W, _ = self.get_weight()
        d_gps = np.asarray(self.data.d_gps, dtype=float).ravel()
        if G_gps_block.shape[0] != d_gps.size:
            raise ValueError('G_gps_block has %d rows but d_gps has %d entries'
                             % (G_gps_block.shape[0], d_gps.size))
        A = W @ G_gps_block
        b = W @ d_gps
        omega, *_ = np.linalg.lstsq(A, b, rcond=None)
        logging.info('linked model: initial block-only omega fitted from %d GPS '
                     'components', d_gps.size)
        return omega

    # ----------------------------------------------------------------- solve
    def _solve(self, nblock, G_gps_block, smo_fact):
        """Build and solve the linear sub-problem for the current s_lt."""
        nf      = self.fault.nf
        nparam  = 3*nblock + nf
        d_gps   = np.asarray(self.data.d_gps, dtype=float).ravel()
        d_lev   = np.asarray(self.data.d_lev, dtype=float).ravel()
        d_sar   = np.asarray(self.data.d_sar, dtype=float).ravel()
        n_gps, n_lev, n_sar = d_gps.size, d_lev.size, d_sar.size
        n_geod  = n_gps + n_lev

        W, WSAR = self.get_weight()
        if n_lev:
            # get_weight() only covers GPS; leveling is unweighted for now
            W = block_diag(W, np.eye(n_lev))
        d_geod = np.hstack([d_gps, d_lev]) if n_lev else d_gps

        if self.green.G.shape[0] != n_geod:
            raise ValueError('green.G has %d rows, expected %d GPS+leveling rows'
                             % (self.green.G.shape[0], n_geod))

        # block rotation: no vertical component, so zero rows for leveling
        B_geod = np.vstack([G_gps_block,
                            np.zeros((n_lev, 3*nblock))])
        Gf_geod = self._fault_columns(self.green.G)
        G_geod  = np.hstack([B_geod, Gf_geod])

        if n_sar:
            if self.green.G_sar.shape[0] != n_sar:
                raise ValueError('green.G_sar has %d rows, expected %d'
                                 % (self.green.G_sar.shape[0], n_sar))
            G_sar = np.hstack([np.zeros((n_sar, 3*nblock)),
                               self._fault_columns(self.green.G_sar)])
        else:
            G_sar = np.zeros((0, nparam))

        Lchi    = self._chi_laplacian()
        Lpad    = np.zeros((nf, nparam))
        Lpad[:, 3*nblock:] = Lchi

        if n_sar:
            G2I = np.vstack([W @ G_geod, WSAR @ G_sar, smo_fact*Lpad])
            G2R = np.vstack([W @ G_geod, WSAR @ G_sar])
            G   = np.vstack([G_geod, G_sar])
            d2I = np.hstack([W @ d_geod, WSAR @ d_sar, np.zeros(nf)])
            d2R = np.hstack([W @ d_geod, WSAR @ d_sar])
        else:
            G2I = np.vstack([W @ G_geod, smo_fact*Lpad])
            G2R = np.vstack([W @ G_geod])
            G   = np.vstack([G_geod])
            d2I = np.hstack([W @ d_geod, np.zeros(nf)])
            d2R = np.hstack([W @ d_geod])

        bl = np.full(nparam, -np.inf)
        bu = np.full(nparam,  np.inf)
        eulerfile = self.dict_bound.get('euler_range')
        if eulerfile:
            euler_val = np.asarray(np.genfromtxt(eulerfile), dtype=float).ravel()
            if euler_val.size != 3*nblock:
                raise ValueError('euler_range has %d values, expected %d'
                                 % (euler_val.size, 3*nblock))
            bl[:3*nblock] = euler_val - 0.001*np.abs(euler_val)
            bu[:3*nblock] = euler_val + 0.001*np.abs(euler_val)

        # chi bound; patches without a defined s_lt are fixed at 0
        bl[3*nblock:] = self.chi_lo
        bu[3*nblock:] = self.chi_hi
        bad = ~self.coupling.mask
        if bad.any():
            bl[3*nblock:][bad] = 0.0
            bu[3*nblock:][bad] = 0.0

        if self.dict_bound.get('bound_switch', False):
            res    = optimize.lsq_linear(G2I, d2I, (bl, bu),
                                         method='bvls', lsmr_tol='auto')
            x      = res.x
            status = res.status
            # scipy's BVLS status: 1/2 = converged, 3 = the unconstrained
            # solution is optimal (no bound is active, also a valid optimum).
            if res.status not in (1, 2, 3):
                logging.warning('linked BVLS did not converge: status=%d (%s)',
                                res.status, getattr(res, 'message', ''))
        else:
            x, *_  = np.linalg.lstsq(G2I, d2I, rcond=None)
            status = 0

        x[3*nblock:] = np.clip(x[3*nblock:], self.chi_lo, self.chi_hi)
        if bad.any():
            x[3*nblock:][bad] = 0.0

        omega = x[:3*nblock]
        chi   = x[3*nblock:]
        slip  = (self.coupling.s_lt_vec * chi[:, None]).reshape(-1)

        dhat = G @ x
        r    = np.zeros(n_geod + n_sar)
        mis_gps = mis_sar = 0.0
        if n_geod:
            r[:n_geod] = d2R[:n_geod] - W @ dhat[:n_geod]
            mis_gps    = r[:n_geod] @ r[:n_geod]
        if n_sar:
            r[n_geod:] = d2R[n_geod:] - WSAR @ dhat[n_geod:]
            mis_sar    = r[n_geod:] @ r[n_geod:]

        lap_chi = Lchi @ chi
        return dict(x=x, omega=omega, chi=chi, slip=slip, r=r, dhat=dhat,
                    G=G, G_geod=G_geod, d2R=d2R,
                    misfit_gps=mis_gps, misfit_sar=mis_sar,
                    misfit=mis_gps + mis_sar,
                    smoothness=np.sum((smo_fact*lap_chi)**2),
                    ruff=np.sum(lap_chi**2), status=status,
                    n_geod=n_geod, len_all=n_geod + n_sar)

    # ------------------------------------------------------------ main loop
    def run_inversion(self, nblock, G_gps_block):
        nf          = self.fault.nf
        nparam      = 3*nblock + nf
        smo_facts   = self.get_smoothing_factors()
        nsmooth     = len(smo_facts)
        omega0      = self._initial_euler(G_gps_block)

        slip           = np.zeros((nsmooth, 3*nf))
        chi_all        = np.zeros((nsmooth, nf))
        euler          = np.zeros((nsmooth, 3*nblock))
        x_all          = np.zeros((nsmooth, nparam))
        r              = np.zeros((nsmooth, len(self.data.d_gps)
                                   + len(self.data.d_lev)
                                   + len(self.data.d_sar)))
        dhat_all       = np.zeros_like(r)
        gps_mod_slip   = np.zeros((nsmooth, len(self.data.d_gps)))
        gps_mod_rot    = np.zeros((nsmooth, len(self.data.d_gps)))
        sar_mod_slip   = np.zeros((nsmooth, len(self.data.d_sar)))
        sar_mod_rot    = np.zeros((nsmooth, len(self.data.d_sar)))
        misfit         = np.zeros(nsmooth)
        misfit_gps     = np.zeros(nsmooth)
        misfit_sar     = np.zeros(nsmooth)
        smoothness     = np.zeros(nsmooth)
        ruff           = np.zeros(nsmooth)
        moment         = np.zeros((nsmooth, 2))
        status         = np.zeros(nsmooth, dtype=int)
        n_iter_done    = np.zeros(nsmooth, dtype=int)

        for i, smo_fact in enumerate(smo_facts):
            logging.info('Linked inversion %d | smoothing factor = %g',
                         i + 1, smo_fact)
            omega     = omega0.copy()
            chi_prev  = None

            for it in range(self.max_iter):
                self.coupling.long_term_slip(
                    omega, require_exact=self.require_exact,
                    min_slt=self.min_slt, verbose=False)
                out       = self._solve(nblock, G_gps_block, smo_fact)
                omega_new = out['omega']
                d_omega   = (np.linalg.norm(omega_new - omega)
                             / max(np.linalg.norm(omega_new), 1e-12))
                d_chi     = (np.inf if chi_prev is None else
                             np.linalg.norm(out['chi'] - chi_prev)
                             / max(np.linalg.norm(out['chi']), 1e-12))
                logging.info('   iter %d: d_omega = %.3e   d_chi = %.3e',
                             it + 1, d_omega, d_chi)
                omega    = omega_new
                chi_prev = out['chi']
                n_iter_done[i] = it + 1
                if it > 0 and d_omega < self.tol and d_chi < self.tol:
                    break

            # final solve with the converged omega
            self.coupling.long_term_slip(
                omega, require_exact=self.require_exact,
                min_slt=self.min_slt, verbose=True)
            out = self._solve(nblock, G_gps_block, smo_fact)

            x_all[i]  = out['x']
            slip[i]   = out['slip']
            chi_all[i] = out['chi']
            euler[i]  = out['omega']
            r[i]      = out['r']
            dhat_all[i] = out['dhat']
            misfit[i]     = out['misfit']
            misfit_gps[i] = out['misfit_gps']
            misfit_sar[i] = out['misfit_sar']
            smoothness[i] = out['smoothness']
            ruff[i]       = out['ruff']
            status[i]     = out['status']
            moment[i]     = np.array(self.fault.moment(
                out['slip'].reshape(-1, 3), shear_modulus=self.green.modulus))

            n_geod = out['n_geod']
            gps_mod_slip[i] = out['G'][:len(self.data.d_gps),
                                       3*nblock:] @ out['chi']
            gps_mod_rot[i]  = out['G'][:len(self.data.d_gps),
                                        :3*nblock] @ out['omega']
            if len(self.data.d_sar):
                sar_mod_slip[i] = out['G'][n_geod:, 3*nblock:] @ out['chi']
                sar_mod_rot[i]  = out['G'][n_geod:, :3*nblock] @ out['omega']

            logging.info('GPS WRSS %10.3f   WRSS/N %8.3f', misfit_gps[i],
                         misfit_gps[i] / max(len(self.data.d_gps), 1))
            if len(self.data.d_sar):
                logging.info('SAR WRSS %10.3f   WRSS/N %8.3f', misfit_sar[i],
                             misfit_sar[i] / max(len(self.data.d_sar), 1))

        # the driver reads sol.coupling for s_lt / chi / save()
        self.coupling.set_chi(chi_all[-1])
        self.coupling.slip = slip[-1].reshape(-1, 3)

        self.ruff            = ruff
        self.misfit_gps      = misfit_gps
        self.misfit_sar      = misfit_sar
        self.misfit          = misfit
        self.slip            = slip
        self.chi             = chi_all
        self.sig_slip        = np.zeros_like(slip)
        self.r               = r
        self.dhat            = dhat_all
        self.smo_facts       = smo_facts
        self.smoothness      = smoothness
        self.moment          = moment
        self.WSAR            = self.get_weight()[1]
        self.euler           = euler
        self.x               = x_all
        self.status          = status
        self.n_iter          = n_iter_done
        self.gps_mod_slip    = gps_mod_slip
        self.gps_mod_rotation = gps_mod_rot
        self.sar_mod_slip    = sar_mod_slip
        self.sar_mod_rotation = sar_mod_rot

        self.report_result()

    def report_result(self):
        nf = self.fault.nf
        chi = self.chi[-1]
        m = self.coupling.mask
        if not m.any():
            logging.warning('linked model: no patch has a defined s_lt')
            return
        c = chi[m]
        on = int(((c <= self.chi_lo + 1e-9) | (c >= self.chi_hi - 1e-9)).sum())
        logging.info('linked chi over %d patches: min %.3f  median %.3f  max %.3f',
                     c.size, c.min(), np.median(c), c.max())
        logging.info('   at a bound: %d (%.0f%%)   zero: %d   one: %d',
                     on, 100.0*on/c.size,
                     int((c <= self.chi_lo + 1e-9).sum()),
                     int((c >= self.chi_hi - 1e-9).sum()))
        logging.info('   iterations per smoothing factor: %s',
                     ', '.join(str(int(v)) for v in self.n_iter))
