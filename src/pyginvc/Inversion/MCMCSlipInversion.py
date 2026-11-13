#!/usr/bin/env python
"""

"""

import os
import time
import logging
import numpy as np
from scipy import linalg
from scipy import sparse as sp
from scipy.sparse import linalg as splinalg
from pyginvc.Export.MultiOutput import MultiOutput

try:  # emcee / corner 
    import emcee
    _HAS_EMCEE = True
except Exception:  
    _HAS_EMCEE = False

try:
    import corner
    _HAS_CORNER = True
except Exception:  # pragma: no cover
    corner = None
    _HAS_CORNER = False

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

class MCMCSlipInversion(object):
    '''
    Bayesian/MCMC inversion of fault slip distribution。

    Input:
        inversion: instance of MultiFaultInversion
        config: mcmc part of the yaml file, e.g. ``config['mcmc']``
        export: export part of the yaml file
    '''

    def __init__(self, inversion, config=None, export=None):
        self.inv    = inversion
        self.cfg    = config or {}
        self.export = export or {}

        out_cfg     = self.cfg.get('output') or {}
        base        = (out_cfg.get('outdir') or self.export.get('outpath', './'))
        subdir      = out_cfg.get('subdir', 'mcmc')
        self.outdir = (base if os.path.isabs(subdir)
                       else os.path.join(base, subdir))
        os.makedirs(self.outdir, exist_ok=True)

        self.prepared = False
        self.result = None

    # ------------------------------------------------------------------ #
    # prepare data
    # ------------------------------------------------------------------ #
    def _ensure_prepared(self):
        if not self.prepared:
            self.prepare()

    def prepare(self):
        """Prepare for MCMC inversion: assemble weighted data, design matrix, reference model, and basis."""
        inv = self.inv
        inv.prepare()

        cfg    = self.cfg
        bcfg   = cfg.get('basis') or {}
        pcfg   = cfg.get('prior') or {}
        scfg   = cfg.get('sampler') or {}
        outcfg = cfg.get('output') or {}

        self.nf3    = inv.nf3
        self.n_geod = inv.n_geod
        self.n_sar  = inv.n_sar
        self.n_data = self.n_geod + self.n_sar
        if self.n_data == 0:
            raise ValueError('No observation data for MCMC inversion.')

        # data
        G_parts, d_parts = [], []
        if self.n_geod:
            G_parts.append(inv._Gw)
            d_parts.append(inv._d_geodw)
        if self.n_sar:
            G_parts.append(inv._Gsarw)
            d_parts.append(inv._d_sarw)
        self.G_w   = np.vstack(G_parts)
        self.d_w   = np.hstack(d_parts)
        self.d_raw = (np.hstack([inv.d_geod, inv.d_sar])
                      if self.n_data else np.array([]))

        # ---- ramp
        self.n_ramp     = int(inv.n_ramp) if inv.n_ramp else 0
        self.ramp_keep  = np.zeros(self.n_ramp, dtype=bool)
        self.G_ramp_w   = np.zeros((self.n_data, 0))
        self.G_ramp_raw = np.zeros((self.n_data, 0))
        self.ramp_sigma = np.zeros(0)
        if self.n_ramp > 0:
            norms = np.linalg.norm(inv.G_ramp, axis=0)
            self.ramp_keep = norms > 0
            nr = int(self.ramp_keep.sum())
            if nr > 0:
                self.G_ramp_w = inv._G_ramp_w[:, self.ramp_keep]
                self.G_ramp_raw = inv.G_ramp[:, self.ramp_keep]
                sigma_full = self._ramp_prior_sigma()
                self.ramp_sigma = np.maximum(
                    sigma_full[self.ramp_keep], 1e-12)
        self.nr = self.G_ramp_w.shape[1]
        self.G_ramp_raw_sar = (self.G_ramp_raw[self.n_geod:]
                               if self.n_sar else np.zeros((0, self.nr)))

        # ---- smoothing factor ---------------------------------------- #
        self._init_smoothing_prior(pcfg.get('smoothing') or {})

        # ---------------------------------------- #
        ref_input = bcfg.get('reference', 'auto')
        zero_ref = (ref_input is False
                    or (isinstance(ref_input, str)
                        and ref_input.lower() in ('zero', 'none', 'false')))
        if isinstance(ref_input, (int, float)) and not isinstance(ref_input, bool):
            ref_smo = float(ref_input)
        else:
            ref_smo = (bcfg.get('ref_smoothing')
                       or self._default_ref_smoothing())
        self.ref_smoothing = float(ref_smo)

        det = inv.solve(self.ref_smoothing)
        s_det = np.asarray(det['slip'], dtype=float).reshape(-1)
        self.slip_scale = max(float(np.max(np.abs(s_det))), 1.0)
        rho0_full = np.asarray(det['ramp'], dtype=float).reshape(-1)
        if self.nr and rho0_full.size:
            self.rho0 = rho0_full[self.ramp_keep]
        else:
            self.rho0 = np.zeros(self.nr)

        # 边界
        self.bound_switch = bool(getattr(inv, 'bound_switch', False))
        if self.bound_switch and inv.bl is not None and inv.bu is not None:
            self.bl = np.asarray(inv.bl, dtype=float)[:self.nf3].copy()
            self.bu = np.asarray(inv.bu, dtype=float)[:self.nf3].copy()
        else:
            self.bl = np.full(self.nf3, -np.inf)
            self.bu = np.full(self.nf3, np.inf)
        self.bound_margin = float((pcfg.get('bound_margin')
                                   if pcfg.get('bound_margin') is not None
                                   else 1e-3))
        self.bound_tol = float(pcfg.get('bound_tol', 1e-6))

        if zero_ref:
            s_ref = np.zeros(self.nf3)
        else:
            s_ref = self._interior(s_det)
        if not self._slip_feasible(s_ref):
            logging.warning('Reference slip is infeasible; using a feasible '
                            'mid-point reference instead.')
            s_ref = self._feasible_reference()
        self.s_ref = s_ref

        self.b0 = self.G_w @ self.s_ref
        # 采样时 d_w - b0 不变，只算一次，避免每次 _log_prob 重复分配。
        self.res0 = self.d_w - self.b0
        self.L_sref = inv.G_lap @ self.s_ref

        # ---- 降维基 -------------------------------------------------- #
        self._build_basis(bcfg)

        # ---- 先验尺度 ------------------------------------------------- #
        self.slip_coeff_sigma = pcfg.get('slip_coeff_sigma')
        if self.slip_coeff_sigma is not None:
            self.slip_coeff_sigma = float(self.slip_coeff_sigma)

        # ---- 采样参数 ------------------------------------------------ #
        self.sampler_method = str(scfg.get('method', 'emcee')).lower()
        if self.sampler_method == 'emcee' and not _HAS_EMCEE:
            logging.warning('emcee is not installed; falling back to the '
                            'built-in Metropolis sampler.')
            self.sampler_method = 'metropolis'
        self.nwalkers = int(scfg.get('nwalkers', 0) or 0)
        self.moves_name = str(scfg.get('moves', 'de')).lower()
        self.burnin = int(scfg.get('burnin', 500))
        self.nsteps = int(scfg.get('nsteps', 2000))
        self.thin = max(1, int(scfg.get('thin', 1)))
        self.seed = scfg.get('seed')
        self.progress = bool(scfg.get('progress', False))
        if self.progress:
            import importlib.util
            if importlib.util.find_spec('tqdm') is None:
                logging.warning('tqdm is not installed; disabling emcee '
                                'progress bar.')
                self.progress = False
        self.metro_cfg = scfg.get('metropolis') or {}

        self.nsamples = int(outcfg.get('nsamples', 500))
        self.save_samples = bool(outcfg.get('save_samples', True))
        self.save_chain = bool(outcfg.get('save_chain', True))
        self.make_plots = bool(outcfg.get('plots', True))
        self.make_corner = bool(outcfg.get('corner', True))

        self.prepared = True
        logging.info('MCMC prepared: n_data=%d, nf3=%d, k=%d, ramp_params=%d, '
                     'lambda in [%.4g, %.4g].',
                     self.n_data, self.nf3, self.k, self.nr,
                     self.lam_lower, self.lam_upper)
        return self

    # ------------------------------------------------------------------ #
    # 先验与参考状态
    # ------------------------------------------------------------------ #
    def _init_smoothing_prior(self, cfg):
        dist = str(cfg.get('dist', 'loguniform')).lower()
        if dist not in ('loguniform', 'uniform', 'fixed'):
            raise ValueError("mcmc.prior.smoothing.dist must be "
                             "'loguniform', 'uniform' or 'fixed'.")
        self.smoothing_dist = dist

        sf = self.inv.config.get('smoothfactor', [0.1, 0.5, 0.1])
        sf = np.asarray(sf, dtype=float).reshape(-1)
        lo0 = float(np.min(sf[[0, min(1, len(sf) - 1)]]))
        hi0 = float(np.max(sf[[0, min(1, len(sf) - 1)]]))
        custom_lo = cfg.get('lower')
        custom_hi = cfg.get('upper')
        lo = float(custom_lo) if custom_lo is not None else max(lo0 / 10., 1e-6)
        hi = float(custom_hi) if custom_hi is not None else max(hi0 * 10., 1e-2)
        if hi <= lo:
            hi = max(lo * 10., lo + 1e-3)
        self.lam_lower = lo
        self.lam_upper = hi
        self.lam_fixed = float(cfg.get('value', np.sqrt(lo * hi)))
        if self.lam_fixed <= 0:
            raise ValueError('Fixed smoothing factor must be positive.')
        self.sample_smoothing = dist != 'fixed'

        # 高斯先验 N(0, alpha^2 (L^T L)^+) 的归一化因子 alpha^(-M)：
        # alpha = 1/lambda，所以对数后验多出 M*log(lambda)。
        # 这是 Fukuda & Johnson (2008) 式 (14)/(16) 中的 (alpha^2)^(-M/2)
        # 项；关掉它等价于只使用“条件式”平滑先验。
        self.smooth_norm = bool(cfg.get('normalization', True))
        n_eff = cfg.get('n_effective')
        if self.smooth_norm:
            self.smooth_power = (int(n_eff) if n_eff is not None
                                 else self._smoothing_rank())
        else:
            self.smooth_power = 0

    def _default_ref_smoothing(self):
        if not self.sample_smoothing:
            return self.lam_fixed
        return float(np.sqrt(self.lam_lower * self.lam_upper))

    def _smoothing_rank(self):
        """L 的秩（高斯平滑先验的有效约束数）。"""
        L = np.asarray(self.inv.G_lap)
        n = L.shape[0]
        if n == 0:
            return 0
        if n <= 1800:
            s = np.linalg.svd(L, compute_uv=False)
            if s.size == 0 or s[0] <= 0:
                return 0
            return int(np.sum(s > s[0] * 1e-10))
        # 大问题避免 O(n^3)：逐断层 3 分量 Laplacian 的零空间为常数模式
        return max(1, n - 3 * len(list(self.inv.faults.pairs())))

    def _ramp_prior_sigma(self):
        """按 ramp 列构造弱高斯先验尺度（数据单位 / 数据单位每 km）。"""
        rcfg = (self.cfg.get('prior') or {}).get('ramp') or {}
        data = self.inv.data
        origin = self.inv.origin

        # 数据幅值与空间跨度，用于自动定先验尺度
        dvals = []
        if self.n_sar:
            dvals.append(np.asarray(data.d_sar, dtype=float).reshape(-1))
        if self.n_geod:
            dvals.append(np.asarray(self.inv.d_geod, dtype=float).reshape(-1))
        dscale = (max(float(np.max(np.abs(v))) for v in dvals if v.size)
                  if any(v.size for v in dvals) else 1.0)
        if not np.isfinite(dscale) or dscale <= 0:
            dscale = 1.0

        xy_all = []
        if self.n_geod and len(data.llh_gps):
            from pyginvc.libs import geotools as gt
            xy_all.append(np.array([gt.llh2utm(p, origin)
                                    for p in data.llh_gps]).reshape(-1, 2))
        if self.n_sar and len(data.llh_sar):
            from pyginvc.libs import geotools as gt
            xy_all.append(np.array([gt.llh2utm(p, origin)
                                    for p in data.llh_sar]).reshape(-1, 2))
        if xy_all:
            xy = np.vstack(xy_all)
            span = float(np.max(np.ptp(xy, axis=0))) if len(xy) > 1 else 0.0
        else:
            span = 0.0
        span = max(span, 1.0)  # km

        offset = rcfg.get('offset_sigma')
        offset = float(offset) if offset is not None else 2.0 * dscale
        if offset <= 0:
            offset = 1.0
        grad = rcfg.get('gradient_sigma')
        grad = float(grad) if grad is not None else offset / span
        if grad <= 0:
            grad = 1e-6
        gps_off = rcfg.get('gps_offset_sigma')
        gps_off = float(gps_off) if gps_off is not None else offset
        gps_grad = rcfg.get('gps_gradient_sigma')
        gps_grad = float(gps_grad) if gps_grad is not None else grad
        gps_rot = rcfg.get('gps_rotation_sigma')
        gps_rot = (float(gps_rot) if gps_rot is not None
                   else offset / 6370.0)

        sig = []
        for gps, g in data.gps_pairs():
            g = g or {}
            if not g.get('ramp', False) or len(gps.llh_gps) == 0:
                continue
            method = int(g.get('ramp_method', 1))
            if method == 2:
                sig.extend([gps_rot] * 3)
            else:
                sig.extend([gps_off, gps_off, gps_grad])
        for sar, s in data.sar_pairs():
            s = s or {}
            if len(sar.llh_sar) == 0:
                continue
            # MultiFaultInversion.build_ramp 对每个非空 track 都估计常数
            # 参考电平；ramp=True 时再加 [east, north] 坡度项。
            if s.get('ramp', False):
                # make_green_sar_ramp = [1, x, y]：常数偏移 + 东西/南北坡度
                sig.extend([grad, grad, grad])
            else:
                sig.extend([offset])
        return np.asarray(sig, dtype=float)

    def _interior(self, s):
        """把参考模型向边界内部推一个很小的 margin。"""
        s = np.asarray(s, dtype=float).copy()
        lb_f = np.isfinite(self.bl)
        ub_f = np.isfinite(self.bu)
        s[lb_f] = np.maximum(s[lb_f], self.bl[lb_f] + self.bound_margin)
        s[ub_f] = np.minimum(s[ub_f], self.bu[ub_f] - self.bound_margin)
        return s

    def _feasible_reference(self):
        mid = np.zeros(self.nf3)
        lb_f = np.isfinite(self.bl)
        ub_f = np.isfinite(self.bu)
        both = lb_f & ub_f
        mid[both] = 0.5 * (self.bl[both] + self.bu[both])
        mid[lb_f & ~ub_f] = self.bl[lb_f & ~ub_f] + self.bound_margin
        mid[ub_f & ~lb_f] = self.bu[ub_f & ~lb_f] - self.bound_margin
        return mid

    def _slip_feasible(self, s):
        s = np.asarray(s)
        return bool(np.all(s >= self.bl - self.bound_tol)
                    and np.all(s <= self.bu + self.bound_tol))

    # ------------------------------------------------------------------ #
    # 降维基
    # ------------------------------------------------------------------ #
    def _build_basis(self, cfg):
        method = str(cfg.get('method', 'laplacian')).lower()
        nf3 = self.nf3
        max_modes = int(cfg.get('max_modes', 200) or 200)
        max_modes = max(2, min(max_modes, nf3))
        nmodes = int(cfg.get('nmodes', 0) or 0)

        if method == 'none':
            V_all = np.eye(nf3)
            self.basis_method = 'none'
        elif method == 'laplacian':
            V_all = self._lap_eigvecs(max_modes)
            self.basis_method = 'laplacian'
        elif method in ('svd', 'data-svd'):
            V_all = self._svd_vecs(max_modes)
            self.basis_method = 'svd'
        elif method in ('posterior', 'pca'):
            lam0 = float(cfg.get('lambda0', self.ref_smoothing))
            V_all = self._posterior_eigvecs(max_modes, lam0)
            self.basis_method = 'posterior'
        else:
            raise ValueError("mcmc.basis.method must be one of "
                             "'laplacian', 'svd', 'posterior', 'none'.")

        if method == 'none':
            k = nf3 if nmodes <= 0 else min(nmodes, nf3)
            V = V_all[:, :k]
            energy_frac = np.ones(k)
        else:
            Gv = self.G_w @ V_all
            energy = np.sum(Gv ** 2, axis=0)
            total = float(np.sum(self.G_w ** 2))
            frac = np.cumsum(energy) / total if total > 0 else np.zeros_like(energy)
            if nmodes > 0:
                k = min(nmodes, V_all.shape[1])
            else:
                thresh = float(cfg.get('energy', 0.99))
                k = int(np.searchsorted(frac, thresh) + 1)
                k = max(k, int(cfg.get('min_modes', 5)))
                k = min(k, V_all.shape[1])
            V = V_all[:, :k]
            energy_frac = frac[:k]

        self.V = np.ascontiguousarray(V)
        self.k = int(self.V.shape[1])
        self.V_all_modes = int(V_all.shape[1])
        self.basis_energy = energy_frac
        self.Gv = self.G_w @ self.V                  # (n_data, k)
        self.LV = self.inv.G_lap @ self.V            # (nf3, k)

        if method != 'none' and nmodes <= 0:
            logging.info('Reduced model: %d -> %d modes (data energy %.3f).',
                         nf3, self.k, float(energy_frac[-1]))
        if self.k >= nf3:
            logging.info('MCMC uses the full model space (%d parameters).',
                         self.k)

    def _lap_eigvecs(self, kmax):
        L = np.asarray(self.inv.G_lap)
        n = L.shape[0]
        if n <= 1:
            return np.eye(n)
        dense_limit = int((self.cfg.get('basis') or {}).get(
            'dense_limit', 1800))
        if n <= dense_limit:
            kmax = max(1, min(kmax, n))
            A = L.T @ L
            _, vecs = linalg.eigh(A, subset_by_index=[0, kmax - 1])
            return vecs
        kmax = max(1, min(kmax, n - 1))
        Ls = sp.csr_matrix(L)
        A = Ls.T @ Ls
        vals, vecs = splinalg.eigsh(A, k=kmax, which='SM')
        order = np.argsort(vals)
        return vecs[:, order]

    def _posterior_eigvecs(self, kmax, lam0):
        L = np.asarray(self.inv.G_lap)
        n = L.shape[0]
        if n <= 1:
            return np.eye(n)
        dense_limit = int((self.cfg.get('basis') or {}).get(
            'dense_limit', 1800))
        if n <= dense_limit:
            kmax = max(1, min(kmax, n))
            A = self.G_w.T @ self.G_w + (lam0 ** 2) * (L.T @ L)
            _, vecs = linalg.eigh(A, subset_by_index=[n - kmax, n - 1])
            return vecs[:, ::-1]
        kmax = max(1, min(kmax, n - 1))
        Ls = sp.csr_matrix(L)
        A = sp.csr_matrix(self.G_w.T @ self.G_w) + (lam0 ** 2) * (Ls.T @ Ls)
        vals, vecs = splinalg.eigsh(A, k=kmax, which='LA')
        order = np.argsort(vals)[::-1]
        return vecs[:, order]

    def _svd_vecs(self, kmax):
        G = self.G_w
        n_data, n = G.shape
        m = min(n_data, n)
        if m <= 1800 or n_data * n <= 4e7:
            kmax = max(1, min(kmax, m))
            _, _, Vt = np.linalg.svd(G, full_matrices=False)
            return Vt[:kmax].T
        kmax = max(1, min(kmax, m - 1))
        _, s, Vt = splinalg.svds(G, k=kmax, which='LM')
        order = np.argsort(s)[::-1]
        return Vt[order].T

    # ------------------------------------------------------------------ #
    # 对数后验（支持 emcee 的 vectorize=True）
    # ------------------------------------------------------------------ #
    def _log_prob(self, theta):
        theta = np.atleast_2d(np.asarray(theta, dtype=float))
        nw = theta.shape[0]
        k, nr = self.k, self.nr
        coeff = theta[:, :k]
        rho = theta[:, k:k + nr] if nr else None

        if self.sample_smoothing:
            loglam = theta[:, -1]
            lam = np.exp(np.clip(loglam, -50.0, 50.0))
            inside = ((loglam >= np.log(self.lam_lower))
                      & (loglam <= np.log(self.lam_upper)))
        else:
            loglam = None
            lam = np.full(nw, self.lam_fixed)
            inside = np.ones(nw, dtype=bool)

        s = self.s_ref[None, :] + coeff @ self.V.T
        res = self.res0[None, :] - coeff @ self.Gv.T
        if nr:
            res = res - rho @ self.G_ramp_w.T
        f_data = np.einsum('ij,ij->i', res, res)

        ls = self.L_sref[None, :] + coeff @ self.LV.T
        f_smooth = lam ** 2 * np.einsum('ij,ij->i', ls, ls)

        logp = -0.5 * (f_data + f_smooth)

        if nr:
            logp -= 0.5 * np.sum((rho / self.ramp_sigma[None, :]) ** 2, axis=1)
        if self.slip_coeff_sigma:
            logp -= 0.5 * np.sum(
                (coeff / float(self.slip_coeff_sigma)) ** 2, axis=1)

        # 平滑先验的归一化因子 alpha^(-M) = lambda^M
        if self.sample_smoothing and self.smooth_power:
            logp = logp + self.smooth_power * loglam

        feasible = np.all((s >= self.bl - self.bound_tol)
                          & (s <= self.bu + self.bound_tol), axis=1)
        logp = np.where(feasible & inside, logp, -np.inf)

        # 以 log(lambda) 为采样坐标时，lambda 上的均匀先验需要 Jacobian
        if self.sample_smoothing and self.smoothing_dist == 'uniform':
            logp = logp + loglam
        return logp

    # ------------------------------------------------------------------ #
    # 初始状态
    # ------------------------------------------------------------------ #
    def _base_theta(self):
        lam0 = self.lam_fixed if not self.sample_smoothing else \
            float(np.clip(self.ref_smoothing, self.lam_lower, self.lam_upper))
        parts = [np.zeros(self.k), self.rho0]
        if self.sample_smoothing:
            parts.append(np.array([np.log(lam0)]))
        return np.hstack(parts) if any(p.size for p in parts) else np.zeros(0)

    def _theta_scales(self):
        k, nr = self.k, self.nr
        scatter = float((self.cfg.get('basis') or {}).get(
            'init_scatter', 0.05))
        a_sigma = scatter * self.slip_scale / np.sqrt(max(k, 1))
        parts = [np.full(k, max(a_sigma, 1e-12))]
        if nr:
            parts.append(np.maximum(self.ramp_sigma, 1e-12) * scatter)
        if self.sample_smoothing:
            parts.append(np.array([0.05]))
        return np.hstack(parts) if any(p.size for p in parts) else np.zeros(0)

    def _feasible_theta(self, theta):
        coeff = np.asarray(theta)[:self.k]
        s = self.s_ref + self.V @ coeff
        return self._slip_feasible(s) and (
            not self.sample_smoothing
            or (self.lam_lower <= np.exp(theta[-1]) <= self.lam_upper))

    def _initial_positions(self, n, rng):
        base = self._base_theta()
        scales = self._theta_scales()
        ndim = base.size
        pos = np.empty((n, ndim))
        for i in range(n):
            scale = scales.copy()
            for _ in range(30):
                cand = base + rng.uniform(-1.0, 1.0, ndim) * scale
                if self._feasible_theta(cand):
                    pos[i] = cand
                    break
                scale *= 0.5
            else:
                cand = base.copy()
                if self.sample_smoothing:
                    cand[-1] += rng.normal(0.0, 1e-4)
                pos[i] = cand
        return pos

    # ------------------------------------------------------------------ #
    # 采样
    # ------------------------------------------------------------------ #
    def run(self):
        self._ensure_prepared()
        rng = np.random.default_rng(self.seed)
        t0 = time.time()

        if self.sampler_method == 'emcee':
            chain, logp, acc, tau, trace = self._run_emcee(rng)
        else:
            chain, logp, acc, tau, trace = self._run_metropolis(rng)

        self.elapsed = time.time() - t0
        logging.info('MCMC finished in %.1f s: %d samples, acceptance %.2f.',
                     self.elapsed, len(chain), acc)

        result = self._postprocess(chain, logp, rng)
        result.update({
            'chain': chain,
            'log_prob': logp,
            'acceptance': float(acc),
            'autocorr_time': tau,
            'trace': trace,
            'elapsed': self.elapsed,
            'names': self._parameter_names(),
        })
        self.result = result
        return result

    def _emcee_moves(self):
        name = self.moves_name
        if name in ('de', 'differential'):
            return emcee.moves.DEMove()
        if name in ('desnooker', 'snooker'):
            return emcee.moves.DESnookerMove()
        if name in ('stretch', 'default'):
            return emcee.moves.StretchMove()
        logging.warning("Unknown mcmc.sampler.moves=%r; using DEMove.", name)
        return emcee.moves.DEMove()

    def _run_emcee(self, rng):
        ndim = self._base_theta().size
        nwalkers = self.nwalkers or max(2 * ndim, 32)
        nwalkers = max(nwalkers, 2 * ndim)
        if nwalkers % 2:
            nwalkers += 1
        p0 = self._initial_positions(nwalkers, rng)

        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, self._log_prob, moves=self._emcee_moves(),
            vectorize=True)
        niter = self.burnin + self.nsteps
        sampler.run_mcmc(p0, niter, progress=self.progress)
        chain = sampler.get_chain(discard=self.burnin, thin=self.thin,
                                  flat=True)
        logp = sampler.get_log_prob(discard=self.burnin, thin=self.thin,
                                    flat=True).reshape(-1)
        acc = float(np.mean(sampler.acceptance_fraction))
        tau = None
        try:
            tau = np.asarray(sampler.get_autocorr_time(tol=0), dtype=float)
        except Exception:
            pass

        trace = {}
        prod = sampler.get_chain(discard=self.burnin, thin=self.thin)
        if self.sample_smoothing:
            trace['log_lambda'] = prod[:, :, -1]
        if self.nr:
            trace['ramp'] = prod[:, :, self.k:self.k + self.nr]
        return chain, logp, acc, tau, trace

    def _run_metropolis(self, rng):
        cfg = self.metro_cfg
        nchains = int(cfg.get('nchains', 4))
        adapt_every = max(1, int(cfg.get('adapt_every', 100)))
        target = float(cfg.get('target_accept', 0.25))
        initial_step = float(cfg.get('initial_step', 0.1))

        ndim = self._base_theta().size
        scales = self._theta_scales()
        widths = np.maximum(scales * initial_step, 1e-12)

        all_chain, all_logp = [], []
        trace_loglam, trace_ramp = [], []
        total_acc = 0.0
        total_trials = 0

        for _c in range(nchains):
            theta = self._initial_positions(1, rng)[0]
            lp = float(self._log_prob(theta[None, :])[0])
            n_total = self.burnin + self.nsteps
            acc_j = np.zeros(ndim)
            trial_j = np.zeros(ndim)

            for it in range(n_total):
                # 逐坐标更新（Fukuda & Johnson 式）：
                # 一次只改一个参数，因此在 s >= lb 这类边界附近不会被
                # 其它坐标的合法更新“连坐”，接受率显著高于整向量提议。
                for j in range(ndim):
                    prop = theta.copy()
                    prop[j] += rng.uniform(-1.0, 1.0) * widths[j]
                    lp_prop = float(self._log_prob(prop[None, :])[0])
                    trial_j[j] += 1.0
                    total_trials += 1
                    if np.isfinite(lp_prop) and (
                            np.log(rng.random()) < lp_prop - lp):
                        theta, lp = prop, lp_prop
                        acc_j[j] += 1.0
                        total_acc += 1.0

                # 只在 burn-in 内按坐标自适应；生产阶段固定提议宽度
                if (it + 1) % adapt_every == 0 and it < self.burnin:
                    rate = acc_j / np.maximum(trial_j, 1.0)
                    factor = np.ones(ndim)
                    factor[rate > target * 1.2] = 1.2
                    factor[rate < target * 0.8] = 0.8
                    widths = np.clip(widths * factor, 1e-14, 1e6)
                    acc_j[:] = 0.0
                    trial_j[:] = 0.0

                if it >= self.burnin and (it - self.burnin) % self.thin == 0:
                    all_chain.append(theta.copy())
                    all_logp.append(lp)
                    if self.sample_smoothing:
                        trace_loglam.append([theta[-1]])
                    if self.nr:
                        trace_ramp.append(theta[self.k:self.k + self.nr][None, :])

        chain = np.asarray(all_chain)
        logp = np.asarray(all_logp)
        acc = total_acc / max(total_trials, 1)
        trace = {}
        if self.sample_smoothing:
            trace['log_lambda'] = np.asarray(trace_loglam)
        if self.nr:
            trace['ramp'] = np.concatenate(trace_ramp, axis=0)
        return chain, logp, acc, None, trace

    # ------------------------------------------------------------------ #
    # 后处理
    # ------------------------------------------------------------------ #
    def _parameter_names(self):
        names = ['a_%03d' % i for i in range(self.k)]
        names += ['ramp_%02d' % i for i in range(self.nr)]
        if self.sample_smoothing:
            names += ['log_lambda']
        return names

    def _postprocess(self, chain, logp, rng):
        k, nr = self.k, self.nr
        ntot = len(chain)
        if ntot == 0:
            raise RuntimeError('MCMC produced no post-burn-in samples; '
                               'increase nsteps or reduce burnin.')

        coeff = chain[:, :k]
        rho = chain[:, k:k + nr] if nr else np.zeros((ntot, 0))
        if self.sample_smoothing:
            lam = np.exp(chain[:, -1])
        else:
            lam = np.full(ntot, self.lam_fixed)

        slip_all = self.s_ref[None, :] + coeff @ self.V.T
        ramp_sar_all = (self.G_ramp_raw_sar @ rho.T).T if (
            self.n_sar and nr) else np.zeros((ntot, self.n_sar))

        ns = min(max(int(self.nsamples), 2), ntot)
        idx = (rng.choice(ntot, size=ns, replace=False)
               if ns < ntot else np.arange(ntot))
        slip_s = slip_all[idx]
        ramp_s = ramp_sar_all[idx]
        lam_s = lam[idx]
        rho_s = rho[idx]

        slip_mean = slip_all.mean(axis=0)
        slip_std = slip_all.std(axis=0)
        slip_p16, slip_p50, slip_p84 = np.percentile(
            slip_s, [16, 50, 84], axis=0)
        rho_mean = rho.mean(axis=0) if nr else np.zeros(0)
        ramp_sar_mean = (self.G_ramp_raw_sar @ rho_mean
                         if (self.n_sar and nr) else np.zeros(self.n_sar))

        # 数据拟合（后验均值）
        fit_geod = self.inv.Gmat @ slip_mean
        fit_sar = self.inv.G_sar @ slip_mean
        if nr:
            fit_geod = fit_geod + self.inv.G_ramp[:self.n_geod] @ \
                self._full_ramp(rho_mean)
            fit_sar = fit_sar + self.inv.G_ramp[self.n_geod:] @ \
                self._full_ramp(rho_mean)
        dhat = np.hstack([fit_geod, fit_sar])
        r = (self.d_raw - dhat) if self.n_data else np.array([])
        res_geod = self.inv.d_geod - fit_geod
        res_sar = self.inv.d_sar - fit_sar
        misfit_gps = (float(np.sum((self.inv.W @ res_geod) ** 2))
                      if self.n_geod else 0.0)
        misfit_sar = (float(np.sum((self.inv.WSAR @ res_sar) ** 2))
                      if self.n_sar else 0.0)
        misfit = misfit_gps + misfit_sar

        # 矩的后验分布
        Mo_s = np.empty(ns)
        Mw_s = np.empty(ns)
        default_modulus = self.inv.config.get('modulus')
        if default_modulus is None:
            default_modulus = ((self.inv.config.get('green') or {})
                               .get('modulus', 3e10))
        for j, m in enumerate(slip_s):
            Mo_s[j], Mw_s[j] = self.inv.faults.total_moment(
                m, default_modulus=default_modulus)

        lam_p16, lam_p50, lam_p84 = np.percentile(lam_s, [16, 50, 84])
        return {
            'slip_samples': slip_s,
            'slip_mean': slip_mean,
            'slip_std': slip_std,
            'slip_p16': slip_p16,
            'slip_p50': slip_p50,
            'slip_p84': slip_p84,
            'ramp_samples': ramp_s,
            'ramp_sar_mean': ramp_sar_mean,
            'rho_samples': rho_s,
            'rho_mean': rho_mean,
            'rho_std': rho.std(axis=0) if nr else np.zeros(0),
            'lambda_samples': lam_s,
            'lambda_mean': float(lam.mean()),
            'lambda_median': float(lam_p50),
            'lambda_p16': float(lam_p16),
            'lambda_p84': float(lam_p84),
            'dhat': dhat,
            'r': r,
            'misfit_gps': misfit_gps,
            'misfit_sar': misfit_sar,
            'misfit': misfit,
            'Mo_samples': Mo_s,
            'Mw_samples': Mw_s,
            'Mo': float(np.mean(Mo_s)),
            'Mw': float(np.mean(Mw_s)),
        }

    def _full_ramp(self, rho_keep):
        """把剔除零列后的 ramp 系数放回原始 n_ramp 长度。"""
        full = np.zeros(self.n_ramp)
        if self.n_ramp:
            full[self.ramp_keep] = rho_keep
        return full

    # ------------------------------------------------------------------ #
    # 输出
    # ------------------------------------------------------------------ #
    def write_outputs(self, result=None):
        self._ensure_prepared()
        result = result if result is not None else self.result
        if result is None:
            raise RuntimeError('Run the sampler before writing outputs.')

        os.makedirs(self.outdir, exist_ok=True)
        written = []
        h5_mcmc = os.path.join(
            self.outdir,
            (self.cfg.get('output') or {}).get('chain_file',
                                               'mcmc_chain.h5'))
        self._write_chain(h5_mcmc, result)
        written.append(h5_mcmc)

        # 后验均值解 -> 复用 MultiOutput 的 HDF5 / 数据 / 模型写出
        sol = self._solution_view(result)
        output = MultiOutput(self.inv.data, self.inv.faults, sol=sol,
                             archdir=self.outdir, outdir=self.outdir)
        # WriteHDF5 joins its ``filename`` onto output.outdir itself, so pass
        # the bare name; handing it a pre-joined path duplicated the subdir.
        h5_sol = output.WriteHDF5(filename='mcmc_solutions.h5')
        written.append(h5_sol)

        written += output.WriteData(dhat=result['dhat'], outdir=self.outdir,
                                    ramp=result['ramp_sar_mean']
                                    if self.nr and self.n_sar else None)
        for name, key in (('mcmc_mean', 'slip_mean'),
                          ('mcmc_p16', 'slip_p16'),
                          ('mcmc_p50', 'slip_p50'),
                          ('mcmc_p84', 'slip_p84'),
                          ('mcmc_std', 'slip_std')):
            written += output.WriteSlipModel(
                slip=result[key], k=0, outdir=self.outdir, prefix=name)
        written += output.WriteGMTSlip(slip=result['slip_mean'],
                                       outdir=self.outdir)

        summary = self._write_summary(result)
        written.append(summary)

        if self.make_plots:
            try:
                written += self._make_plots(result)
            except Exception as exc:  # 绘图失败不应毁掉采样结果
                logging.warning('MCMC plot failed: %s', exc)

        self.written = written
        logging.info('%d MCMC output file(s) written to %s.',
                     len(written), self.outdir)
        return written

    def _solution_view(self, result):
        sol = _SolutionView()
        sol.origin = self.inv.origin
        sol.slip = np.asarray(result['slip_mean']).reshape(-1)
        sol.r = np.asarray(result['r']).reshape(-1)
        sol.misfit = np.asarray([result['misfit']])
        sol.dhat = np.asarray(result['dhat']).reshape(-1)
        sol.ramp_sar = (result['ramp_sar_mean']
                        if self.nr and self.n_sar else None)
        sol.sar_ramp_flag = self.inv.sar_ramp_flag
        sol.smo_facts = np.asarray([result['lambda_median']])
        sol.misfit_gps = np.asarray([result['misfit_gps']])
        sol.misfit_sar = np.asarray([result['misfit_sar']])
        sol.moment = np.asarray([[result['Mo'], result['Mw']]])
        return sol

    def _write_chain(self, path, result):
        import h5py
        chain = result['chain']
        names = np.array(result['names'], dtype='S40')
        with h5py.File(path, 'w') as f:
            g = f.create_group('mcmc')
            if self.save_chain:
                g.create_dataset('chain', data=chain.astype('f4'),
                                 compression='gzip')
                g.create_dataset('log_prob',
                                 data=np.asarray(result['log_prob'],
                                                 dtype='f8'))
            g.create_dataset('parameter_names', data=names)
            g.create_dataset('basis', data=self.V.astype('f8'),
                             compression='gzip')
            g.create_dataset('slip_reference', data=self.s_ref)
            g.create_dataset('lambda', data=result['lambda_samples'])
            g.create_dataset('ramp_coefficients',
                             data=result['rho_samples'])
            g.create_dataset('ramp_keep', data=self.ramp_keep)
            g.create_dataset('slip_mean', data=result['slip_mean'])
            g.create_dataset('slip_std', data=result['slip_std'])
            g.create_dataset('slip_p16', data=result['slip_p16'])
            g.create_dataset('slip_p50', data=result['slip_p50'])
            g.create_dataset('slip_p84', data=result['slip_p84'])
            g.create_dataset('moment_Mo', data=result['Mo_samples'])
            g.create_dataset('moment_Mw', data=result['Mw_samples'])
            if self.save_samples:
                g.create_dataset('slip_samples',
                                 data=result['slip_samples'].astype('f4'),
                                 compression='gzip')
                g.create_dataset('ramp_samples',
                                 data=result['ramp_samples'].astype('f4'),
                                 compression='gzip')
            g.attrs['basis_method'] = self.basis_method
            g.attrs['nmodes'] = self.k
            g.attrs['n_fault_params'] = self.nf3
            g.attrs['n_data'] = self.n_data
            g.attrs['n_ramp_params'] = self.nr
            g.attrs['sampler'] = self.sampler_method
            g.attrs['acceptance'] = result['acceptance']
            g.attrs['burnin'] = self.burnin
            g.attrs['nsteps'] = self.nsteps
            g.attrs['thin'] = self.thin
            g.attrs['smoothing_dist'] = self.smoothing_dist
            g.attrs['smoothing_lower'] = self.lam_lower
            g.attrs['smoothing_upper'] = self.lam_upper
            g.attrs['reference_smoothing'] = self.ref_smoothing
            g.attrs['elapsed_seconds'] = result['elapsed']
            if result['autocorr_time'] is not None:
                g.attrs['autocorr_time'] = result['autocorr_time']
        logging.info('MCMC chain is written to %s.' % path)
        return path

    def _write_summary(self, result):
        path = os.path.join(self.outdir, 'mcmc_summary.txt')
        tau = result['autocorr_time']
        lines = [
            'MCMCSlipInversion summary',
            '=========================',
            'mode               : %s' % self.sampler_method,
            'elapsed (s)        : %.2f' % result['elapsed'],
            'posterior samples  : %d (chain %d, thin %d)'
            % (len(result['chain']), self.nsteps, self.thin),
            'acceptance         : %.3f' % result['acceptance'],
            'data / model / ramp: %d / %d / %d'
            % (self.n_data, self.nf3, self.nr),
            'basis              : %s, k = %d (%d modes available)'
            % (self.basis_method, self.k, self.V_all_modes),
            'reference smoothing: %.6g' % self.ref_smoothing,
            'smoothing prior    : %s [%.6g, %.6g]'
            % (self.smoothing_dist, self.lam_lower, self.lam_upper),
            'lambda (16/50/84)  : %.6g  %.6g  %.6g'
            % (result['lambda_p16'], result['lambda_median'],
               result['lambda_p84']),
            'lambda mean/std    : %.6g  %.6g'
            % (result['lambda_mean'], np.std(result['lambda_samples'])),
            'misfit (gps/sar)   : %.6g  %.6g'
            % (result['misfit_gps'], result['misfit_sar']),
            'Mo mean / Mw mean  : %.6e  %.4f' % (result['Mo'], result['Mw']),
        ]
        if result['Mw_samples'].size:
            lines.append('Mw (16/50/84)      : %.4f  %.4f  %.4f'
                         % tuple(np.percentile(result['Mw_samples'],
                                               [16, 50, 84])))
        if self.nr:
            lines.append('ramp coefficients (posterior mean +/- std):')
            for i in range(self.nr):
                lines.append('  ramp_%02d  % .6e  % .6e'
                             % (i, result['rho_mean'][i],
                                result['rho_std'][i]))
        if tau is not None:
            lines.append('autocorr time      : %s'
                         % np.array2string(np.asarray(tau), precision=2))
            try:
                ess = len(result['chain']) / np.max(np.asarray(tau))
                lines.append('min effective size : %.1f' % ess)
            except Exception:
                pass
        text = '\n'.join(lines) + '\n'
        with open(path, 'w') as fid:
            fid.write(text)
        logging.info('\n' + text)
        return path

    def _make_plots(self, result):
        import matplotlib.pyplot as plt
        written = []
        trace = result.get('trace') or {}

        # lambda trace
        if trace.get('log_lambda') is not None:
            arr = np.asarray(trace['log_lambda'])
            fig, ax = plt.subplots(figsize=(8, 3.2))
            if arr.ndim == 1:
                ax.plot(arr, lw=0.8)
            else:
                ax.plot(arr, lw=0.6, alpha=0.6)
            ax.set_xlabel('iteration')
            ax.set_ylabel(r'$\log\lambda$')
            ax.set_title('Smoothing factor trace')
            fig.tight_layout()
            path = os.path.join(self.outdir, 'mcmc_lambda_trace.pdf')
            fig.savefig(path)
            plt.close(fig)
            written.append(path)

        # ramp trace
        if trace.get('ramp') is not None and self.nr:
            arr = np.asarray(trace['ramp'])
            fig, axes = plt.subplots(min(3, self.nr), 1,
                                     figsize=(8, 2.2 * min(3, self.nr)),
                                     squeeze=False)
            for i in range(min(3, self.nr)):
                ax = axes[i, 0]
                y = arr[:, :, i] if arr.ndim == 3 else arr[:, i]
                ax.plot(y, lw=0.5, alpha=0.6)
                ax.set_ylabel('ramp_%02d' % i)
            axes[-1, 0].set_xlabel('iteration')
            fig.tight_layout()
            path = os.path.join(self.outdir, 'mcmc_ramp_trace.pdf')
            fig.savefig(path)
            plt.close(fig)
            written.append(path)

        # corner: lambda + ramp + 前几个滑动系数
        if self.make_corner and _HAS_CORNER:
            chain = result['chain']
            cols, labels = [], []
            for i in range(min(6, self.k)):
                cols.append(i)
                labels.append(r'$a_{%d}$' % i)
            for i in range(self.nr):
                cols.append(self.k + i)
                labels.append(r'$\rho_{%d}$' % i)
            if self.sample_smoothing:
                cols.append(chain.shape[1] - 1)
                labels.append(r'$\log\lambda$')
            if len(cols) >= 2:
                samples = chain[:, cols]
                fig = corner.corner(samples, labels=labels,
                                    show_titles=True, quantiles=[0.16, 0.5, 0.84])
                path = os.path.join(self.outdir, 'mcmc_corner.pdf')
                fig.savefig(path)
                plt.close(fig)
                written.append(path)
        return written
