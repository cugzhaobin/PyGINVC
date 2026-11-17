#!/usr/bin/env python
"""Slip inversion for rectangular, triangular, or mixed fault models.

The inversion is geometry agnostic once each fault has a Green's function
block: every patch contributes three columns ``(strike-slip, dip-slip,
opening)`` and the full model is the horizontal concatenation of those
blocks.  Rectangular patches use Okada/EDCMP Green's functions and
``RectPlaneLap`` smoothing; triangular patches use Nikkhoo/Meade
Green's functions and ``TriLap`` smoothing.  The normal equations never need
to know which kind of fault produced a block.

The class is intentionally standalone.  It does not inherit from
``MultiFaultInversion`` and it does not write through ``Export.Output``,
whose writers assume a single flat fault type.
"""

import hashlib
import os
import logging

import numpy as np
from scipy import linalg, optimize

from pyginvc.Export.OutputUtils import rotate_slip
from pyginvc.Geometry.UnifiedFaultCollection import UnifiedFaultCollection
from pyginvc.Greens.BaseGreen import BaseGreen
from pyginvc.Greens.EDCMP import EDCMP
from pyginvc.Greens.Meade import Meade
from pyginvc.Greens.Nikkhoo import Nikkhoo
from pyginvc.Greens.Okada import Okada
from pyginvc.Inversion.BaseInversion import BaseInversion
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap
from pyginvc.Laplacian.TriLap import TriLap
from pyginvc.Unified.UnifiedOutput import UnifiedOutput
from pyginvc.libs import geotools as gt

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")


class _TriangleGreenMixin(object):
    """Shared cache handling for triangle Green's function classes.

    ``BaseGreen._cache_fingerprint`` only hashes rectangular geometry arrays,
    so a triangle mesh rebuilt with the same number of elements would reuse a
    stale cache.  It also saves an empty ``(0, 0)`` GPS/SAR matrix when that
    observation type is absent, which the cache loader rejects.  This mixin
    fixes both problems without touching the existing Green classes.
    """

    def _cache_fingerprint(self):
        h = hashlib.sha1()
        flt = self.flt
        arrays = [getattr(flt, 'vertex_enu', None),
                  getattr(flt, 'vertex_llh', None),
                  getattr(flt, 'element', None)]
        for name in ('llh_gps', 'llh_lev', 'llh_sar', 'unit'):
            arrays.append(getattr(self.data, name, None))
        for arr in arrays:
            if arr is None:
                continue
            a = np.ascontiguousarray(np.asarray(arr, dtype=float))
            h.update(str(a.shape).encode())
            h.update(a.tobytes())
        # The enabled slip components change which columns are populated.
        h.update(str(self.dict_green.get('greentype', [])).encode())
        h.update(str(self.dict_green.get('nu', 0.25)).encode())
        return h.hexdigest()

    def _normalize_empty_blocks(self):
        """Give absent GPS/SAR matrices a (0, 3*nf) shape for the cache."""
        nf3 = 3 * int(self.flt.nf)
        for name in ('G', 'G_sar'):
            M = np.asarray(getattr(self, name, np.empty((0, 0))))
            if M.size == 0:
                setattr(self, name, np.zeros((0, nf3)))

    def _generate_triangles(self):
        self.generate_greens()

    def build_greens(self):
        if not self.load_greens():
            self._generate_triangles()
            self._normalize_empty_blocks()
            self.save_greens()
        else:
            self._normalize_empty_blocks()
        # The existing triangle classes do not rotate their Green's functions;
        # the rectangular classes do.  Apply the same rake rotation here so
        # both fault kinds share one parameter convention.
        self.rotate_greens()
        logging.info("Building triangle Green's functions finished.")


class UnifiedNikkhoo(_TriangleGreenMixin, Nikkhoo):
    pass


class UnifiedMeade(_TriangleGreenMixin, Meade):
    pass


TRIANGLE_GREENS = {
    'nikkhoo': UnifiedNikkhoo,
    'meade': UnifiedMeade,
}


class UnifiedSlipInversion(BaseInversion):
    """Slip inversion over any mix of rectangular and triangular faults.

    Parameters
    ----------
    data : GeoDataCollection
        Aggregated GPS, leveling and InSAR observations.
    faults : UnifiedFaultCollection
        Loaded rectangular and/or triangular faults.
    origin : sequence, optional
        ``[lat, lon]`` reference origin.  Defaults to ``faults.origin``.
    config : dict, optional
        ``global`` section of the unified YAML configuration.
    export : dict, optional
        ``export`` section of the unified YAML configuration.
    """

    def __init__(self, data, faults, origin=None, config=None, export=None):
        config = config or {}
        dict_weight = {
            'smoothfactor': config.get('smoothfactor', [0.1, 0.5, 0.1]),
        }
        dict_bound = {
            'bound_switch': config.get('bound_switch', True),
        }
        super().__init__(data, dict_weight, dict_bound)

        if not isinstance(faults, UnifiedFaultCollection):
            raise TypeError('faults must be a UnifiedFaultCollection.')
        self.faults = faults
        self.origin = list(origin) if origin is not None else list(faults.origin)
        self.config = config
        self.export = export or {}
        self.outpath = self.export.get('outpath') or './'

        self.prepared = False
        self.Gmat = self.G_sar = None
        self.G_lap = None
        self.G_data = None
        self.G_ramp = None
        self.n_ramp = 0
        self.sar_ramp_flag = np.array([], dtype=bool)
        self.ramp = np.array([])
        self.ramp_sar = np.array([])
        self.bl = self.bu = None
        self.d_geod = np.array([])
        self.d_sar = np.array([])
        self.n_geod = 0
        self.n_sar = 0
        self.nf3 = 0
        self._Gw = self._Gsarw = None
        self._d2R = self._d2I = None
        self._G_ramp_w = None

    # ------------------------------------------------------------------ #
    # Preparation
    # ------------------------------------------------------------------ #
    def _ensure_prepared(self):
        if not self.prepared:
            self.prepare()

    def prepare(self):
        """Assemble observations, Green's functions, smoothing and bounds."""
        data = self.data
        faults = self.faults

        geod_parts = [x for x in (data.d_gps, data.d_lev) if len(x) > 0]
        self.d_geod = np.hstack(geod_parts) if geod_parts else np.array([])
        self.d_sar = np.asarray(data.d_sar).reshape(-1)
        self.n_geod = len(self.d_geod)
        self.n_sar = len(self.d_sar)
        self.nf3 = faults.nf3
        if self.n_geod + self.n_sar == 0:
            raise ValueError('No observations are configured.')

        self.Gmat, self.G_sar = self.build_green()
        self.G_data = np.vstack([self.Gmat, self.G_sar])
        self.G_lap = self.build_laps()

        self.W, self.WSAR = self.data.get_weights()
        self.G_ramp, self.n_ramp = self.build_ramp()

        self._Gw = (self.W @ self.Gmat if self.n_geod
                    else np.zeros((0, self.nf3)))
        self._Gsarw = (self.WSAR @ self.G_sar if self.n_sar
                       else np.zeros((0, self.nf3)))
        self._d_geodw = (self.W @ self.d_geod if self.n_geod
                         else np.array([]))
        self._d_sarw = (self.WSAR @ self.d_sar if self.n_sar
                        else np.array([]))
        self._d2R = np.hstack([self._d_geodw, self._d_sarw])
        self._d2I = np.hstack([self._d2R, np.zeros(self.nf3)])

        if self.n_ramp > 0:
            ramp_parts = []
            if self.n_geod:
                ramp_parts.append(self.W @ self.G_ramp[:self.n_geod])
            if self.n_sar:
                ramp_parts.append(self.WSAR @ self.G_ramp[self.n_geod:])
            self._G_ramp_w = (np.vstack(ramp_parts) if ramp_parts
                              else np.zeros((0, self.n_ramp)))
        else:
            self._G_ramp_w = np.zeros((0, 0))

        self.bound_switch = self.dict_bound.get('bound_switch', True)
        if self.bound_switch:
            self.bl, self.bu = self.build_bounds()

        self.prepared = True
        return self

    # ------------------------------------------------------------------ #
    # Green's functions
    # ------------------------------------------------------------------ #
    def build_green(self):
        """Build one Green block per fault and concatenate along columns."""
        os.makedirs(self.outpath, exist_ok=True)
        G_blocks, Gsar_blocks = [], []
        for i, flt, cfg, kind, name in self.faults.groups():
            g = self._green_config(i, cfg, kind)
            nf3 = 3 * int(flt.nf)
            if kind == 'rectangle':
                green = (EDCMP(flt, self.data, g) if g['method'] == 'wang'
                         else Okada(flt, self.data, g))
                green.build_greens()
            else:
                green = TRIANGLE_GREENS[g['method']](flt, self.data, g)
                green.build_greens()
            G_blocks.append(self._green_block(
                green.G, self.n_geod, nf3, name, 'geodetic'))
            Gsar_blocks.append(self._green_block(
                green.G_sar, self.n_sar, nf3, name, 'InSAR'))
            logging.info("Green's function block ready for fault %d (%s, %s).",
                         i, name, kind)
        return np.hstack(G_blocks), np.hstack(Gsar_blocks)

    def _green_config(self, index, fault_cfg, kind):
        """Merge global/per-fault Green settings and fill defaults."""
        cfg = fault_cfg or {}
        g = dict(self.config.get('green') or {})
        # rake_beta belongs to one fault only.  A global value would rotate
        # every block while the output layer rotated per-fault slips.
        g.pop('rake_beta', None)
        g.update(cfg.get('green') or {})

        g.setdefault('greentype', [1, 1, 0])
        g.setdefault('nu', 0.25)
        g.setdefault('modulus', 3e10)
        # Ramp columns are assembled once for all faults by build_ramp().
        # Leaving the per-fault switches on would add a second set of ramps.
        g['gps_ramp'] = False
        g['sar_ramp'] = False
        g['greenfile'] = self._green_file(index, g.get('greenfile'), kind)

        method = g.get('method') or g.get('grnmethod')
        method = 'okada' if method is None else str(method).lower()
        if kind == 'rectangle':
            if method in ('okada', 'rectokada'):
                g['method'] = 'okada'
            elif method in ('wang', 'edcmp', 'layerwang'):
                g['method'] = 'wang'
                # Layered-model Green's functions are generated by edgrn and
                # consumed by edcmp.
                g.setdefault('grndir', './')
                g.setdefault('grn_ssfile', 'edgrnhs.ss')
                g.setdefault('grn_dsfile', 'edgrnhs.ds')
                g.setdefault('grn_clfile', 'edgrnhs.cl')
                g.setdefault('edcmp_exe', 'edcmp')
            else:
                raise NotImplementedError(
                    "Green's function method %r is not implemented for "
                    "rectangles; use 'okada' or 'wang'." % method)
        else:
            if method not in TRIANGLE_GREENS:
                raise NotImplementedError(
                    "Green's function method %r is not implemented for "
                    "triangles; use 'nikkhoo' or 'meade'.  The local Poly3D "
                    "wrapper is incomplete, so it is not offered here."
                    % method)
            g['method'] = method
            greentype = list(g['greentype'])
            if len(greentype) != 3 or any(int(v) not in (0, 1) for v in greentype):
                raise ValueError(
                    "greentype must contain three 0/1 values, got %r"
                    % (g['greentype'],))
            g['greentype'] = greentype
            g['nu'] = float(g['nu'])
            g['modulus'] = float(g['modulus'])
        return g

    def _green_file(self, index, explicit, kind):
        """Return a unique cache path per fault, or ``'NONE'``.

        A per-fault ``greenfile`` wins over ``global.greenfile``:
        ``'NONE'`` disables caching, ``'SAVE'``/empty writes
        ``outpath/green_XX.h5``, and any other path gets an ``_XX.h5`` suffix.
        """
        requested = explicit
        if requested is None or str(requested).strip() == '':
            requested = self.config.get('greenfile')
        if requested is None or str(requested).strip() == '':
            return os.path.join(self.outpath, 'green_%02d.h5' % index)

        key = str(requested).strip().upper()
        if key == 'NONE':
            return 'NONE'
        if key == 'SAVE':
            return os.path.join(self.outpath, 'green_%02d.h5' % index)

        base = str(requested)
        root, ext = os.path.splitext(base)
        if ext.lower() != '.h5':
            if ext:
                logging.warning(
                    'Green cache %s is not HDF5; using .h5 instead.', base)
            ext = '.h5'
        return '%s_%02d%s' % (root, index, ext)

    @staticmethod
    def _green_block(M, nrows, ncols, name, label):
        """Validate one fault's Green block; allow (0, ncols) when no data."""
        M = np.asarray(M)
        if M.size == 0:
            if nrows == 0:
                return np.zeros((0, ncols))
            raise ValueError(
                "Green's function block for fault %s (%s) is empty but %d "
                "observations are expected; delete the cached Green file or "
                "change greenfile." % (name, label, nrows))
        if M.shape != (nrows, ncols):
            raise ValueError(
                "Green's function block for fault %s (%s) has shape %s, "
                "expected %s; delete the cached Green file or change "
                "greenfile." % (name, label, M.shape, (nrows, ncols)))
        return M

    # ------------------------------------------------------------------ #
    # Smoothing
    # ------------------------------------------------------------------ #
    def _lap_method(self, cfg, kind):
        """Resolve a per-fault Laplacian method with a global fallback."""
        value = (cfg or {}).get('lapmethod')
        global_lap = self.config.get('lapmethod')
        if value is None and isinstance(global_lap, dict):
            value = global_lap.get(kind)
        if value is None and isinstance(global_lap, (str, int, float)):
            value = global_lap
        if value is None:
            value = '1' if kind == 'rectangle' else '3'
        return str(value).strip()

    def build_laps(self):
        """Build a block-diagonal smoothing matrix in fault order."""
        laps = []
        for i, flt, cfg, kind, name in self.faults.groups():
            cfg = cfg or {}
            method = self._lap_method(cfg, kind)
            weight = cfg.get('smoothing_weight', 1.0)
            weight = float(weight)
            if not np.isfinite(weight) or weight < 0:
                raise ValueError(
                    'smoothing_weight for fault %s must be finite and '
                    'non-negative, got %r.' % (name, weight))

            if kind == 'rectangle':
                if method not in ('1', '4'):
                    raise ValueError(
                        "Unsupported lapmethod %r for rectangle fault %s; "
                        "use '1' or '4'." % (method, name))
                # RectPlaneLap cannot build the finite-difference Laplacian
                # for a single patch along either direction.
                if method == '1' and (int(flt.nsegs) < 2 or int(flt.ndeps) < 2):
                    logging.info(
                        'Fault %s has fewer than two patches along strike or '
                        'dip; using Tikhonov smoothing.', name)
                    method = '4'
                bcs = cfg.get('bcs', [0, 0, 0, 0])
                if len(bcs) != 4:
                    raise ValueError(
                        'bcs for fault %s must have four values, got %r.'
                        % (name, bcs))
                lap = RectPlaneLap(flt.nsegs, flt.ndeps, bcs, method)
            else:
                if method not in ('1', '2', '3', '4'):
                    raise ValueError(
                        "Unsupported lapmethod %r for triangle fault %s; "
                        "use '1', '2', '3' or '4'." % (method, name))
                lap = TriLap(flt.vertex_enu, flt.element, method)

            lap.build_laps()
            L = np.asarray(lap.G_lap, dtype=float)
            if kind == 'triangle' and not np.isfinite(L).all():
                # The geometric TriLap modes divide by a neighbour-distance
                # sum that is zero for one or two isolated elements.  A
                # Tikhonov matrix is well defined for any mesh size.
                logging.warning(
                    'Laplacian method %s produced non-finite values for fault '
                    '%s; falling back to method 4.', method, name)
                lap = TriLap(flt.vertex_enu, flt.element, '4')
                lap.build_laps()
                L = np.asarray(lap.G_lap, dtype=float)
                method = '4'
            expected = (3 * int(flt.nf), 3 * int(flt.nf))
            if L.shape != expected:
                raise ValueError(
                    'Laplacian for fault %s has shape %s, expected %s.'
                    % (name, L.shape, expected))
            laps.append(L if weight == 1.0 else weight * L)
            logging.info('Laplacian ready for fault %d (%s, %s, method %s).',
                         i, name, kind, method)
        return linalg.block_diag(*laps) if laps else np.zeros((self.nf3, self.nf3))

    # ------------------------------------------------------------------ #
    # Bounds
    # ------------------------------------------------------------------ #
    def build_bounds(self):
        """Build lower/upper bounds for every fault plus the ramp columns."""
        bl_parts, bu_parts = [], []
        for i, flt, cfg, kind, name in self.faults.groups():
            cfg = cfg or {}
            b = dict(cfg.get('bound') or {})
            if kind == 'rectangle':
                b.setdefault('ss_range', [-np.inf, np.inf])
                b.setdefault('ds_range', [-np.inf, np.inf])
                b.setdefault('op_range', [-np.inf, np.inf])
                bu, bl = BaseInversion.set_slip_bounds(
                    flt.nsegs, flt.ndeps, 0, **b)
            else:
                bu, bl = self._triangle_bounds(flt, b, name)
            bl = np.asarray(bl, dtype=float).reshape(-1)
            bu = np.asarray(bu, dtype=float).reshape(-1)
            expected = 3 * int(flt.nf)
            if bl.size != expected or bu.size != expected:
                raise ValueError(
                    'Bounds for fault %s have lengths %d/%d, expected %d.'
                    % (name, bl.size, bu.size, expected))
            if np.any(bl > bu):
                raise ValueError(
                    'Lower bound exceeds upper bound for fault %s.' % name)
            bl_parts.append(bl)
            bu_parts.append(bu)

        bl = np.concatenate(bl_parts) if bl_parts else np.array([])
        bu = np.concatenate(bu_parts) if bu_parts else np.array([])
        if self.n_ramp > 0:
            bl = np.hstack([bl, np.full(self.n_ramp, -np.inf)])
            bu = np.hstack([bu, np.full(self.n_ramp, np.inf)])
        return bl, bu

    @staticmethod
    def _triangle_bounds(flt, bound_cfg, name):
        """Build bounds for a triangle mesh from ranges or slip files.

        Triangle slip files carry ``(strike, dip, opening)`` in columns 4-6
        (zero-based 3-5), matching the element format used by ``Triangle``.
        """
        nf = int(flt.nf)
        nparam = 3 * nf
        bl = np.full(nparam, -np.inf)
        bu = np.full(nparam, np.inf)

        for component, key in enumerate(('ss_range', 'ds_range', 'op_range')):
            value = bound_cfg.get(key)
            if value is None:
                continue
            value = np.asarray(value, dtype=float).reshape(-1)
            if value.size != 2:
                raise ValueError(
                    '%s for triangle fault %s must have two values, got %r.'
                    % (key, name, value.tolist()))
            bl[component::3] = value[0]
            bu[component::3] = value[1]

        slip_lb = bound_cfg.get('slip_lb', '')
        slip_ub = bound_cfg.get('slip_ub', '')
        has_lb = bool(str(slip_lb or '').strip())
        has_ub = bool(str(slip_ub or '').strip())
        if has_lb != has_ub:
            raise ValueError(
                'slip_lb and slip_ub must be given together for fault %s; '
                'got %r and %r.' % (name, slip_lb, slip_ub))
        if has_lb and has_ub:
            lb = np.genfromtxt(slip_lb, comments='#', usecols=[3, 4, 5])
            ub = np.genfromtxt(slip_ub, comments='#', usecols=[3, 4, 5])
            lb = np.asarray(lb, dtype=float).reshape(-1)
            ub = np.asarray(ub, dtype=float).reshape(-1)
            if lb.size != nparam or ub.size != nparam:
                raise ValueError(
                    'slip_lb/slip_ub for triangle fault %s must contain %d '
                    'values (%d patches x 3 components); got %d and %d.'
                    % (name, nparam, nf, lb.size, ub.size))
            bl[:] = lb
            bu[:] = ub
            logging.info(
                'Triangle slip bounds for %s read from %s and %s.',
                name, slip_lb, slip_ub)

        supported = {'ss_range', 'ds_range', 'op_range', 'slip_lb', 'slip_ub'}
        ignored = sorted(set(bound_cfg) - supported)
        if ignored:
            logging.warning(
                'Ignoring rectangle-only bound keys for triangle fault %s: %s',
                name, ', '.join(ignored))

        return bu, bl

    # ------------------------------------------------------------------ #
    # Ramps
    # ------------------------------------------------------------------ #
    def build_ramp(self):
        """Build per-file GPS and per-track InSAR ramp columns.

        Row order matches ``[GPS (ndim), leveling (1), InSAR (1)]``.  Every
        non-empty InSAR track gets at least a constant reference level; the
        east/north tilt columns are added only when its config sets
        ``ramp: true``.
        """
        data = self.data
        origin = self.origin
        ndim = data.ndim

        gps_entries = list(data.gps_pairs())
        sar_entries = list(data.sar_pairs())

        gps_need = any(bool((g or {}).get('ramp', False)) and
                       len(gps.llh_gps) > 0 for gps, g in gps_entries)
        sar_need = any(len(sar.llh_sar) > 0 for sar, _ in sar_entries)

        xy_gps = gt.llh2utm_batch(data.llh_gps, origin) if gps_need else None
        xy_sar = gt.llh2utm_batch(data.llh_sar, origin) if sar_need else None

        blocks = []
        n_ramp = 0

        sta0 = 0
        for gps, g in gps_entries:
            g = g or {}
            n = len(gps.llh_gps)
            if g.get('ramp', False) and n > 0:
                helper = BaseGreen(None, gps, {'greenfile': 'NONE'})
                M = helper.make_green_gps_ramp(
                    xy_gps[sta0:sta0 + n], ndim, g.get('ramp_method', 1))
                blocks.append(M)
                n_ramp += M.shape[1]
            else:
                blocks.append(np.zeros((n * ndim, 0)))
            sta0 += n

        # Leveling has no ramp; it only contributes rows.
        blocks.append(np.zeros((len(data.d_lev), 0)))

        sta0 = 0
        sar_ramp_flag = []
        for sar, s in sar_entries:
            s = s or {}
            n = len(sar.llh_sar)
            if n > 0:
                M = BaseGreen.make_green_sar_ramp(xy_sar[sta0:sta0 + n])
                if not s.get('ramp', False):
                    M = M[:, :1]
                blocks.append(M)
                n_ramp += M.shape[1]
                sar_ramp_flag.append(bool(s.get('ramp', False)))
            else:
                blocks.append(np.zeros((0, 0)))
                sar_ramp_flag.append(False)
            sta0 += n
        self.sar_ramp_flag = np.asarray(sar_ramp_flag, dtype=bool)

        if n_ramp == 0:
            return np.empty((0, 0)), 0
        return linalg.block_diag(*blocks), n_ramp

    # ------------------------------------------------------------------ #
    # Solve
    # ------------------------------------------------------------------ #
    def _solver(self):
        value = self.config.get('solver', 'bvls')
        value = 'bvls' if value is None else str(value).strip().lower()
        if value in ('lsq', 'least_squares', 'linear'):
            value = 'lstsq'
        if value not in ('bvls', 'trf', 'lstsq'):
            raise ValueError(
                "Unsupported solver %r; use 'bvls', 'trf' or 'lstsq'."
                % (value,))
        return value

    def _solve_linear(self, A, b):
        """Solve ``A x = b`` with the configured method and bounds."""
        solver = self._solver()
        use_bounds = (self.bound_switch and self.bl is not None and
                      not (np.isneginf(self.bl).all() and
                           np.isposinf(self.bu).all()))
        if solver == 'lstsq' or not use_bounds:
            if solver != 'lstsq' and self.bound_switch and self.bl is not None:
                logging.info(
                    'All parameter bounds are infinite; using unconstrained '
                    'least squares.')
            x, *_ = np.linalg.lstsq(A, b, rcond=None)
            logging.info('Unconstrained LSQ finished.')
            return x

        lower = np.array(self.bl, dtype=float, copy=True)
        upper = np.array(self.bu, dtype=float, copy=True)
        # lsq_linear requires strictly ordered bounds.  An equal pair means a
        # parameter pinned to one value (for example opening fixed at zero);
        # nudge the upper bound by one representable step.
        fixed = lower == upper
        if fixed.any():
            upper[fixed] = np.nextafter(lower[fixed], np.inf)

        res = optimize.lsq_linear(
            A, b, (lower, upper), method=solver, max_iter=20000,
            lsmr_tol='auto')
        logging.info('Constrained LSQ (%s) finished: %s.',
                     solver, res.message)
        return res.x

    def solve(self, smo):
        """Invert for one smoothing factor and return the solution record."""
        self._ensure_prepared()

        G2I_slip = np.vstack([self._Gw,
                              self._Gsarw,
                              smo * self.G_lap])
        d2I = self._d2I

        if self.n_ramp > 0:
            G2I = np.hstack([G2I_slip,
                             np.vstack([self._G_ramp_w,
                                        np.zeros((self.nf3, self.n_ramp))])])
        else:
            G2I = G2I_slip

        x = self._solve_linear(G2I, d2I)
        if self.bound_switch and self.bl is not None:
            x = np.clip(x, self.bl, self.bu)
            fixed = np.asarray(self.bl) == np.asarray(self.bu)
            if fixed.any():
                # Pin closed intervals exactly; the solver only sees a
                # one-ulp-wide interval to satisfy its bounding rules.
                x[fixed] = np.asarray(self.bl)[fixed]
        slip = x[:self.nf3]

        ramp_hat = x[self.nf3:] if self.n_ramp > 0 else np.array([])
        if self.n_ramp > 0 and self.n_sar > 0:
            ramp_sar = self.G_ramp[self.n_geod:] @ ramp_hat
        else:
            ramp_sar = np.zeros(self.n_sar)

        pred = self.G_data @ slip
        if self.n_ramp > 0:
            pred = pred + self.G_ramp @ ramp_hat
        r_vec, _, misfit_gps, misfit_sar, _ = self.compute_residuals(
            self._d2R, None, x, self.W, self.WSAR,
            self.n_geod, self.n_sar, dhat=pred)

        lap_slip = self.G_lap @ slip
        roughness = float(lap_slip @ lap_slip)
        Mo, Mw = self.faults.total_moment(
            slip, default_modulus=self._default_modulus())

        return {
            'smo': smo,
            'x': x,
            'slip': slip,
            'dhat': pred,
            'r': r_vec,
            'misfit_gps': float(misfit_gps),
            'misfit_sar': float(misfit_sar),
            'roughness': roughness,
            'Mo': Mo,
            'Mw': Mw,
            'ramp': ramp_hat,
            'ramp_sar': ramp_sar,
        }

    def slip_by_fault(self, slip):
        """Return per-fault ``(nf, 3)`` slip arrays with rake restored."""
        slip = np.asarray(slip, dtype=float).reshape(-1)
        if slip.size != self.nf3:
            raise ValueError('slip has %d values, expected %d.'
                             % (slip.size, self.nf3))
        out = []
        for (i, flt, cfg, kind, name), beta in zip(
                self.faults.groups(), self.faults.rake_betas()):
            offset = self.faults.offsets[i]
            block = slip[offset:offset + 3 * int(flt.nf)].reshape(-1, 3)
            out.append((name, kind, rotate_slip(block, beta)))
        return out

    def _default_modulus(self):
        mod = self.config.get('modulus')
        if mod is None:
            mod = (self.config.get('green') or {}).get('modulus', 3e10)
        return float(mod)

    # ------------------------------------------------------------------ #
    # Driver
    # ------------------------------------------------------------------ #
    def run(self, smo_facts=None, output=None):
        """Run every smoothing factor, write output and keep the last model."""
        self._ensure_prepared()

        if smo_facts is None:
            smo_facts = self.get_smoothing_factors()
        smo_facts = np.asarray(smo_facts).reshape(-1)
        if smo_facts.size == 0:
            raise ValueError('No smoothing factors to invert for.')

        outpath = self.export.get('outpath') or './'
        os.makedirs(outpath, exist_ok=True)

        if output is None:
            output = UnifiedOutput(self.data,
                                   self.faults,
                                   sol=self,
                                   outdir=outpath,
                                   config=self.config,
                                   export=self.export)

        results = []
        slip, misfit_gps, misfit_sar = [], [], []
        rough, moment = [], []
        last = None
        for k, smo in enumerate(smo_facts):
            r = self.solve(smo)
            output.WriteSmoothing(k, r)
            last = r

            results.append((r['smo'], r['misfit_gps'], r['misfit_sar'],
                            r['roughness'], r['Mo'], r['Mw']))
            slip.append(r['slip'])
            misfit_gps.append(r['misfit_gps'])
            misfit_sar.append(r['misfit_sar'])
            rough.append(r['roughness'])
            moment.append((r['Mo'], r['Mw']))
            logging.info('sm=%.3f  WRSS_geod=%10.3f  WRSS_sar=%10.3f  Mw=%.2f',
                         smo, r['misfit_gps'], r['misfit_sar'], r['Mw'])

        self.smo_facts = np.asarray(smo_facts)
        self.slip = np.asarray(slip)
        self.misfit_gps = np.asarray(misfit_gps)
        self.misfit_sar = np.asarray(misfit_sar)
        self.misfit = self.misfit_gps + self.misfit_sar
        self.ruff = np.asarray(rough)
        self.smoothness = self.ruff * self.smo_facts ** 2
        self.moment = np.asarray(moment)
        self.dhat = np.asarray(last['dhat'])
        self.r = np.asarray(last['r'])
        self.ramp = np.asarray(last['ramp'])
        self.ramp_sar = np.asarray(last['ramp_sar'])

        output.WriteHDF5()
        output.WriteSummary()
        output.WriteLCurve()
        output.WriteGMTSlip()
        return results
