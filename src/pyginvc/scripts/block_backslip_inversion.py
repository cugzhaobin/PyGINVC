#!/usr/bin/env python
"""
Block-model backslip interseismic inversion, driven by a YAML config file.

    python -m pyginvc.scripts.block_backslip_inversion --cfgfile block_backslip.yaml

The config picks the model.  With one meshed fault and dict_coupling.linked
true (the default), LinkedBlockBackslipInversion fits block Euler vectors and
one bounded, smoothed locking coefficient per patch, with the slip deficit tied
to the long-term relative block motion:

    v = B omega - G (chi * s_lt)

Set dict_coupling.linked: false for the older free-slip model
(BlockBackslipInversion).  A `dict_segment` section adds a second, unsmoothed
fault and selects GenBlockBackslipInversion; the linked model is not available
for that configuration.

Sections the driver reads:
    dict_data      gpsfile / sarfile / levfile / gfiletype
    dict_block     blockfile, assign
    dict_fault     faultfile, nsegs, ndeps, doSubFault, origin
    dict_segment   same shape as dict_fault; optional
    dict_lap       method
    dict_coupling  compute, linked, chi_range, max_iter, tol, min_slt,
                   require_exact, offsets_km, footwall_block, block_pair,
                   nseg_mean
    dict_export    outpath, view, vecscale, write_forward,
                   write_param_contrib

Sections forwarded to the library verbatim (the driver never inspects them, so
library options can be used without touching this file):
    dict_green, dict_weight, dict_bound

Outputs (last smoothing factor)
------------------------------
Model files
    euler_est.dat             block, wx, wy, wz (buildG convention)
    euler_pole.dat            equivalent Euler pole (lat, lon, deg/Myr, mas/yr)
    slipmodel.faultgeom       geometry + fitted slip deficit + s_lt + chi
    slipmodel_segment.faultgeom
    coupling.faultgeom        same coupling columns for the mesh
    chi_profile.dat

Station-level forward decomposition (mm/yr)
    gps_forward.dat           lon lat obs sig rot backslip total resid wres station
    gps_mod_rotation.gmtvec   rigid block rotation contribution
    gps_mod_slip.gmtvec       fault backslip (negative dislocation) contribution
    gps_mod_total.gmtvec      rotation + backslip
    gps_residual.gmtvec       raw and weighted residuals
    gps_fit_stats.dat         RMS/correlation/WRSS summary

Per-parameter surface contributions (mm/yr)
    gps_contrib_euler.dat     one E/N pair per Euler component (wx, wy, wz)
    gps_contrib_block.gmtvec  one E/N pair per block (sum over its wx, wy, wz)
    gps_contrib_fault.dat     one E/N pair per fault patch (sum over ss/ds/op,
                              or the fitted chi in the linked model)
    param_contrib_summary.dat value + RMS/max contribution of every fitted
                              parameter, and of derived physical slip components
    param_contrib_ranked.dat  the same table sorted by surface contribution
    contributions.h5          full G and contribution matrices at every station
                              (optional; written when write_param_contrib is true)
"""
import argparse
import logging
import os
import sys

import numpy as np
import yaml
from scipy.linalg import block_diag

from pyginvc.Blocks.Blocks import BlockModel
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Patch import Fault
from pyginvc.Greens.Okada import Okada
from pyginvc.Inversion.Coupling import Coupling
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap
from pyginvc.Inversion.BlockBackslipInversion import (
    BlockBackslipInversion, LinkedBlockBackslipInversion)
from pyginvc.Inversion.GenBlockBackslipInversion import GenBlockBackslipInversion

REQUIRED_SECTIONS = ['dict_data', 'dict_block', 'dict_fault',
                     'dict_green', 'dict_weight', 'dict_bound']

SLIPMODEL_HEADER = """\
# Block backslip inversion -- slip deficit model
#
#   1 width_km     2 depth_km    3 dip_deg      4 lat0        5 lon0
#   6 lat1         7 lon1
#   8 ss_mm_yr     9 ds_mm_yr   10 op_mm_yr
#  11 s_lt_ss     12 s_lt_ds    13 |s_lt|      14 chi
#
# ss/ds/op : fitted slip DEFICIT rates (backslip), not slip rates.  In the
#            linked model (dict_coupling.linked: true) this is chi * s_lt.
#
# s_lt_ss / s_lt_ds : long-term strike-slip and dip-slip rates of the hanging
#            wall relative to the footwall, projected onto the fault's slip
#            vector and sign-calibrated against the Okada Green's function.
#            Evaluated at the SURFACE trace of each patch's along-strike
#            column, from the relative rigid rotation of the two flanking
#            blocks.  Units mm/yr.
#
# |s_lt|   : sqrt(s_lt_ss^2 + s_lt_ds^2).
#
# chi      : locking coefficient.  In the linked model it is FITTED, bounded
#            by dict_coupling.chi_range (default [0, 1]), and smoothed.  In
#            the free-slip fallback it is DERIVED as the projection
#            chi = s . s_lt / |s_lt|^2 and is not clipped.
#
# nan in columns 11-14 means |s_lt| is below dict_coupling.min_slt, no block
# pair was found for that patch, or the coupling stage was skipped.  The
# linked model instead fixes those patches at s_lt = chi = 0 (they carry no
# slip deficit and do not smooth against their neighbours).
"""



class ConfigError(Exception):
    """Raised for anything wrong in the config file, so main() can report it
    as a one-line message instead of a traceback."""


_MISSING = object()


def cfg_get(cfg, *keys, **kw):
    """
    cfg['dict_fault']['nsegs'], but the error names the yaml path that is
    missing rather than raising a bare KeyError.  Pass default=... to make a
    key optional.
    """
    default = kw.pop('default', _MISSING)
    if kw:
        raise TypeError('unexpected keyword arguments: %s' % list(kw))
    node = cfg
    for i, k in enumerate(keys):
        if not isinstance(node, dict) or k not in node or node[k] is None:
            if default is not _MISSING:
                return default
            raise ConfigError('missing required key: %s' % '.'.join(keys[:i + 1]))
        node = node[k]
    return node


def cfg_section(cfg, name, required=True):
    """A whole section, forwarded to the library untouched."""
    if name not in cfg or cfg[name] is None:
        if required:
            raise ConfigError('missing required section: %s' % name)
        return None
    if not isinstance(cfg[name], dict):
        raise ConfigError('%s must be a mapping, got %s' % (name, type(cfg[name]).__name__))
    return cfg[name]


def load_config(path):
    if not os.path.isfile(path):
        raise ConfigError('config file not found: %s' % path)
    with open(path, 'r') as fid:
        cfg = yaml.load(fid.read(), Loader=yaml.FullLoader)
    if not isinstance(cfg, dict):
        raise ConfigError('%s does not contain a YAML mapping' % path)
    return cfg


class BlockBackslipRun:
    """One config-driven run: build the model, invert, export, report locking."""

    def __init__(self, args):
        self.cfg  = load_config(args.cfgfile)
        self.args = args
        self._apply_overrides()

        for name in REQUIRED_SECTIONS:
            cfg_section(self.cfg, name)

        self.outpath  = cfg_get(self.cfg, 'dict_export', 'outpath', default='.')
        self.view     = bool(cfg_get(self.cfg, 'dict_export', 'view', default=False))
        self.vecscale = float(cfg_get(self.cfg, 'dict_export', 'vecscale', default=200))
        self.coupling_on = bool(cfg_get(self.cfg, 'dict_coupling', 'compute', default=False))
        # The linked model (s = chi * s_lt) is the default whenever coupling is
        # on.  Set dict_coupling.linked: false to get the old free-s model.
        self.linked = bool(cfg_get(self.cfg, 'dict_coupling', 'linked', default=True))

        # validate everything cheap up front -- a typo should not have to wait
        # for the geodetic data and the Green's functions to load
        self.assign = cfg_get(self.cfg, 'dict_block', 'assign', default='contains')
        if self.assign not in ('contains', 'nearest'):
            raise ConfigError("dict_block.assign must be 'contains' or 'nearest', got %r"
                              % self.assign)

    # ------------------------------------------------------------- overrides
    def _set(self, section, key, value):
        node = self.cfg.get(section)
        if not isinstance(node, dict):
            node = {}
            self.cfg[section] = node
        node[key] = value

    def _apply_overrides(self):
        """Command-line options win over the config file."""
        a = self.args
        if a.gpsfile:
            self._set('dict_data', 'gpsfile', a.gpsfile)
        if a.outpath:
            self._set('dict_export', 'outpath', a.outpath)
        if a.smoothfactor:
            self._set('dict_weight', 'smoothfactor', list(a.smoothfactor))
        if a.no_coupling:
            self._set('dict_coupling', 'compute', False)

    # -------------------------------------------------------------- pipeline
    def run(self):
        # _coupling writes coupling.faultgeom before _export runs, so the
        # directory has to exist before either of them.
        if self.outpath and not os.path.isdir(self.outpath):
            os.makedirs(self.outpath)
            logging.info('created output directory %s',
                         os.path.abspath(self.outpath))

        data   = self._data()
        blocks = self._blocks(data)
        flts, greens, lap = self._model(data)
        sol    = self._invert(data, blocks, flts, greens, lap)
        # coupling before export: slipmodel*.faultgeom carries s_lt and chi
        cpls   = (self._coupling(sol, blocks, flts, greens)
                  if self.coupling_on else (None, None))
        self._export(sol, data, blocks, flts, greens, cpls)
        if self.view:
            self._plot(sol, data, blocks)
        return sol

    def _data(self):
        c = self.cfg
        data = GeoData(cfg_get(c, 'dict_data', 'gpsfile'),
                       cfg_get(c, 'dict_data', 'sarfile', default=''),
                       cfg_get(c, 'dict_data', 'levfile', default=''),
                       cfg_get(c, 'dict_data', 'gfiletype', default='GMT2D'))
        data.load_data()
        return data

    def _blocks(self, data):
        c        = self.cfg
        fallback = self.assign

        blocks = BlockModel(cfg_get(c, 'dict_block', 'blockfile'))
        n_gps  = len(data.llh_gps)
        gps_block_id = blocks.assign_point(data.llh_gps[:, [1, 0]], 'gps', fallback=fallback)

        n_out = int((gps_block_id < 0).sum())
        if n_out:
            logging.warning('%d of %d GPS sites lie inside no block polygon and are '
                            'DROPPED from the inversion (dict_block.assign=%s).  Set '
                            "assign: 'nearest' to keep them.", n_out, n_gps, fallback)
        data.process_data(gps_block_id, [])
        logging.info('%d GPS sites in the inversion (%d dropped), %d blocks',
                     len(data.llh_gps), n_out, blocks.nblock)
        return blocks

    def _make_fault(self, section):
        c   = self.cfg
        flt = Fault(cfg_get(c, section, 'faultfile'),
                    int(cfg_get(c, section, 'nsegs')),
                    int(cfg_get(c, section, 'ndeps')),
                    bool(cfg_get(c, section, 'doSubFault', default=False)),
                    origin=cfg_get(c, section, 'origin', default=[]))
        flt.load_fault()
        return flt

    def _model(self, data):
        c          = self.cfg
        dict_green = cfg_section(c, 'dict_green')

        flt_mesh  = self._make_fault('dict_fault')
        seg_cfg   = c.get('dict_segment')
        flt_seg   = self._make_fault('dict_segment') if isinstance(seg_cfg, dict) and seg_cfg else None

        green_mesh = Okada(flt_mesh, data, dict_green)
        green_mesh.build_greens()
        green_seg  = None
        if flt_seg is not None:
            green_seg = Okada(flt_seg, data, dict_green)
            green_seg.build_greens()

        if 'bcs' not in dict_green:
            raise ConfigError('missing required key: dict_green.bcs')
        lap = RectPlaneLap(flt_mesh.nsegs, flt_mesh.ndeps, dict_green['bcs'],
                           method=cfg_get(c, 'dict_lap', 'method', default='1'))
        lap.build_laps()

        return (flt_mesh, flt_seg), (green_mesh, green_seg), lap

    def _invert(self, data, blocks, flts, greens, lap):
        c      = self.cfg
        flt_mesh, flt_seg     = flts
        green_mesh, green_seg = greens
        dict_weight   = cfg_section(c, 'dict_weight')
        dict_bound    = cfg_section(c, 'dict_bound')
        dict_coupling = cfg_section(c, 'dict_coupling', required=False) or {}

        G_gps_block = block_diag(*blocks.buildG_block(data.ndim, block_type='gps'))
        nblock      = blocks.nblock

        if self.linked and not self.coupling_on:
            logging.warning('dict_coupling.linked is true but compute is false; '
                            'falling back to the free-slip model because chi is '
                            'not being estimated')
            linked = False
        elif self.linked and flt_seg is not None:
            logging.warning('the linked (s = chi * s_lt) model currently supports '
                            'one meshed fault only; falling back to the free-slip '
                            'model for the mesh + segment configuration')
            linked = False
        else:
            linked = self.linked

        if flt_seg is None:
            if linked:
                logging.info('model: mesh only (%d patches) + linked chi, %d blocks',
                             flt_mesh.nf, nblock)
                sol = LinkedBlockBackslipInversion(
                    data, blocks, flt_mesh, green_mesh, lap,
                    dict_weight, dict_bound, dict_coupling)
            else:
                logging.info('model: mesh only (%d patches) + free slip, %d blocks',
                             flt_mesh.nf, nblock)
                sol = BlockBackslipInversion(data, flt_mesh, green_mesh, lap,
                                             dict_weight, dict_bound)
        else:
            logging.info('model: mesh (%d patches) + segment (%d patches), %d blocks',
                         flt_mesh.nf, flt_seg.nf, nblock)
            sol = GenBlockBackslipInversion(data, flt_mesh, flt_seg, green_mesh,
                                            green_seg, lap, dict_weight, dict_bound)
        sol.run_inversion(nblock, G_gps_block)
        return sol

    def _path(self, name):
        return os.path.join(self.outpath, name)

    def _export(self, sol, data, blocks, flts, greens, cpls):
        """
        Write the last smoothing factor's model AND its forward decomposition.

        The model is

            v(x) = B(x) omega          (rigid block rotation)
                 - G(x) s              (fault backslip / negative dislocation)

        with s = chi * s_lt in the linked model.  Everything below is computed
        from the same unweighted kernels the inversion used, so the pieces sum
        to the prediction and observed - prediction is the residual:

            rotation + backslip = dhat
            (observed - dhat) / sigma = sol.r   (the weighted residual)

        The per-parameter files answer "how much does parameter j move the
        surface": column j of G times x_j.  In the linked model the fitted
        fault parameters are chi; the physical slip components chi*s_lt are
        also exported and marked as derived.
        """
        last     = -1
        flt_mesh, flt_seg = flts
        green_mesh, green_seg = greens
        nf_mesh  = flt_mesh.nf
        nf_seg   = flt_seg.nf if flt_seg is not None else 0
        nf_total = nf_mesh + nf_seg
        nblock   = blocks.nblock

        write_forward = bool(cfg_get(self.cfg, 'dict_export', 'write_forward',
                                     default=True))
        write_contrib = bool(cfg_get(self.cfg, 'dict_export',
                                     'write_param_contrib', default=True))

        # ---------------------------------------------------- model values
        # sol.slip, not sol.x: in the linked model the fault part of x is chi,
        # the slip deficit is chi * s_lt.
        omega = np.asarray(sol.euler[last], dtype=float).ravel()
        slip  = np.asarray(sol.slip[last], dtype=float).ravel()
        euler = omega.reshape(nblock, 3)
        slip3 = slip.reshape(nf_total, 3)
        linked = getattr(sol, 'chi', None) is not None
        chi = (np.asarray(sol.chi[last], dtype=float).ravel()
               if linked else None)

        self._euler_files(omega)
        self._slipmodel('slipmodel.faultgeom', flt_mesh,
                        slip[:3*nf_mesh], cpls[0])
        if flt_seg is not None:
            self._slipmodel('slipmodel_segment.faultgeom', flt_seg,
                            slip[3*nf_mesh:], cpls[1])

        if not (write_forward or write_contrib):
            logging.info('wrote euler/slipmodel files to %s (forward and '
                         'contribution outputs disabled)',
                         os.path.abspath(self.outpath))
            return

        # ------------------------------------------- unweighted GPS kernels
        n_gps = len(data.llh_gps)
        n_obs = len(data.d_gps)
        if n_obs != 2 * n_gps:
            raise ValueError('expected two GPS components per station, got '
                             '%d components for %d stations' % (n_obs, n_gps))

        G_rot = np.asarray(block_diag(
            *blocks.buildG_block(data.ndim, block_type='gps')), dtype=float)
        G_slip_mesh = -np.asarray(green_mesh.G[:n_obs], dtype=float)
        if flt_seg is not None:
            G_slip_seg = -np.asarray(green_seg.G[:n_obs], dtype=float)
            G_slip = np.hstack([G_slip_mesh, G_slip_seg])
        else:
            G_slip = G_slip_mesh
        if G_rot.shape != (n_obs, 3*nblock):
            raise ValueError('rotation kernel has shape %r, expected (%d, %d)'
                             % (G_rot.shape, n_obs, 3*nblock))
        if G_slip.shape != (n_obs, 3*nf_total):
            raise ValueError('backslip kernel has shape %r, expected (%d, %d)'
                             % (G_slip.shape, n_obs, 3*nf_total))

        # Required total per station: rotation + backslip.
        contrib_rot  = G_rot * omega[None, :]
        pred_rot     = contrib_rot.sum(axis=1)
        contrib_slip = G_slip * slip[None, :]
        pred_slip    = contrib_slip.sum(axis=1)
        total        = pred_rot + pred_slip

        obs   = np.asarray(data.d_gps, dtype=float).ravel()
        sigma = np.asarray(data.get_sigma(), dtype=float).ravel()
        sigma = np.where(sigma > 0, sigma, 1.0)
        resid = obs - total
        wres  = resid / sigma

        # Consistency checks against the inversion's own arrays.  A mismatch
        # here means the station order or the kernel convention changed, and
        # the decomposition would be misleading.
        dhat = np.asarray(sol.dhat[last], dtype=float)[:n_obs]
        rinv = np.asarray(sol.r[last], dtype=float)[:n_obs]
        err_pred  = float(np.max(np.abs(total - dhat)))
        err_resid = float(np.max(np.abs(wres - rinv)))
        if err_pred > 1.0e-6:
            logging.warning('forward decomposition mismatch: max|rotation + '
                            'backslip - dhat| = %.3e mm/yr', err_pred)
        if err_resid > 1.0e-6:
            logging.warning('weighted residual mismatch: max|(obs-total)/sigma '
                            '- sol.r| = %.3e', err_resid)
        logging.info('forward check: max|rotation + backslip - dhat| = %.3e mm/yr; '
                     'max|weighted residual difference| = %.3e',
                     err_pred, err_resid)

        # chi parameter contribution (linked model): sum_k -G(:,3p+k)*s_lt(3p+k)*chi
        G_chi = chi_contrib = None
        if linked:
            if chi.size != nf_mesh:
                raise ValueError('chi has %d entries, expected %d mesh patches'
                                 % (chi.size, nf_mesh))
            if cpls[0] is None:
                raise ValueError('linked model has no Coupling object to supply s_lt')
            s_lt = np.asarray(cpls[0].s_lt_vec, dtype=float)
            if s_lt.shape != (nf_mesh, 3):
                raise ValueError('s_lt_vec has shape %r, expected (%d, 3)'
                                 % (s_lt.shape, nf_mesh))
            # Physical slip-component coefficients, then the chi coefficient.
            G_phys = G_slip_mesh * s_lt.reshape(-1)[None, :]
            G_chi  = G_phys.reshape(n_obs, nf_mesh, 3).sum(axis=2)
            chi_contrib = G_chi * chi[None, :]
            err_chi = float(np.max(np.abs(
                chi_contrib.sum(axis=1) - contrib_slip[:, :3*nf_mesh].sum(axis=1))))
            if err_chi > 1.0e-6:
                logging.warning('chi contribution does not reproduce the mesh '
                                'backslip: max difference %.3e mm/yr', err_chi)

        lon = np.asarray(data.llh_gps[:, 1], dtype=float)
        lat = np.asarray(data.llh_gps[:, 0], dtype=float)
        stations = np.asarray(data.station_gps).astype(str)
        if stations.size != n_gps:
            stations = np.array(['sta%05d' % i for i in range(n_gps)])

        euler_names = ['blk%02d_%s' % (b, c)
                       for b in range(nblock) for c in ('wx', 'wy', 'wz')]
        patch_names = (['patch%03d' % p for p in range(nf_mesh)]
                       + ['seg%03d' % p for p in range(nf_seg)])
        slip_names = []
        for p in range(nf_mesh):
            slip_names += ['patch%03d_%s' % (p, c) for c in ('ss', 'ds', 'op')]
        for p in range(nf_seg):
            slip_names += ['seg%03d_%s' % (p, c) for c in ('ss', 'ds', 'op')]

        # Patch-aggregated backslip contribution (sum over ss/ds/op).
        contrib_patch = contrib_slip.reshape(n_obs, nf_total, 3).sum(axis=2)

        # ------------------------------------------------- station tables
        if write_forward:
            self._gmtvec('gps_mod_rotation.gmtvec', lon, lat,
                         pred_rot[0::2], pred_rot[1::2],
                         'rigid block rotation contribution (mm/yr)')
            self._gmtvec('gps_mod_slip.gmtvec', lon, lat,
                         pred_slip[0::2], pred_slip[1::2],
                         'fault backslip / negative dislocation contribution '
                         '(mm/yr)')
            self._gmtvec('gps_mod_total.gmtvec', lon, lat,
                         total[0::2], total[1::2],
                         'rotation + backslip forward model (mm/yr)')
            # Four-column GMT readers use the first four columns; the two
            # weighted-residual columns are appended for diagnostics.
            with open(self._path('gps_residual.gmtvec'), 'w') as f:
                f.write('# observed - total residual (mm/yr) and weighted residual\n')
                f.write('# lon lat res_east res_north wres_east wres_north\n')
                np.savetxt(f, np.column_stack([
                    lon, lat, resid[0::2], resid[1::2],
                    wres[0::2], wres[1::2]]), fmt='%14.6f')

            self._gps_forward_table(lon, lat, stations, obs, sigma, pred_rot,
                                    pred_slip, total, resid, wres)
            self._fit_stats(obs, pred_rot + pred_slip, resid, wres,
                            sol, n_obs, linked, nf_mesh, nf_seg)

        if write_contrib:
            self._gps_contrib_euler(lon, lat, contrib_rot, euler_names, nblock)
            self._gps_contrib_patch(lon, lat, contrib_patch,
                                    patch_names[:nf_mesh], 'gps_contrib_fault.dat',
                                    'per-patch backslip contribution')
            if flt_seg is not None:
                self._gps_contrib_patch(lon, lat,
                                        contrib_patch[:, nf_mesh:],
                                        patch_names[nf_mesh:],
                                        'gps_contrib_fault_segment.dat',
                                        'per-patch segment backslip contribution')
            self._param_summary(omega, slip, chi, contrib_rot, contrib_slip,
                                chi_contrib, euler_names, slip_names,
                                patch_names, total, linked)
            self._write_contrib_h5(
                lon, lat, stations, obs, sigma, pred_rot, pred_slip, total,
                resid, wres, euler, omega, G_rot, contrib_rot, G_slip,
                contrib_slip, slip3, slip_names, patch_names,
                G_chi, chi_contrib, chi, linked, nf_mesh, nf_seg,
                float(sol.smo_facts[last]))

        # ------------------------------------------------------------- SAR
        if len(data.d_sar) and write_forward:
            self._sar_forward(sol, data, last)

        nfiles = ('gps_forward.dat, gps_mod_rotation.gmtvec, '
                  'gps_mod_slip.gmtvec, gps_mod_total.gmtvec, '
                  'gps_residual.gmtvec, gps_fit_stats.dat')
        if write_contrib:
            nfiles += (', gps_contrib_euler.dat, gps_contrib_block.gmtvec, '
                       'gps_contrib_fault.dat, param_contrib_summary.dat, '
                       'param_contrib_ranked.dat, contributions.h5')
        logging.info('wrote %s to %s', nfiles, os.path.abspath(self.outpath))

    # ---------------------------------------------------------- forward I/O
    def _euler_files(self, omega):
        """Raw buildG-convention vector, block table, and physical Euler pole."""
        euler = omega.reshape(-1, 3)
        np.savetxt(self._path('euler_est.dat'), omega.reshape(1, -1),
                   fmt='%22.12e',
                   header=('block-major Euler vectors [wx wy wz] per block; '
                           'v = B @ omega in the velocity units'),
                   comments='# ')
        with open(self._path('euler_block.dat'), 'w') as f:
            f.write('# block              wx                 wy                 wz\n')
            for b, row in enumerate(euler):
                f.write('%6d %18.10e %18.10e %18.10e\n'
                        % (b, row[0], row[1], row[2]))

        R_MM = 6.371e9                 # mean Earth radius, mm
        mag  = np.linalg.norm(euler, axis=1)
        safe = np.where(mag > 0, mag, 1.0)
        pole_lat = np.degrees(np.arcsin(np.clip(euler[:, 2] / safe, -1.0, 1.0)))
        pole_lon = np.degrees(np.arctan2(euler[:, 1], euler[:, 0]))
        rate_deg = np.degrees(mag / R_MM) * 1.0e6      # deg/Myr
        rate_mas = mag / R_MM * np.degrees(1.0) * 3600.0 * 1.0e3
        pole_lat[mag <= 0] = np.nan
        pole_lon[mag <= 0] = np.nan
        with open(self._path('euler_pole.dat'), 'w') as f:
            f.write('# Euler poles converted from omega under the buildG convention,\n'
                    '# assuming the input velocities are mm/yr (R = 6371 km):\n'
                    '# block  pole_lat_deg  pole_lon_deg  rate_deg_Myr  rate_mas_yr'
                    '              wx                 wy                 wz\n')
            for b in range(euler.shape[0]):
                f.write('%6d %14.6f %13.6f %13.6f %14.6f %18.10e %18.10e %18.10e\n'
                        % (b, pole_lat[b], pole_lon[b], rate_deg[b], rate_mas[b],
                           euler[b, 0], euler[b, 1], euler[b, 2]))

    def _gmtvec(self, name, lon, lat, east, north, title):
        with open(self._path(name), 'w') as f:
            f.write('# %s\n' % title)
            f.write('# lon lat east north\n')
            np.savetxt(f, np.column_stack([lon, lat, east, north]),
                       fmt='%14.6f')

    def _gps_forward_table(self, lon, lat, stations, obs, sigma, rot, slip,
                           total, resid, wres):
        header = ('# GPS forward decomposition (mm/yr)\n'
                  '# lon lat obs_E obs_N sig_E sig_N rot_E rot_N '
                  'slip_E slip_N total_E total_N res_E res_N wres_E wres_N station\n'
                  '# rot = rigid block rotation; slip = fault backslip '
                  '(negative dislocation);\n'
                  '# res = observed - total; wres = res / sigma\n')
        with open(self._path('gps_forward.dat'), 'w') as f:
            f.write(header)
            for i in range(len(lon)):
                vals = (obs[2*i], obs[2*i+1], sigma[2*i], sigma[2*i+1],
                        rot[2*i], rot[2*i+1], slip[2*i], slip[2*i+1],
                        total[2*i], total[2*i+1], resid[2*i], resid[2*i+1],
                        wres[2*i], wres[2*i+1])
                f.write(('%13.6f %13.6f' + ' %12.5f' * 14 + '  %s\n')
                        % ((lon[i], lat[i]) + vals + (stations[i],)))

    def _fit_stats(self, obs, pred, resid, wres, sol, n_obs, linked,
                   nf_mesh, nf_seg):
        def comp(sl):
            o, p, r = obs[sl], pred[sl], resid[sl]
            n = o.size
            ro = np.sqrt(np.mean(o*o)) if n else np.nan
            rp = np.sqrt(np.mean(p*p)) if n else np.nan
            rr = np.sqrt(np.mean(r*r)) if n else np.nan
            corr = (np.corrcoef(o, p)[0, 1]
                    if n > 1 and np.std(o) > 0 and np.std(p) > 0 else np.nan)
            return n, ro, rp, rr, np.mean(r), np.std(r), corr

        wrss = float(wres @ wres)
        n_fit = 3 * (len(sol.euler[-1]) // 3) + (nf_mesh if linked
                                                 else 3 * (nf_mesh + nf_seg))
        with open(self._path('gps_fit_stats.dat'), 'w') as f:
            f.write('# GPS fit statistics; velocities and residuals in mm/yr\n')
            f.write('# component       N    rms_obs   rms_pred  rms_resid'
                    '   mean_resid   std_resid        corr\n')
            for name, sl in (('east', slice(0, None, 2)),
                             ('north', slice(1, None, 2)),
                             ('both', slice(None))):
                n, ro, rp, rr, mr, sr, cc = comp(sl)
                f.write('%-10s %7d %10.4f %10.4f %10.4f %12.4f %11.4f %11.4f\n'
                        % (name, n, ro, rp, rr, mr, sr, cc))
            f.write('#\n')
            f.write('# WRSS (sum of weighted residual squared) : %.6f\n' % wrss)
            f.write('# WRSS / N                                : %.6f\n'
                    % (wrss / max(n_obs, 1)))
            f.write('# fitted parameters (omega, chi/slip)      : %d\n' % n_fit)
            f.write('# GPS components                          : %d\n' % n_obs)
            f.write('# note: the model is regularised, so WRSS/N is not a reduced\n'
                    '#       chi-square and N < n_fit is possible.\n')

    def _gps_contrib_euler(self, lon, lat, contrib, names, nblock):
        cols  = [lon, lat]
        label = ['lon', 'lat']
        for j, name in enumerate(names):
            cols.append(contrib[0::2, j])
            cols.append(contrib[1::2, j])
            label += [name + '_E', name + '_N']
        with open(self._path('gps_contrib_euler.dat'), 'w') as f:
            f.write('# per-Euler-parameter contribution at each GPS site (mm/yr)\n')
            f.write('# ' + ' '.join(label) + '\n')
            np.savetxt(f, np.column_stack(cols), fmt='%14.6f')

        # Aggregate the three components of each block into one E/N vector.
        cb = contrib.reshape(-1, nblock, 3).sum(axis=2)
        cols  = [lon, lat]
        label = ['lon', 'lat']
        for b in range(nblock):
            cols.append(cb[0::2, b])
            cols.append(cb[1::2, b])
            label += ['blk%02d_E' % b, 'blk%02d_N' % b]
        with open(self._path('gps_contrib_block.gmtvec'), 'w') as f:
            f.write('# rigid-rotation contribution of each block (sum of wx,wy,wz; mm/yr)\n')
            f.write('# ' + ' '.join(label) + '\n')
            np.savetxt(f, np.column_stack(cols), fmt='%14.6f')

    def _gps_contrib_patch(self, lon, lat, contrib, names, filename, title):
        cols  = [lon, lat]
        label = ['lon', 'lat']
        for p, name in enumerate(names):
            cols.append(contrib[0::2, p])
            cols.append(contrib[1::2, p])
            label += [name + '_E', name + '_N']
        with open(self._path(filename), 'w') as f:
            f.write('# %s (mm/yr; linked model: contribution of the fitted chi)\n'
                    % title)
            f.write('# ' + ' '.join(label) + '\n')
            np.savetxt(f, np.column_stack(cols), fmt='%14.6f')

    @staticmethod
    def _contrib_stats(c):
        e = c[0::2]
        n = c[1::2]
        return (np.sqrt(np.mean(e * e)), np.sqrt(np.mean(n * n)),
                np.max(np.abs(e)), np.max(np.abs(n)),
                np.sqrt(np.mean(e * e + n * n)))

    def _param_summary(self, omega, slip, chi, contrib_rot, contrib_slip,
                       chi_contrib, euler_names, slip_names, patch_names,
                       total, linked):
        model_mag = np.sqrt(np.mean(total[0::2] ** 2 + total[1::2] ** 2))
        rows = []
        for j, name in enumerate(euler_names):
            rows.append(('euler', name, omega[j], contrib_rot[:, j]))
        if linked:
            for p in range(chi.size):
                rows.append(('chi', patch_names[p], chi[p], chi_contrib[:, p]))
            # Physical slip components are derived from chi * s_lt, not
            # independent parameters, but their contributions are what the
            # backslip forward model is made of.
            for j, name in enumerate(slip_names):
                rows.append(('slip*', name, slip[j], contrib_slip[:, j]))
        else:
            for j, name in enumerate(slip_names):
                rows.append(('slip', name, slip[j], contrib_slip[:, j]))

        stats = [(kind, name, value, self._contrib_stats(c))
                 for kind, name, value, c in rows]
        intro = ('# Euler values use the buildG convention; chi is dimensionless;\n'
                 '# slip is mm/yr.  slip* = derived chi * s_lt (linked model).\n'
                 '# frac_model = rms_mag / RMS magnitude of the total model; it is\n'
                 '# a magnitude ratio, not a variance decomposition (contributions\n'
                 '# from different parameters are correlated).\n')
        columns = ('# kind   parameter              value        rms_E        rms_N'
                   '       max_E       max_N     rms_mag  frac_model\n')

        def write_rows(f, data):
            for kind, name, value, (re, rn, me, mn, rm) in data:
                frac = rm / model_mag if model_mag > 0 else np.nan
                f.write('%-6s %-20s %14.6g %12.5f %12.5f %11.5f %11.5f %11.5f %10.4f\n'
                        % (kind, name, value, re, rn, me, mn, rm, frac))

        with open(self._path('param_contrib_summary.dat'), 'w') as f:
            f.write('# Contribution of each model parameter to the GPS velocities '
                    '(mm/yr)\n')
            f.write(intro)
            f.write(columns)
            write_rows(f, stats)
        with open(self._path('param_contrib_ranked.dat'), 'w') as f:
            f.write('# Same table, sorted by rms_mag (largest surface effect first)\n')
            f.write(intro)
            f.write(columns)
            write_rows(f, sorted(stats, key=lambda r: r[3][4], reverse=True))

    def _write_contrib_h5(self, lon, lat, stations, obs, sigma, pred_rot,
                          pred_slip, total, resid, wres, euler, omega,
                          G_rot, contrib_rot, G_slip, contrib_slip, slip3,
                          slip_names, patch_names, G_chi, chi_contrib, chi,
                          linked, nf_mesh, nf_seg, smo_fact):
        try:
            import h5py
        except ImportError:
            logging.warning('h5py is not available; skipping contributions.h5')
            return
        n_gps = len(lon)
        n_obs = 2 * n_gps
        nf_total = nf_mesh + nf_seg
        dt = h5py.string_dtype(encoding='utf-8')
        with h5py.File(self._path('contributions.h5'), 'w') as h5:
            h5.create_dataset('station', data=np.asarray(stations, dtype=object),
                              dtype=dt)
            h5['lon'] = lon
            h5['lat'] = lat
            h5['obs'] = obs.reshape(n_gps, 2)
            h5['sigma'] = sigma.reshape(n_gps, 2)
            h5['pred_rotation'] = pred_rot.reshape(n_gps, 2)
            h5['pred_backslip'] = pred_slip.reshape(n_gps, 2)
            h5['pred_total'] = total.reshape(n_gps, 2)
            h5['residual'] = resid.reshape(n_gps, 2)
            h5['residual_weighted'] = wres.reshape(n_gps, 2)

            h5['euler'] = euler
            h5['omega_raw'] = omega
            h5.create_dataset('euler_names', data=np.asarray(
                ['blk%02d_%s' % (b, c)
                 for b in range(euler.shape[0]) for c in ('wx', 'wy', 'wz')],
                dtype=object), dtype=dt)
            h5['G_rotation'] = G_rot.reshape(n_gps, 2, -1)
            h5['contrib_rotation'] = contrib_rot.reshape(n_gps, 2, -1)

            h5['slip'] = slip3
            h5.create_dataset('slip_names',
                              data=np.asarray(slip_names, dtype=object), dtype=dt)
            h5['G_backslip'] = G_slip.reshape(n_gps, 2, -1)
            h5['contrib_backslip'] = contrib_slip.reshape(n_gps, 2, -1)
            h5['G_backslip_patch'] = G_slip.reshape(
                n_gps, 2, nf_total, 3).sum(axis=3)
            h5['contrib_backslip_patch'] = contrib_slip.reshape(
                n_gps, 2, nf_total, 3).sum(axis=3)
            h5.create_dataset('patch_names',
                              data=np.asarray(patch_names, dtype=object), dtype=dt)

            if linked:
                h5['chi'] = chi
                h5['G_chi'] = G_chi.reshape(n_gps, 2, nf_mesh)
                h5['contrib_chi'] = chi_contrib.reshape(n_gps, 2, nf_mesh)
                h5.create_dataset(
                    'chi_names',
                    data=np.asarray(['patch%03d' % p for p in range(nf_mesh)],
                                    dtype=object), dtype=dt)

            h5.attrs['units'] = ('velocities and contributions mm/yr; chi '
                                 'dimensionless')
            h5.attrs['station_order'] = 'block-ordered, same as gps_forward.dat'
            h5.attrs['model'] = 'linked' if linked else 'free-slip'
            h5.attrs['nblock'] = euler.shape[0]
            h5.attrs['nf_mesh'] = nf_mesh
            h5.attrs['nf_segment'] = nf_seg
            h5.attrs['smoothing_factor'] = smo_fact
            h5.attrs['array_layout'] = ('G_* and contrib_* are (station, '
                                        'component, parameter); component 0=E, 1=N')

    def _sar_forward(self, sol, data, last):
        obs   = np.asarray(data.d_sar, dtype=float).ravel()
        rot   = np.asarray(sol.sar_mod_rotation[last], dtype=float).ravel()
        slip  = np.asarray(sol.sar_mod_slip[last], dtype=float).ravel()
        total = rot + slip
        resid = obs - total
        llh   = np.asarray(data.llh_sar, dtype=float)
        with open(self._path('sar_forward.dat'), 'w') as f:
            f.write('# SAR forward decomposition (LOS units of the input data)\n')
            f.write('# lon lat obs rot slip total res\n')
            out = np.column_stack([llh[:, 1], llh[:, 0], obs, rot, slip,
                                   total, resid])
            np.savetxt(f, out, fmt='%14.6f')

    def _slipmodel(self, name, flt, slip_flat, cpl):
        """
        geom(7) + slip(3) + s_lt_ss + s_lt_ds + |s_lt| + chi.

        The s_lt/chi block is nan when the coupling stage did not run for this
        fault, so the column count never depends on the config.
        """
        slip = np.asarray(slip_flat, dtype=float).reshape(-1, 3)
        if cpl is None:
            s_lt_ss = np.full(flt.nf, np.nan)
            s_lt_ds = np.full(flt.nf, np.nan)
            s_lt    = np.full(flt.nf, np.nan)
            chi     = np.full(flt.nf, np.nan)
        else:
            s_lt_ss = cpl.s_lt_vec[:, 0]
            s_lt_ds = cpl.s_lt_vec[:, 1]
            s_lt    = cpl.s_lt
            chi     = cpl.chi

        out = np.column_stack([flt.geom_grid[:, 0:7], slip,
                               s_lt_ss, s_lt_ds, s_lt, chi])
        with open(self._path(name), 'w') as f:
            f.write(SLIPMODEL_HEADER)
            np.savetxt(f, out, fmt=14*' %10.3f')

    # ------------------------------------------------------------- coupling
    def _coupling(self, sol, blocks, flts, greens):
        """
        Linked model: chi was fitted, so reuse sol.coupling as-is.
        Free-slip fallback: derive chi as the projection of s onto s_lt.

        Returns (cpl_mesh, cpl_seg); the second is None when the model has no
        dict_segment.
        """
        c       = self.cfg
        last    = -1
        dict_coupling = cfg_section(c, 'dict_coupling', required=False) or {}
        min_slt = float(dict_coupling.get('min_slt', 1.0))
        exact   = bool(dict_coupling.get('require_exact', False))

        flt_mesh, flt_seg = flts
        green_mesh, green_seg = greens
        nf_mesh = flt_mesh.nf
        # sol.slip is [mesh_slip, segment_slip] for the Gen model, mesh only otherwise
        slip      = sol.slip[last]
        slip_mesh = slip[:3*nf_mesh]

        if hasattr(sol, 'coupling'):
            logging.info('coupling: linked model already fitted chi; reusing it')
            cpl_mesh = sol.coupling
            cpl_mesh.set_chi(sol.chi[last])
            cpl_mesh.slip = np.asarray(slip_mesh).reshape(-1, 3)
            # sign_agreement is tautological in linked mode (s = chi * s_lt
            # with chi >= 0), so report the chi distribution without it.
            cpl_mesh.report()
        else:
            cpl_mesh = self._one_coupling(blocks, flt_mesh, green_mesh,
                                          sol.euler[last], slip_mesh,
                                          min_slt, exact, dict_coupling)
            cpl_mesh.report(slip_mesh)
        cpl_mesh.save(self._path('coupling.faultgeom'))
        self._chi_profile(cpl_mesh, flt_mesh)

        cpl_seg = None
        if flt_seg is not None:
            cpl_seg = self._one_coupling(blocks, flt_seg, green_seg,
                                         sol.euler[last], slip[3*nf_mesh:],
                                         min_slt, exact, dict_coupling)
            logging.info('segment fault: s_lt = %.3f mm/yr, chi = %.3f',
                         cpl_seg.s_lt[0], cpl_seg.chi[0])

        logging.info('wrote coupling.faultgeom and chi_profile.dat to %s',
                     os.path.abspath(self.outpath))
        return cpl_mesh, cpl_seg

    @staticmethod
    def _one_coupling(blocks, flt, green, euler, slip, min_slt, exact,
                      dict_coupling):
        cpl = Coupling(blocks, flt, green=green, dict_coupling=dict_coupling)
        cpl.local_basis()
        cpl.long_term_slip(euler, require_exact=exact, min_slt=min_slt)
        cpl.chi_from_slip(slip)
        return cpl

    def _chi_profile(self, cpl, flt):
        prof  = cpl.depth_profile()
        depth = flt.geom_grid[:, 1].reshape(flt.ndeps, flt.nsegs)[:, 0]
        with open(self._path('chi_profile.dat'), 'w') as f:
            f.write('# depth_km  n_patch   chi_p25   chi_median  chi_p75\n')
            for iz in range(flt.ndeps):
                row = prof[iz][np.isfinite(prof[iz])]
                if row.size:
                    q25, q50, q75 = np.percentile(row, [25, 50, 75])
                    f.write('%9.3f %8d %10.3f %11.3f %10.3f\n'
                            % (depth[iz], row.size, q25, q50, q75))
                else:
                    f.write('%9.3f %8d %10s\n' % (depth[iz], 0, 'nan'))

        nseg = int(cfg_get(self.cfg, 'dict_coupling', 'nseg_mean', default=8))
        seg_mean, seg_edges = cpl.segment_mean(nseg=nseg)
        logging.info('segment-mean chi over %d along-strike groups:', nseg)
        for k in range(len(seg_mean)):
            logging.info('   columns %2d-%2d   chi = %s', seg_edges[k], seg_edges[k+1]-1,
                         '%6.2f' % seg_mean[k] if np.isfinite(seg_mean[k]) else '  n/a')

    # ------------------------------------------------------------------ view
    def _plot(self, sol, data, blocks):
        import matplotlib.pyplot as plt   # only pay for it when asked
        last = -1
        for k, val in ((1, sol.r[last]), (2, sol.dhat[last])):
            plt.subplot(1, 2, k)
            ax = plt.gca()
            blocks.gdf.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1)
            data.plot_gps_vector(ax, val, scale=self.vecscale, color='red')
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Block-model backslip interseismic inversion, driven by a YAML config file.')
    parser.add_argument('--cfgfile', type=str, required=True,
                        help='YAML config file; see block_backslip.yaml')
    parser.add_argument('--gpsfile', type=str, default='',
                        help='override dict_data.gpsfile')
    parser.add_argument('--outpath', type=str, default='',
                        help='override dict_export.outpath')
    parser.add_argument('--smoothfactor', type=float, nargs=3,
                        metavar=('START', 'STOP', 'STEP'), default=None,
                        help='override dict_weight.smoothfactor; evaluated as '
                             'START, START+STEP, ..., STOP (STOP included when '
                             'it lies on the grid)')
    parser.add_argument('--no-coupling', dest='no_coupling', action='store_true',
                        help='skip the locking-coefficient stage')
    args = parser.parse_args()

    try:
        BlockBackslipRun(args).run()
    except ConfigError as err:
        sys.exit('config error: %s' % err)


if __name__ == '__main__':
    main()
