#!/usr/bin/env python
"""
Block-model backslip interseismic inversion, driven by a YAML config file.

    python -m pyginvc.scripts.block_backslip_inversion --cfgfile block_backslip.yaml

The config picks the model.  A `dict_segment` section adds a second, unsmoothed
fault to the inversion and selects GenBlockBackslipInversion; without it the
meshed fault alone is inverted with BlockBackslipInversion.

Sections the driver reads:
    dict_data      gpsfile / sarfile / levfile / gfiletype
    dict_block     blockfile, assign
    dict_fault     faultfile, nsegs, ndeps, doSubFault, origin
    dict_segment   same shape as dict_fault; optional
    dict_lap       method
    dict_coupling  compute, min_slt, require_exact, nseg_mean
    dict_export    outpath, view, vecscale

Sections forwarded to the library verbatim (the driver never inspects them, so
library options can be used without touching this file):
    dict_green, dict_weight, dict_bound
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
from pyginvc.Inversion.BlockBackslipInversion import BlockBackslipInversion
from pyginvc.Inversion.Coupling import Coupling
from pyginvc.Inversion.GenBlockBackslipInversion import GenBlockBackslipInversion
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap

REQUIRED_SECTIONS = ['dict_data', 'dict_block', 'dict_fault',
                     'dict_green', 'dict_weight', 'dict_bound']

SLIPMODEL_HEADER = """\
# Block backslip inversion -- slip deficit model
#
#   1 width_km     2 depth_km    3 dip_deg      4 lat0        5 lon0
#   6 lat1         7 lon1
#   8 ss_mm_yr     9 ds_mm_yr   10 op_mm_yr
#  11 s_lt_mm_yr  12 chi
#
# ss/ds/op : fitted slip DEFICIT rates (backslip), not slip rates.
#
# s_lt     : long-term along-strike rate at the patch, from the rigid relative
#            rotation of the two blocks flanking it, using the same B and omega
#            the inversion fitted.  It is NOT an independent measurement -- it
#            inherits the omega <-> s trade-off, so chi is only as well
#            constrained as omega is.
#
# chi      : ss / s_lt, the locking coefficient.  DERIVED, not fitted, and NOT
#            clipped to [0, 1].  A value outside that range says the data want
#            more (or less) slip deficit than the block model's long-term rate
#            allows.
#
# nan in columns 11-12 means |s_lt| is below dict_coupling.min_slt, or no block
# pair was found for that patch, or the coupling stage was skipped.  It does NOT
# mean chi = 0.
#
# Note that ss is bounded by dict_bound.ss_range while chi is not: when many
# patches sit exactly on a bound, chi is partly a picture of where the bound
# bit, not only of the data.
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
        data   = self._data()
        blocks = self._blocks(data)
        flts, greens, lap = self._model(data)
        sol    = self._invert(data, blocks, flts, greens, lap)
        # coupling before export: slipmodel*.faultgeom carries s_lt and chi
        cpls   = self._coupling(sol, blocks, flts) if self.coupling_on else (None, None)
        self._export(sol, data, blocks, flts, cpls)
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
        dict_weight = cfg_section(c, 'dict_weight')
        dict_bound  = cfg_section(c, 'dict_bound')

        G_gps_block = block_diag(*blocks.buildG_block(data.ndim, block_type='gps'))
        nblock      = blocks.nblock

        if flt_seg is None:
            logging.info('model: mesh only (%d patches), %d blocks', flt_mesh.nf, nblock)
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

    def _export(self, sol, data, blocks, flts, cpls):
        """
        Write the last smoothing factor's model.

        sol.slip / sol.euler / sol.x / sol.r / sol.dhat / sol.gps_mod_* are all
        (nsmooth, ...); [-1] is the last smoothing factor.

        The slip vector for the Gen model is [mesh_slip, segment_slip], so the
        mesh slice has to be taken explicitly -- slicing to the end and
        reshaping would give flt_mesh.nf + flt_seg.nf rows and column_stack
        against the mesh geometry would fail.
        """
        last     = -1
        flt_mesh = flts[0]
        nblock   = blocks.nblock
        nf_mesh  = flt_mesh.nf

        if not os.path.isdir(self.outpath):
            os.makedirs(self.outpath)

        np.savetxt(self._path('euler_est.dat'), sol.euler)

        outdata = np.column_stack([data.llh_gps[:, 1], data.llh_gps[:, 0],
                                   sol.gps_mod_slip[last].reshape(-1, 2)])
        np.savetxt(self._path('gps_mod_slip.gmtvec'), outdata)

        self._slipmodel('slipmodel.faultgeom', flt_mesh,
                        sol.x[last, 3*nblock:3*nblock + 3*nf_mesh], cpls[0])

        if flts[1] is not None:
            self._slipmodel('slipmodel_segment.faultgeom', flts[1],
                            sol.x[last, 3*nblock + 3*nf_mesh:], cpls[1])

        logging.info('wrote euler_est.dat, gps_mod_slip.gmtvec, slipmodel.faultgeom%s to %s',
                     ' and slipmodel_segment.faultgeom' if flts[1] is not None else '',
                     os.path.abspath(self.outpath))

    def _slipmodel(self, name, flt, slip_flat, cpl):
        """
        geom(7) + slip(3) + s_lt + chi.

        s_lt and chi are nan when the coupling stage did not run for this fault,
        so the column count never depends on the config.
        """
        slip = np.asarray(slip_flat, dtype=float).reshape(-1, 3)
        if cpl is None:
            s_lt = np.full(flt.nf, np.nan)
            chi  = np.full(flt.nf, np.nan)
        else:
            s_lt, chi = cpl.s_lt, cpl.chi

        out = np.column_stack([flt.geom_grid[:, 0:7], slip, s_lt, chi])
        with open(self._path(name), 'w') as f:
            f.write(SLIPMODEL_HEADER)
            np.savetxt(f, out, fmt=12*' %10.3f')

    # ------------------------------------------------------------- coupling
    def _coupling(self, sol, blocks, flts):
        """
        Mode A: chi = s_ss / s_lt, free, not clipped to [0, 1].

        Returns (cpl_mesh, cpl_seg); the second is None when the model has no
        dict_segment.
        """
        c       = self.cfg
        last    = -1
        min_slt = float(cfg_get(c, 'dict_coupling', 'min_slt', default=1.0))
        exact   = bool(cfg_get(c, 'dict_coupling', 'require_exact', default=False))

        flt_mesh, flt_seg = flts
        nf_mesh = flt_mesh.nf
        # sol.slip is [mesh_slip, segment_slip] for the Gen model, mesh only otherwise
        slip      = sol.slip[last]
        slip_mesh = slip[:3*nf_mesh]

        cpl_mesh = self._one_coupling(blocks, flt_mesh, sol.euler[last], slip_mesh,
                                      min_slt, exact)
        cpl_mesh.report(slip_mesh)
        cpl_mesh.save(self._path('coupling.faultgeom'))
        self._chi_profile(cpl_mesh, flt_mesh)

        cpl_seg = None
        if flt_seg is not None:
            cpl_seg = self._one_coupling(blocks, flt_seg, sol.euler[last],
                                         slip[3*nf_mesh:], min_slt, exact)
            logging.info('segment fault: s_lt = %.3f mm/yr, chi = %.3f',
                         cpl_seg.s_lt[0], cpl_seg.chi[0])

        logging.info('wrote coupling.faultgeom and chi_profile.dat to %s',
                     os.path.abspath(self.outpath))
        return cpl_mesh, cpl_seg

    @staticmethod
    def _one_coupling(blocks, flt, euler, slip, min_slt, exact):
        cpl = Coupling(blocks, flt)
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
                             'np.arange(START, STOP, STEP), STOP excluded')
    parser.add_argument('--no-coupling', dest='no_coupling', action='store_true',
                        help='skip the locking-coefficient stage')
    args = parser.parse_args()

    try:
        BlockBackslipRun(args).run()
    except ConfigError as err:
        sys.exit('config error: %s' % err)


if __name__ == '__main__':
    main()
