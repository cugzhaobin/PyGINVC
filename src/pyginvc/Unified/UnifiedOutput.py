#!/usr/bin/env python
"""Output writer for mixed rectangular/triangular slip models.

The existing ``Export.Output`` writer stores one flat fault type under
``/flt``.  That cannot represent a model containing both a rectangle and a
triangle mesh, so this class writes an explicit per-fault HDF5 layout:

``/flt/schema_version = 2``
``/flt/origin``
``/flt/faults/<NN>/{kind,name,nf,slip_offset,rake_beta,...}``

Each fault group carries either ``dis_geom``/``geom`` (rectangle) or
``vertex_enu``/``vertex_llh``/``element`` (triangle).  Observations and
solutions use the same column order as the legacy files.

A homogeneous model can additionally be written in the legacy flat layout
by setting ``export.legacy_h5: true``.  The legacy file is named
``solutions_legacy.h5`` so that a mixed run can never be mistaken for a
legacy one.
"""

import os
import logging
import tempfile

import h5py
import numpy as np

from pyginvc.Export import OutputUtils as ou

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(filename)s:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")


class UnifiedOutput(object):
    """Write models, observations and reports for a unified inversion.

    Parameters
    ----------
    data : GeoDataCollection
        Loaded observation collection.
    faults : UnifiedFaultCollection
        Loaded rectangular/triangular fault collection.
    sol : object, optional
        Solution object exposing ``slip``, ``dhat``, ``r``, ``misfit`` and
        the smoothing metadata.
    outdir : str
        Output directory.
    config, export : dict, optional
        Configuration sections used for optional switches.
    """

    def __init__(self, data, faults, sol=None, outdir='.',
                 config=None, export=None):
        self.data = data
        self.faults = faults
        self.sol = sol
        self.outdir = outdir or '.'
        self.config = config or {}
        self.export = export or {}
        os.makedirs(self.outdir, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Observation metadata
    # ------------------------------------------------------------------ #
    @staticmethod
    def _gps_sigma(gps, n, ndim):
        """Per-station formal uncertainty from a GPS covariance array."""
        cov = getattr(gps, 'cov_gps', None)
        if cov is None or np.size(cov) == 0:
            return None
        cov = np.asarray(cov)
        var = cov if cov.ndim == 1 else np.diag(cov)
        return np.sqrt(var).reshape(n, ndim)

    @staticmethod
    def _file_names(data, kind):
        """Base names of the source files for one observation kind."""
        objects = getattr(data, {'gps': 'gps_data', 'lev': 'lev_data',
                                 'sar': 'sar_data'}[kind], None) or []
        attr = {'gps': 'gpsfile', 'lev': 'levfile', 'sar': 'sarfile'}[kind]
        return [os.path.basename(str(getattr(obj, attr, '') or ''))
                for obj in objects]

    def _obs_dict(self):
        data = self.data
        obs = {
            'llh_gps': data.llh_gps,
            'llh_lev': data.llh_lev,
            'llh_sar': data.llh_sar,
            'unit': data.unit,
            'd_gps': data.d_gps,
            'd_lev': data.d_lev,
            'd_sar': data.d_sar,
            'ndim': int(data.ndim),
            'n_sar': np.asarray(data.n_sar, dtype=int),
            'n_gps': np.asarray(data.n_gps, dtype=int),
            'n_lev': np.asarray(data.n_lev, dtype=int),
            'gps_files': self._file_names(data, 'gps'),
            'lev_files': self._file_names(data, 'lev'),
            'sar_files': self._file_names(data, 'sar'),
        }
        if self.sol is not None:
            ramp = getattr(self.sol, 'ramp_sar', None)
            flags = getattr(self.sol, 'sar_ramp_flag', None)
            if ramp is not None and np.size(ramp) > 0:
                obs['d_sar_ramp'] = np.asarray(ramp)
            if flags is not None and np.size(flags) > 0:
                obs['sar_ramp'] = np.asarray(flags, dtype=bool)
        return obs

    # ------------------------------------------------------------------ #
    # HDF5
    # ------------------------------------------------------------------ #
    def WriteHDF5(self, filename='solutions.h5'):
        """Write the v2 mixed-fault HDF5 file and optional legacy file."""
        if self.sol is None:
            raise ValueError('No solution is available for HDF5 output.')

        path = os.path.join(self.outdir, filename)
        self._write_h5_v2(path)
        logging.info('Unified solution is written to %s.', path)

        if bool(self.export.get('legacy_h5', False)):
            legacy = os.path.join(self.outdir, 'solutions_legacy.h5')
            if self._write_legacy_h5(legacy):
                logging.info('Legacy-compatible solution is written to %s.',
                             legacy)
        return path

    def _write_h5_v2(self, path):
        obs = self._obs_dict()
        sol = self.sol
        faults = self.faults

        outdir = os.path.dirname(os.path.abspath(path))
        fd, tmpfile = tempfile.mkstemp(prefix='.solutions_', suffix='.h5',
                                       dir=outdir)
        os.close(fd)
        try:
            with h5py.File(tmpfile, 'w') as h5:
                grp = h5.create_group('obs')
                for key in ('llh_gps', 'llh_lev', 'llh_sar', 'unit',
                            'd_gps', 'd_lev', 'd_sar'):
                    grp.create_dataset(key, data=obs[key], compression='gzip')
                grp['ndim'] = obs['ndim']
                for key in ('n_sar', 'n_gps', 'n_lev'):
                    if np.size(obs[key]):
                        grp.create_dataset(key, data=obs[key])
                name_dtype = h5py.string_dtype(encoding='utf-8')
                for key in ('gps_files', 'lev_files', 'sar_files'):
                    if obs[key]:
                        grp.create_dataset(key, data=obs[key], dtype=name_dtype)
                if obs.get('d_sar_ramp') is not None:
                    grp.create_dataset('d_sar_ramp',
                                       data=obs['d_sar_ramp'],
                                       compression='gzip')
                if obs.get('sar_ramp') is not None:
                    grp.create_dataset('sar_ramp', data=obs['sar_ramp'])

                grp = h5.create_group('flt')
                grp.create_dataset('schema_version', data=2)
                grp.create_dataset('origin', data=np.asarray(faults.origin))
                kinds = faults.kinds()
                flttype = kinds[0] if len(set(kinds)) == 1 else 'mixed'
                grp.create_dataset('flttype', data=flttype)
                grp.create_dataset('fault_count', data=len(kinds))
                fgrp = grp.create_group('faults')

                for i, flt, cfg, kind, name in faults.groups():
                    fg = fgrp.create_group('%02d' % i)
                    fg.create_dataset('kind', data=kind)
                    fg.create_dataset('name', data=name, dtype=name_dtype)
                    fg.create_dataset('nf', data=int(flt.nf))
                    fg.create_dataset('slip_offset',
                                      data=int(faults.offsets[i]))
                    beta = float(faults.rake_betas()[i])
                    fg.create_dataset('rake_beta', data=beta)
                    fg.create_dataset(
                        'smoothing_weight',
                        data=float((cfg or {}).get('smoothing_weight', 1.0)))
                    if kind == 'rectangle':
                        fg.create_dataset('dis_geom',
                                          data=np.asarray(flt.dis_geom_grid),
                                          compression='gzip')
                        fg.create_dataset('geom',
                                          data=np.asarray(flt.geom_grid),
                                          compression='gzip')
                    else:
                        fg.create_dataset('vertex_enu',
                                          data=np.asarray(flt.vertex_enu),
                                          compression='gzip')
                        fg.create_dataset('vertex_llh',
                                          data=np.asarray(flt.vertex_llh),
                                          compression='gzip')
                        fg.create_dataset('element',
                                          data=np.asarray(flt.element),
                                          compression='gzip')

                grp = h5.create_group('sol')
                grp.create_dataset('slip', data=np.asarray(sol.slip),
                                   compression='gzip')
                grp.create_dataset('r', data=np.asarray(sol.r),
                                   compression='gzip')
                grp.create_dataset('misfit', data=np.asarray(sol.misfit),
                                   compression='gzip')
                grp.create_dataset('dhat', data=np.asarray(sol.dhat),
                                   compression='gzip')
                grp.create_dataset('smo_facts',
                                   data=np.asarray(sol.smo_facts),
                                   compression='gzip')
                grp.create_dataset('ruff', data=np.asarray(sol.ruff),
                                   compression='gzip')
                grp.create_dataset('moment', data=np.asarray(sol.moment),
                                   compression='gzip')
                nmodels = np.asarray(sol.slip).shape[0]
                grp.create_dataset(
                    'chosen_index',
                    data=np.array([max(nmodels - 1, 0)], dtype=int))
                ramp = getattr(sol, 'ramp_sar', None)
                if ramp is not None and np.size(ramp) > 0:
                    grp.create_dataset('ramp_sar', data=np.asarray(ramp),
                                       compression='gzip')
                flags = getattr(sol, 'sar_ramp_flag', None)
                if flags is not None and np.size(flags) > 0:
                    grp.create_dataset('sar_ramp_flag',
                                       data=np.asarray(flags, dtype=bool))
            os.replace(tmpfile, path)
        finally:
            if os.path.exists(tmpfile):
                os.remove(tmpfile)

    def _write_legacy_h5(self, path):
        """Write the flat legacy layout when the model is representable."""
        kinds = self.faults.kinds()
        unique = set(kinds)
        groups = list(self.faults.groups())
        if len(unique) > 1:
            logging.warning(
                'legacy_h5 requested for a mixed model; the legacy format '
                'cannot represent it, so only solutions.h5 is written.')
            return False
        if unique == {'triangle'} and len(groups) != 1:
            logging.warning(
                'legacy_h5 requested for several triangle meshes; the legacy '
                'format stores only one mesh, so it is skipped.')
            return False

        if unique == {'rectangle'}:
            fault = {
                'origin': np.asarray(self.faults.origin),
                'flttype': 'rectangle',
                'dis_geom': np.vstack([np.asarray(flt.dis_geom_grid)
                                       for _, flt, _, _, _ in groups]),
                'geom': np.vstack([np.asarray(flt.geom_grid)
                                   for _, flt, _, _, _ in groups]),
            }
        else:
            _, flt, _, _, _ = groups[0]
            fault = {
                'origin': np.asarray(self.faults.origin),
                'flttype': 'triangle',
                'vertex_enu': np.asarray(flt.vertex_enu),
                'vertex_llh': np.asarray(flt.vertex_llh),
                'element': np.asarray(flt.element),
            }

        sol = self.sol
        solution = {
            'slip': np.asarray(sol.slip),
            'r': np.asarray(sol.r),
            'misfit': np.asarray(sol.misfit),
            'dhat': np.asarray(sol.dhat),
            'smo_facts': np.asarray(sol.smo_facts),
            'ruff': np.asarray(sol.ruff),
            'ramp_sar': np.asarray(getattr(sol, 'ramp_sar', [])),
            'sar_ramp_flag': np.asarray(
                getattr(sol, 'sar_ramp_flag', []), dtype=bool),
            'chosen_index': np.array(
                [max(np.asarray(sol.slip).shape[0] - 1, 0)], dtype=int),
        }
        ou.write_solution_h5(path, self._obs_dict(), fault, solution)
        return True

    # ------------------------------------------------------------------ #
    # Slip models
    # ------------------------------------------------------------------ #
    def WriteSlipModel(self, slip=None, k=None, outdir=None):
        """Write per-fault slip model files.

        ``slip`` may be a single flattened model or a ``(nmodel, 3*nf)``
        array.  ``k`` selects one zero-based model index.
        """
        if slip is None:
            if self.sol is None:
                raise ValueError('No slip given and the solution has none.')
            slip = self.sol.slip
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        arr = np.asarray(slip, dtype=float)
        single = arr.ndim == 1
        if single:
            arr = arr[None, :]
        if arr.shape[1] != self.faults.nf3:
            raise ValueError('slip has %d values, expected %d.'
                             % (arr.shape[1], self.faults.nf3))

        if k is not None:
            if k < 0:
                raise IndexError('model index %d is out of range.' % k)
            if not single and k >= arr.shape[0]:
                raise IndexError('model index %d is out of range.' % k)
            indices = [(0, k)] if single else [(k, k)]
        else:
            indices = [(row, row) for row in range(arr.shape[0])]

        written = []
        for row, label in indices:
            for (i, flt, cfg, kind, name), (_, _, block) in zip(
                    self.faults.groups(),
                    self.sol.slip_by_fault(arr[row])):
                prefix = 'slipmodel_%s_%03d' % (name, label + 1)
                if kind == 'rectangle':
                    path = os.path.join(prefix + '.faultgeom')
                    path = os.path.join(outdir, path)
                    ou.write_faultgeom(path, np.asarray(flt.geom_grid), block)
                else:
                    path = os.path.join(outdir, prefix + '.tri')
                    ou.write_tri_model(path, np.asarray(flt.element), block)
                written.append(path)

            self._write_combined_model(outdir, label, arr[row])
        return written

    def _write_combined_model(self, outdir, row, slip):
        """Write one long table with all faults' slip values."""
        path = os.path.join(outdir, 'slipmodel_all_%03d.tsv' % (row + 1))
        lines = ['# fault\tkind\tpatch\tstrike_slip_mm\tdip_slip_mm\t'
                 'opening_mm\ttotal_mm']
        for (i, flt, cfg, kind, name), (_, _, block) in zip(
                self.faults.groups(), self.sol.slip_by_fault(slip)):
            mag = np.linalg.norm(block[:, :2], axis=1)
            for j in range(int(flt.nf)):
                lines.append('%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f'
                             % (name, kind, j + 1, block[j, 0], block[j, 1],
                                block[j, 2], mag[j]))
        with open(path, 'w') as fid:
            fid.write('\n'.join(lines) + '\n')
        return path

    def WriteGMTSlip(self, slip=None, outdir=None, model_index=-1):
        """Write ss/ds/total slip polygons in GMT multi-segment format."""
        if slip is None:
            if self.sol is None:
                raise ValueError('No slip given and the solution has none.')
            slip = self.sol.slip
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        arr = np.asarray(slip, dtype=float)
        if arr.ndim == 1:
            row = arr
        else:
            row = arr[model_index]

        polygons, values = [], []
        for (i, flt, cfg, kind, name), (_, _, block) in zip(
                self.faults.groups(), self.sol.slip_by_fault(row)):
            vertices = self._vertices(flt, kind)
            # Rectangles have four corners and triangles three, so the
            # combined file cannot be a single (nf, ncorner, 3) array.
            polygons.extend(list(vertices))
            values.append(block)
        slip3 = np.vstack(values)
        slip_mag = np.linalg.norm(slip3[:, :2], axis=1)

        header = '# Unified fault slip model'
        fmt = '%10.4f %10.4f %10.4f'
        written = []
        for fname, value in (('ss_slip_llh.gmtlin', slip3[:, 0]),
                             ('ds_slip_llh.gmtlin', slip3[:, 1]),
                             ('ts_slip_llh.gmtlin', slip_mag)):
            path = os.path.join(outdir, fname)
            self._write_gmt_variable(path, polygons, value, header, fmt)
            written.append(path)
        return written

    @staticmethod
    def _write_gmt_variable(path, polygons, values, header, fmt):
        """Write GMT multi-segment polygons with a varying corner count."""
        values = np.asarray(values).reshape(-1)
        lines = [header]
        for polygon, value in zip(polygons, values):
            lines.append('>-Z%.3f' % value)
            for vertex in np.asarray(polygon):
                lines.append(fmt % (vertex[0], vertex[1], vertex[2]))
        np.savetxt(path, lines, fmt='%s')
        return path

    @staticmethod
    def _vertices(flt, kind):
        """Return ``(nf, ncorner, 3)`` ``[lon, lat, depth]`` patch corners."""
        if kind == 'triangle':
            v = np.asarray(flt.vertex_llh)[:, [1, 0, 2]]
            return v[np.asarray(flt.element, dtype=int) - 1]

        felem = flt.FaultGeom2AllVertex()
        corners = [('lb_lon', 'lb_lat', 'lb_dep'),
                   ('rb_lon', 'rb_lat', 'rb_dep'),
                   ('rt_lon', 'rt_lat', 'rt_dep'),
                   ('lt_lon', 'lt_lat', 'lt_dep')]
        return np.stack(
            [np.column_stack([np.asarray(felem[a]), np.asarray(felem[b]),
                              np.asarray(felem[c])])
             for a, b, c in corners], axis=1)

    # ------------------------------------------------------------------ #
    # Observed / modeled / residual data
    # ------------------------------------------------------------------ #
    def WriteData(self, dhat=None, outdir=None, ramp=None, parts=None):
        """Write per-file observed, modeled and residual data."""
        if dhat is None:
            if self.sol is None or self.sol.dhat is None:
                raise ValueError('No dhat given and the solution has none.')
            dhat = self.sol.dhat
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)
        dhat = np.asarray(dhat, dtype=float).reshape(-1)
        data = self.data
        ndim = int(data.ndim)
        parts = None if parts is None else {str(p) for p in parts}

        written = []
        offset = 0
        for i, (gps, _) in enumerate(data.gps_pairs()):
            n = len(gps.llh_gps)
            if n == 0:
                continue
            ndata = n * ndim
            obs = np.asarray(gps.d_gps).reshape(n, ndim)
            mod = dhat[offset:offset + ndata].reshape(n, ndim)
            if parts is None or parts & {'obs', 'mod', 'res'}:
                written += ou.write_gps_data(
                    gps.llh_gps, obs, mod,
                    sigma=self._gps_sigma(gps, n, ndim),
                    station=getattr(gps, 'station_gps', None),
                    outdir=outdir, prefix='gps%d' % (i + 1), parts=parts)
            offset += ndata

        levs = list(data.lev_pairs())
        for i, (lev, _) in enumerate(levs):
            n = len(np.asarray(lev.llh_lev))
            if n == 0:
                continue
            obs = np.asarray(lev.d_lev).reshape(-1)
            mod = dhat[offset:offset + n].reshape(-1)
            prefix = 'lev%d' % (i + 1) if len(levs) > 1 else 'lev'
            written += ou.write_lev_data(
                lev.llh_lev, obs, mod, outdir=outdir, prefix=prefix,
                parts=parts)
            offset += n

        n_sar = len(np.asarray(data.llh_sar))
        if n_sar:
            obs = np.asarray(data.d_sar).reshape(-1)
            mod = dhat[offset:offset + n_sar].reshape(-1)
            written += ou.write_sar_data(
                data.llh_sar, obs, mod, n_sar=data.n_sar,
                outdir=outdir, prefix='sar', parts=parts)
            if parts is None or any(p in ('mod', 'res') for p in parts):
                if ramp is None and self.sol is not None:
                    ramp = getattr(self.sol, 'ramp_sar', None)
                if ramp is not None and np.size(ramp) > 0:
                    flags = getattr(self.sol, 'sar_ramp_flag', None)
                    written += ou.write_sar_ramp_data(
                        data.llh_sar, ramp, n_sar=data.n_sar,
                        outdir=outdir, prefix='sar_ramp', flags=flags)
        return written

    def WriteSmoothing(self, k, result, outdir=None):
        """Write the model and fit for one smoothing factor."""
        outdir = outdir or self.outdir
        self.WriteSlipModel(result['slip'], k=k, outdir=outdir)
        if k == 0:
            self.WriteData(result['dhat'], outdir=outdir, parts=('obs',))

        subdir = os.path.join(outdir, 'sm%03d' % k)
        os.makedirs(subdir, exist_ok=True)
        self.WriteData(result['dhat'], outdir=subdir,
                       ramp=result.get('ramp_sar'), parts=('mod', 'res'))
        return subdir

    # ------------------------------------------------------------------ #
    # Summary / L-curve / archive
    # ------------------------------------------------------------------ #
    def WriteSummary(self, filename='summary'):
        """Write a summary that includes every fault's kind and size."""
        if self.sol is None:
            raise ValueError('No solution is available for summary output.')
        path = os.path.join(self.outdir, filename)
        sol = self.sol
        data = self.data
        faults = self.faults
        ou.write_summary(
            path,
            ndim=int(data.ndim),
            n_gps=int(np.sum(data.n_gps)),
            n_sar=int(np.sum(data.n_sar)),
            n_lev=int(np.sum(data.n_lev)),
            flttype='unified',
            nstrike=[], ndip=[],
            smo_facts=np.asarray(sol.smo_facts),
            misfit_gps=np.asarray(sol.misfit_gps),
            misfit_sar=np.asarray(sol.misfit_sar),
            n_d_gps=max(len(np.asarray(data.d_gps).reshape(-1)), 1),
            n_d_sar=max(len(np.asarray(data.d_sar).reshape(-1)), 1),
            moment=np.asarray(sol.moment))

        slip_last = np.asarray(sol.slip)
        slip_last = slip_last[-1] if slip_last.ndim > 1 else slip_last
        rows = faults.summary_rows(slip_last)
        with open(path, 'a') as fid:
            fid.write('\n#*#*#*#*#*#*#*#*#Unified Fault Model#*#*#*#*#*#*#*#*\n')
            fid.write('%-6s %-16s %-10s %8s %10s %10s %12s\n'
                      % ('Fault', 'Name', 'Kind', 'Patches', 'AlongStr',
                         'AlongDip', 'MeanSlip(mm)'))
            for row in rows:
                fid.write('%-6d %-16s %-10s %8d %10d %10d %12.4f\n'
                          % (row['index'] + 1, row['name'], row['kind'],
                             row['nf'], row['nstrike'], row['ndip'],
                             row.get('mean_slip', 0.0)))
        logging.info('Summary is written to %s.', path)
        return path

    def WriteLCurve(self, filename='ruff_misfit'):
        """Write the roughness/misfit trade-off table."""
        if self.sol is None:
            raise ValueError('No solution is available for L-curve output.')
        path = os.path.join(self.outdir, filename)
        ruff = np.asarray(self.sol.ruff, dtype=float).reshape(-1)
        misfit = np.asarray(self.sol.misfit, dtype=float).reshape(-1)
        smo = np.asarray(self.sol.smo_facts, dtype=float).reshape(-1)
        if smo.size < 2:
            return False
        moment = np.asarray(self.sol.moment, dtype=float)
        if moment.size == 2 * smo.size:
            moment = moment.reshape(smo.size, 2)
            mo, mw = moment[:, 0], moment[:, 1]
        else:
            mo = moment.reshape(-1)
            mw = np.full(smo.size, np.nan)
        # Keep the legacy column order so Export.View can read this file.
        table = np.column_stack([ruff, misfit, smo, mo, mw])
        np.savetxt(
            path, table,
            header=' Ruffness\t Misfit\t Smoothing_factor\t Mo\t Mw',
            fmt='%14.6e\t%14.6e\t%10.4f\t%12.3e\t%10.3f')
        return True

    # ------------------------------------------------------------------ #
    # HTML report
    # ------------------------------------------------------------------ #
    FIGURE_PATTERNS = (
        ('GPS observed / modeled', 'gps_mod_obs*.png'),
        ('InSAR observed / modeled', 'sar_mod_obs*.png'),
        # Plain lon/lat fallback used when cartopy is unavailable.
        ('Observed / modeled data (fallback)', 'obs_mod*.png'),
        ('InSAR ramp', 'sar_ramp*.png'),
        ('3D slip distribution', 'slip_3d*.png'),
        ('Trade-off (L-) curve', 'lcurve*.png'),
    )

    DOWNLOAD_PATTERNS = (
        'solutions.h5', 'solutions_legacy.h5', 'summary', 'ruff_misfit',
        'slipmodel*.faultgeom', 'slipmodel*.tri', 'slipmodel_all_*.tsv',
        'ss_slip_llh.gmtlin', 'ds_slip_llh.gmtlin', 'ts_slip_llh.gmtlin',
        'gps*.gmtvec', 'sar*.xyz', 'lev*.xyz',
    )

    def WriteReport(self, filename='report.html', embed=False, cfgfile=None,
                    title='PyGINVC report'):
        """Write an HTML summary report (metrics, tables, figures, downloads).

        Parameters
        ----------
        filename : str
            Output HTML file name inside ``outdir``.
        embed : bool
            Inline images as base64 for a single portable file.
        cfgfile : str, optional
            YAML config file shown in an appendix.
        title : str
            Page title.

        Returns
        -------
        str or None
            Path of the written HTML file, or None when it could not be built.
        """
        from datetime import datetime
        from pyginvc.Export.Report import HtmlReport

        try:
            os.makedirs(self.outdir, exist_ok=True)
            data = self.data
            faults = self.faults
            kinds = faults.kinds()

            n_gps = int(np.sum(np.asarray(data.n_gps, dtype=int)))
            n_sar = int(np.sum(np.asarray(data.n_sar, dtype=int)))
            n_lev = int(np.sum(np.asarray(data.n_lev, dtype=int)))

            rep = HtmlReport(self.outdir, title=title, embed=embed)
            rep.set_subtitle('Run finished %s'
                             % datetime.now().strftime('%Y-%m-%d %H:%M'))
            rep.add_meta('Fault type',
                         kinds[0] if len(set(kinds)) == 1 else 'mixed')
            rep.add_meta('Faults', len(kinds))
            rep.add_metric('GPS stations', n_gps)
            rep.add_metric('InSAR points', n_sar)
            rep.add_metric('Leveling points', n_lev)
            rep.add_metric('Fault patches', faults.nf)

            sol = self.sol
            if sol is not None:
                smo = np.asarray(sol.smo_facts).reshape(-1)
                mgps = np.asarray(sol.misfit_gps).reshape(-1)
                msar = np.asarray(sol.misfit_sar).reshape(-1)
                ruff = np.asarray(sol.ruff).reshape(-1)
                moment = np.asarray(sol.moment, dtype=float)
                if moment.size == 2 * smo.size:
                    moment = moment.reshape(smo.size, 2)

                headers = ['Smoothing factor']
                if n_gps:
                    headers.append('WRSS GPS')
                if n_sar:
                    headers.append('WRSS InSAR')
                if ruff.size:
                    headers.append('Roughness')
                if moment.ndim == 2 and moment.shape[0] == smo.size:
                    headers += ['Mo (N m)', 'Mw']

                rows = []
                for i in range(smo.size):
                    row = [float(smo[i])]
                    if n_gps:
                        row.append(float(mgps[i]) if mgps.size else '')
                    if n_sar:
                        row.append(float(msar[i]) if msar.size else '')
                    if ruff.size:
                        row.append(float(ruff[i]))
                    if moment.ndim == 2 and moment.shape[0] == smo.size:
                        row += [float(moment[i, 0]), float(moment[i, 1])]
                    rows.append(row)
                if rows:
                    rep.add_table('Inversion summary', headers, rows,
                                  note='One row per tested smoothing factor.')

                fit_rows = self._fit_rows(sol)
                if fit_rows:
                    rep.add_table(
                        'Data fit by file',
                        ['Data file', 'Type', 'N', 'WRMS (mm)',
                         'WRSS/N', 'Variance reduction'],
                        fit_rows,
                        note='WRMS and variance reduction use unweighted '
                             'residuals; WRSS/N uses the inversion-weighted '
                             'residual.')

                slip_last = np.asarray(sol.slip)
                slip_last = (slip_last[-1] if slip_last.ndim > 1 else slip_last)
                fault_rows = [
                    [r['name'], r['kind'], r['nf'], r.get('mean_slip', 0.0)]
                    for r in faults.summary_rows(slip_last)]
                if fault_rows:
                    rep.add_table('Fault model',
                                  ['Fault', 'Kind', 'Patches', 'Mean slip (mm)'],
                                  fault_rows)

            rep.add_figures_from_dir(patterns=self.FIGURE_PATTERNS)
            rep.add_downloads_from_dir(patterns=self.DOWNLOAD_PATTERNS)
            if cfgfile and os.path.isfile(cfgfile):
                with open(cfgfile, 'r', encoding='utf-8') as fid:
                    rep.add_config(fid.read())
            return rep.write(filename)
        except Exception as exc:      # never fail an inversion over the report
            logging.warning('HTML report could not be written: %s', exc)
            return None

    def _fit_rows(self, sol):
        """Per-file fit table with variance reduction formatted as percent."""
        try:
            rows = ou.fit_table(
                self.data.d_gps, int(self.data.ndim),
                self.data.d_lev, self.data.d_sar, sol.dhat,
                residual=getattr(sol, 'r', None),
                n_gps=self.data.n_gps,
                n_lev=self.data.n_lev,
                n_sar=self.data.n_sar,
                gps_names=self._file_names(self.data, 'gps'),
                lev_names=self._file_names(self.data, 'lev'),
                sar_names=self._file_names(self.data, 'sar'))
        except Exception as exc:
            logging.warning('Fit table could not be built: %s', exc)
            return []
        out = []
        for label, typ, n, wrms, wrss_n, vr in rows:
            out.append([label, typ, n, wrms, wrss_n,
                        ('%.2f%%' % (100.0 * vr)) if np.isfinite(vr) else '-'])
        return out
