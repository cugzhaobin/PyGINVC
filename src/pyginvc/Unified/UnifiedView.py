#!/usr/bin/env python
"""Plotting and result reading for unified (mixed) slip models.

``UnifiedView`` reads both the new per-fault HDF5 layout
(``/flt/schema_version == 2``) and the historical flat layout written by
``Export.Output``.  The legacy layout stores only one fault type, so it is
converted into the same in-memory group list before plotting.  This keeps
old result files readable without changing ``Export.View``.
"""

import os
import logging

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from pyginvc.Export.OutputUtils import rotate_slip
from pyginvc.Geometry.Patch import Fault

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(filename)s:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")


def _decode(value):
    """Decode an HDF5 string scalar or array into a Python string."""
    if value is None:
        return ''
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return ''
        return _decode(value.ravel()[0])
    return str(value)


def _read_optional(grp, key, default=None):
    return grp[key][()] if key in grp else default


class UnifiedView(object):
    """Read and plot a unified or legacy inversion result file."""

    def __init__(self, h5file):
        self.h5file = h5file
        self.outdir = os.path.dirname(os.path.abspath(h5file))
        self.result = self.read(h5file)

    # ------------------------------------------------------------------ #
    # Reading
    # ------------------------------------------------------------------ #
    @staticmethod
    def read(h5file):
        """Read a v2 or legacy HDF5 file into a common dictionary."""
        with h5py.File(h5file, 'r') as h5:
            obs = h5['obs'] if 'obs' in h5 else None
            flt = h5['flt'] if 'flt' in h5 else None
            sol = h5['sol'] if 'sol' in h5 else None
            if obs is None or flt is None or sol is None:
                raise ValueError(
                    '%s is not a PyGINVC solution file.' % h5file)

            result = {
                'obs': {
                    'llh_gps': _read_optional(obs, 'llh_gps',
                                              np.empty((0, 2))),
                    'llh_lev': _read_optional(obs, 'llh_lev',
                                              np.empty((0, 2))),
                    'llh_sar': _read_optional(obs, 'llh_sar',
                                              np.empty((0, 2))),
                    'unit': _read_optional(obs, 'unit', np.empty((0, 3))),
                    'd_gps': _read_optional(obs, 'd_gps', np.array([])),
                    'd_lev': _read_optional(obs, 'd_lev', np.array([])),
                    'd_sar': _read_optional(obs, 'd_sar', np.array([])),
                    'ndim': int(_read_optional(obs, 'ndim', 2)),
                    'n_sar': np.asarray(
                        _read_optional(obs, 'n_sar', []), dtype=int),
                    'n_gps': np.asarray(
                        _read_optional(obs, 'n_gps', []), dtype=int),
                    'n_lev': np.asarray(
                        _read_optional(obs, 'n_lev', []), dtype=int),
                    'd_sar_ramp': _read_optional(obs, 'd_sar_ramp', None),
                    'sar_ramp': _read_optional(obs, 'sar_ramp', None),
                },
                'sol': {
                    'slip': np.asarray(sol['slip'][()]),
                    'dhat': _read_optional(sol, 'dhat', None),
                    'r': _read_optional(sol, 'r', None),
                    'misfit': _read_optional(sol, 'misfit', None),
                    'smo_facts': _read_optional(sol, 'smo_facts', None),
                    'ruff': _read_optional(sol, 'ruff', None),
                    'moment': _read_optional(sol, 'moment', None),
                    'ramp_sar': _read_optional(sol, 'ramp_sar', None),
                    'sar_ramp_flag': _read_optional(sol, 'sar_ramp_flag', None),
                    'chosen_index': _read_optional(sol, 'chosen_index', None),
                },
                'origin': np.asarray(_read_optional(flt, 'origin',
                                                    [0.0, 0.0])),
                'faults': [],
            }

            schema = _read_optional(flt, 'schema_version', 1)
            try:
                schema = int(np.asarray(schema).ravel()[0])
            except (TypeError, ValueError, IndexError):
                schema = 1
            result['schema_version'] = schema

            if schema >= 2 and 'faults' in flt:
                base = flt['faults']
                names = sorted(base.keys(), key=lambda s: int(s))
                for name in names:
                    fg = base[name]
                    kind = _decode(_read_optional(fg, 'kind', 'rectangle'))
                    group = {
                        'name': _decode(_read_optional(fg, 'name', name)),
                        'kind': kind,
                        'nf': int(_read_optional(fg, 'nf', 0)),
                        'slip_offset': int(
                            _read_optional(fg, 'slip_offset', 0)),
                        'rake_beta': float(
                            _read_optional(fg, 'rake_beta', 0.0)),
                    }
                    if kind == 'rectangle':
                        group['dis_geom'] = np.asarray(fg['dis_geom'][()])
                        group['geom'] = np.asarray(fg['geom'][()])
                        group['nf'] = len(group['geom'])
                    else:
                        group['vertex_enu'] = np.asarray(fg['vertex_enu'][()])
                        group['vertex_llh'] = np.asarray(fg['vertex_llh'][()])
                        group['element'] = np.asarray(fg['element'][()])
                        group['nf'] = len(group['element'])
                    result['faults'].append(group)
            else:
                # Legacy flat layout: one rectangle table or one triangle mesh.
                kind = _decode(_read_optional(flt, 'flttype', 'rectangle'))
                group = {'name': 'fault01', 'kind': kind, 'slip_offset': 0,
                         'rake_beta': 0.0}
                if kind == 'rectangle':
                    group['dis_geom'] = np.asarray(flt['dis_geom'][()])
                    group['geom'] = np.asarray(flt['geom'][()])
                    group['nf'] = len(group['geom'])
                else:
                    group['vertex_enu'] = np.asarray(flt['vertex_enu'][()])
                    group['vertex_llh'] = np.asarray(flt['vertex_llh'][()])
                    group['element'] = np.asarray(flt['element'][()])
                    group['nf'] = len(group['element'])
                result['faults'].append(group)

        return result

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    def model_row(self, index=None):
        """Return one flattened slip model (default: chosen/last)."""
        slip = np.asarray(self.result['sol']['slip'], dtype=float)
        if slip.ndim == 1:
            return slip
        if index is None:
            chosen = self.result['sol'].get('chosen_index')
            if chosen is not None and np.size(chosen):
                index = int(np.asarray(chosen).reshape(-1)[-1])
            else:
                index = slip.shape[0] - 1
        return slip[int(index)]

    def _fault_slip(self, group, row):
        start = group['slip_offset']
        nf = int(group['nf'])
        block = np.asarray(row[start:start + 3 * nf]).reshape(nf, 3)
        # The solver works in the per-fault rake frame; undo the rotation so
        # plots match the slipmodel files written by UnifiedOutput.
        beta = float(group.get('rake_beta', 0.0) or 0.0)
        if beta:
            block = rotate_slip(block, beta)
        return block

    @staticmethod
    def _rectangle_corners(group, origin):
        """Return ENU polygon corners for a rectangle group."""
        flt = Fault('', group.get('nsegs', 1), int(group['nf']), False,
                    origin=origin)
        flt.nf = int(group['nf'])
        flt.dis_geom_grid = np.asarray(group['dis_geom'])
        flt.geom_grid = np.asarray(group['geom'])
        df = flt.FaultGeom2AllVertex()
        keys = [('lb_e', 'lb_n', 'lb_u'), ('rb_e', 'rb_n', 'rb_u'),
                ('rt_e', 'rt_n', 'rt_u'), ('lt_e', 'lt_n', 'lt_u')]
        return np.stack(
            [np.column_stack([np.asarray(df[a]), np.asarray(df[b]),
                              np.asarray(df[c])]) for a, b, c in keys],
            axis=1)

    @staticmethod
    def _triangle_corners(group):
        vertex = np.asarray(group['vertex_enu'])
        element = np.asarray(group['element'], dtype=int)
        return vertex[element - 1]

    def _polygons(self, group):
        if group['kind'] == 'triangle':
            corners = self._triangle_corners(group)
        else:
            corners = self._rectangle_corners(group, self.result['origin'])
        # vertex_enu uses depth positive down; plot z positive up.
        return [np.column_stack([c[:, 0], c[:, 1], -c[:, 2]])
                for c in corners]

    @staticmethod
    def _save_fig(fig, prefix, show=False):
        pdf = prefix + '.pdf'
        png = prefix + '.png'
        fig.savefig(pdf, bbox_inches='tight')
        fig.savefig(png, dpi=200, bbox_inches='tight')
        if show:
            plt.show()
        plt.close(fig)
        return pdf, png

    # ------------------------------------------------------------------ #
    # Slip distribution
    # ------------------------------------------------------------------ #
    def plot_slip_3d(self, depth_exaggeration=1.0, outfile=None, show=False,
                     index=None):
        """Plot all fault patches in one 3D view, coloured by slip."""
        row = self.model_row(index)
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        all_polys, all_vals = [], []
        for group in self.result['faults']:
            polys = self._polygons(group)
            slip = self._fault_slip(group, row)
            mag = np.linalg.norm(slip[:, :2], axis=1)
            all_polys.extend(polys)
            all_vals.extend(mag.tolist())

        all_vals = np.asarray(all_vals)
        vmin = float(np.nanmin(all_vals)) if all_vals.size else 0.0
        vmax = float(np.nanmax(all_vals)) if all_vals.size else 1.0
        norm = matplotlib.colors.Normalize(
            vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1.0)
        cmap = plt.get_cmap('jet')

        collection = Poly3DCollection(
            all_polys, facecolors=cmap(norm(all_vals)),
            edgecolors='k', linewidths=0.25, alpha=0.95)
        ax.add_collection3d(collection)

        pts = np.vstack(all_polys)
        ax.set_xlim(pts[:, 0].min(), pts[:, 0].max())
        ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
        z = pts[:, 2]
        zmid = 0.5 * (z.min() + z.max())
        span = max(z.max() - z.min(), 1e-6) / max(depth_exaggeration, 1e-6)
        ax.set_zlim(zmid - span, zmid + span)
        ax.set_xlabel('East (km)')
        ax.set_ylabel('North (km)')
        ax.set_zlabel('Up (km)')
        ax.set_title('Slip distribution (mm)')

        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.08, label='Slip (mm)')

        prefix = outfile or os.path.join(self.outdir, 'slip_3d')
        return self._save_fig(fig, prefix, show=show)

    # ------------------------------------------------------------------ #
    # Observations
    # ------------------------------------------------------------------ #
    def _split_dhat(self):
        obs = self.result['obs']
        dhat = self.result['sol']['dhat']
        if dhat is None:
            raise ValueError('The solution file has no modeled data (dhat).')
        dhat = np.asarray(dhat).reshape(-1)
        ndim = int(obs['ndim'])
        n_gps = obs['n_gps']
        if n_gps.size == 0:
            n_gps = np.array([len(obs['llh_gps'])], dtype=int)
        n_lev = obs['n_lev']
        if n_lev.size == 0:
            n_lev = np.array([len(obs['llh_lev'])], dtype=int)
        n_sar = obs['n_sar']
        if n_sar.size == 0:
            n_sar = np.array([len(obs['llh_sar'])], dtype=int)
        return dhat, ndim, n_gps, n_lev, n_sar

    def _fault_polygons_llh(self):
        """Per-patch lon/lat polygons and total slip for every fault.

        Returns ``(polygons, totals)`` where ``polygons`` is a list of
        ``(ncorner, 2)`` lon/lat arrays (4 corners for rectangles, 3 for
        triangles) and ``totals`` is the matching total slip in mm.
        """
        row = self.model_row()
        origin = self.result['origin']
        polygons, totals = [], []
        for group in self.result['faults']:
            slip = self._fault_slip(group, row)
            totals.append(np.linalg.norm(slip, axis=1))
            if group['kind'] == 'triangle':
                vertex = np.asarray(group['vertex_llh'], dtype=float)
                element = np.asarray(group['element'], dtype=int) - 1
                corners = vertex[element]          # (nf, 3, [lat, lon, dep])
                polygons.extend(np.stack(
                    [corners[:, :, 1], corners[:, :, 0]], axis=-1))
            else:
                flt = Fault.__new__(Fault)
                flt.nf = int(group['nf'])
                flt.geom_grid = np.asarray(group['geom'], dtype=float)
                dis = np.asarray(group['dis_geom'], dtype=float).copy()
                dis[:, 7:10] = slip
                flt.dis_geom_grid = dis
                flt.origin = np.asarray(origin, dtype=float)
                df = flt.FaultGeom2AllVertex()
                polygons.extend(np.stack(
                    [np.column_stack([df[tag + '_lon'], df[tag + '_lat']])
                     for tag in ('lb', 'rb', 'rt', 'lt')], axis=1))
        totals = (np.concatenate(totals) if totals
                  else np.array([], dtype=float))
        return polygons, totals

    @staticmethod
    def _polygon_extent(polygons):
        """Bounding ``[lon0, lon1, lat0, lat1]`` of a polygon list."""
        lon = np.concatenate([p[:, 0] for p in polygons])
        lat = np.concatenate([p[:, 1] for p in polygons])
        return [float(lon.min()), float(lon.max()),
                float(lat.min()), float(lat.max())]

    @staticmethod
    def _draw_outline(ax, polygons, transform):
        """Draw fault patch outlines on a cartopy GeoAxes."""
        if not polygons:
            return
        ax.add_collection(PolyCollection(
            polygons, transform=transform, facecolors='none',
            edgecolors='black', linewidths=0.5))

    @staticmethod
    def _draw_slip(ax, polygons, totals, transform):
        """Fill fault patches by total slip and return the mappable."""
        if not polygons or totals.size == 0:
            return None
        coll = PolyCollection(polygons, transform=transform,
                              cmap=plt.cm.jet, edgecolors='black',
                              linewidths=0.3)
        coll.set_array(np.asarray(totals, dtype=float))
        coll.set_clim(0.0, max(float(np.nanmax(totals)), 1e-6))
        ax.add_collection(coll)
        return coll

    def plot_obs_mod(self, scale=None, show=False, outfile=None):
        """Plot GPS and InSAR observed versus modeled data on maps.

        Uses cartopy ``PlateCarree`` panels (coastlines, borders, gridlines)
        and overlays the mixed-fault geometry: outlines on the observation /
        prediction panels and a slip-coloured fill on the residual panel.
        Falls back to plain lon/lat axes when cartopy is unavailable.
        """
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
        except ImportError:
            logging.warning('cartopy is not installed; falling back to plain '
                            'lon/lat axes for observed/modeled plots.')
            return self._plot_obs_mod_plain(scale=scale, show=show,
                                            outfile=outfile)
        from pyginvc.Export import View as EView

        obs = self.result['obs']
        dhat, ndim, n_gps, n_lev, n_sar = self._split_dhat()
        transform = ccrs.PlateCarree()
        polygons, totals = self._fault_polygons_llh()
        polygons_ext = (self._polygon_extent(polygons) if polygons
                        else None)

        gps_llh = np.asarray(obs['llh_gps'])
        gps_obs = np.asarray(obs['d_gps']).reshape(-1)
        nsta = len(gps_llh)
        has_gps = nsta > 0 and (int(np.sum(n_gps)) if n_gps.size else 0) > 0
        has_sar = np.asarray(obs['llh_sar']).size > 0

        if not has_gps and not has_sar:
            logging.warning('No GPS or InSAR data to plot.')
            return []

        written = []
        offset = 0

        if has_gps:
            # All GPS files share one overview map.
            obs2 = gps_obs[:nsta * ndim].reshape(nsta, ndim)
            mod2 = dhat[offset:offset + nsta * ndim].reshape(nsta, ndim)
            offset += nsta * ndim
            lon, lat = gps_llh[:, 1], gps_llh[:, 0]
            ext = EView._data_extent(lon, lat)
            if polygons_ext is not None:
                ext = EView._union_extent(ext, polygons_ext)
            qscale = scale
            if qscale is None:
                qscale = EView._auto_quiver_scale(
                    lon, lat, (obs2[:, 0], obs2[:, 1]),
                    (mod2[:, 0], mod2[:, 1]))

            figsize, rects, _, vcbar = EView._geo_layout(
                ext, ncols=2, vcbar=True)
            fig = plt.figure(figsize=figsize)
            ax1 = fig.add_axes(rects[0], projection=ccrs.PlateCarree())
            ax2 = fig.add_axes(rects[1], projection=ccrs.PlateCarree())
            cax = fig.add_axes(vcbar)
            EView._setup_map(ax1, ext, ccrs, cfeature, left_labels=True)
            EView._setup_map(ax2, ext, ccrs, cfeature, left_labels=False)
            self._draw_outline(ax1, polygons, transform)
            slip_coll = self._draw_slip(ax2, polygons, totals, transform)

            ax1.quiver(lon, lat, obs2[:, 0], obs2[:, 1], transform=transform,
                       color='b', scale=qscale, scale_units='xy')
            ax1.quiver(lon, lat, mod2[:, 0], mod2[:, 1], transform=transform,
                       color='r', scale=qscale, scale_units='xy')
            if ndim == 3:
                ax1.scatter(lon, lat, c=obs2[:, 2], cmap=plt.cm.jet, s=50,
                            linewidths=0.1, edgecolors='black',
                            transform=transform)
                ax1.scatter(lon, lat, c=mod2[:, 2], cmap=plt.cm.jet, s=15,
                            linewidths=0.1, edgecolors='black',
                            transform=transform)
            ax1.set_title('GPS observation (blue) vs prediction (red)')

            res = obs2 - mod2
            ax2.quiver(lon, lat, res[:, 0], res[:, 1], transform=transform,
                       color='b', scale=qscale, scale_units='xy')
            if ndim == 3:
                ax2.scatter(lon, lat, c=res[:, 2], cmap=plt.cm.jet,
                            linewidths=0.1, edgecolors='black',
                            transform=transform)
            ax2.set_title('GPS residual')
            if slip_coll is not None:
                cb = fig.colorbar(slip_coll, cax=cax)
                cb.set_label('Total slip (mm)')

            if outfile:
                prefix = outfile
            else:
                prefix = os.path.join(self.outdir, 'gps_mod_obs')
            written.extend(self._save_fig(fig, prefix, show=show))

        # Leveling rows sit between GPS and InSAR in the data vector.
        offset += int(np.sum(n_lev))

        if has_sar:
            sar_llh = np.asarray(obs['llh_sar'])
            d_obs = np.asarray(obs['d_sar']).reshape(-1)
            n_model = min(len(sar_llh), len(dhat) - offset)
            d_obs = d_obs[:n_model]
            d_mod = dhat[offset:offset + n_model]

            counts = np.asarray(n_sar, dtype=int)
            if counts.size == 0 or counts.sum() != len(sar_llh):
                counts = np.array([len(sar_llh)], dtype=int)
            multi = len(counts) > 1

            start = 0
            for itrack, count in enumerate(counts):
                count = int(count)
                end = start + count
                if count <= 0:
                    start = end
                    continue
                ob = d_obs[start:end]
                mb = d_mod[start:end]
                rb = ob - mb
                lon = sar_llh[start:end, 1]
                lat = sar_llh[start:end, 0]
                tag = ' Track %d' % (itrack + 1) if multi else ''

                ext = EView._data_extent(lon, lat)
                if polygons_ext is not None:
                    ext = EView._union_extent(ext, polygons_ext)
                norm_fit = EView._div_norm(np.concatenate((ob, mb)))
                norm_res = EView._div_norm(rb)

                figsize, rects, hcbars, vcbar = EView._geo_layout(
                    ext, ncols=3, vcbar=True, hcbar=True)
                fig = plt.figure(figsize=figsize)
                axes = [fig.add_axes(rects[k],
                                     projection=ccrs.PlateCarree())
                        for k in range(3)]
                cbars = [fig.add_axes(hcbars[k]) for k in range(3)]
                cax_slip = fig.add_axes(vcbar)
                EView._setup_map(axes[0], ext, ccrs, cfeature,
                                 left_labels=True)
                EView._setup_map(axes[1], ext, ccrs, cfeature,
                                 left_labels=False)
                EView._setup_map(axes[2], ext, ccrs, cfeature,
                                 left_labels=False)
                self._draw_outline(axes[0], polygons, transform)
                self._draw_outline(axes[1], polygons, transform)
                slip_coll = self._draw_slip(axes[2], polygons, totals,
                                            transform)

                for k, (values, norm, title, label) in enumerate((
                        (ob, norm_fit,
                         'InSAR%s LOS Observation' % tag,
                         'Line of sight changes (mm)'),
                        (mb, norm_fit,
                         'InSAR%s LOS Prediction' % tag,
                         'Line of sight changes (mm)'),
                        (rb, norm_res,
                         'InSAR%s LOS Residual' % tag,
                         'Line of sight residual (mm)'))):
                    ax = axes[k]
                    sc = ax.scatter(lon, lat, c=values, cmap=plt.cm.RdBu_r,
                                    norm=norm, s=2, transform=transform)
                    ax.set_title(title)
                    cb = fig.colorbar(sc, cax=cbars[k],
                                      orientation='horizontal')
                    cb.set_label(label)
                if slip_coll is not None:
                    cb = fig.colorbar(slip_coll, cax=cax_slip)
                    cb.set_label('Total slip (mm)')

                fname = ('sar_mod_obs_%d' % (itrack + 1) if multi
                         else 'sar_mod_obs')
                if outfile:
                    # Do not overwrite the GPS figure when both are asked
                    # for through one explicit prefix.
                    prefix = outfile if not has_gps else outfile + '_' + fname
                else:
                    prefix = os.path.join(self.outdir, fname)
                written.extend(self._save_fig(fig, prefix, show=show))
                start = end

        return written

    def _plot_obs_mod_plain(self, scale=None, show=False, outfile=None):
        """Plain lon/lat fallback used when cartopy is not installed."""
        obs = self.result['obs']
        dhat, ndim, n_gps, n_lev, n_sar = self._split_dhat()
        offset = 0

        gps_llh = np.asarray(obs['llh_gps'])
        gps_obs = np.asarray(obs['d_gps']).reshape(-1)
        n_gps_total = int(np.sum(n_gps)) if n_gps.size else 0
        has_gps = gps_llh.size > 0 and n_gps_total > 0
        has_sar = np.asarray(obs['llh_sar']).size > 0

        nrows = int(has_gps) + int(has_sar)
        if nrows == 0:
            logging.warning('No GPS or InSAR data to plot.')
            return []
        fig, axes = plt.subplots(nrows, 2 if not has_sar else 3,
                                 figsize=(12, 4.0 * nrows), squeeze=False)

        r = 0
        if has_gps:
            nsta_total = len(gps_llh)
            obs2 = gps_obs[:nsta_total * ndim].reshape(nsta_total, ndim)
            mod2 = dhat[offset:offset + nsta_total * ndim].reshape(
                nsta_total, ndim)
            offset += nsta_total * ndim
            lon, lat = gps_llh[:, 1], gps_llh[:, 0]
            for j, (title, data) in enumerate(
                    (('Observed GPS', obs2), ('Modeled GPS', mod2))):
                ax = axes[r, j]
                ax.quiver(lon, lat, data[:, 0], data[:, 1], color='C%d' % j,
                          angles='xy', scale_units='xy', scale=scale)
                ax.set_title(title)
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')
                ax.set_aspect('equal', adjustable='datalim')
            ax = axes[r, 2] if axes.shape[1] > 2 else None
            if ax is not None:
                if ndim == 3:
                    ax.scatter(lon, obs2[:, 2], c='k', s=10, label='obs')
                    ax.scatter(lon, mod2[:, 2], c='r', s=10, marker='x',
                               label='model')
                    ax.legend()
                ax.set_title('GPS vertical')
            r += 1

        offset += int(np.sum(n_lev))

        if has_sar:
            sar_llh = np.asarray(obs['llh_sar'])
            d_obs = np.asarray(obs['d_sar']).reshape(-1)
            n_sar_total = len(sar_llh)
            n_model = min(n_sar_total, len(dhat) - offset)
            d_mod = dhat[offset:offset + n_model]
            d_obs = d_obs[:n_model]
            lon, lat = sar_llh[:n_model, 1], sar_llh[:n_model, 0]
            for j, (title, value) in enumerate(
                    (('Observed InSAR', d_obs),
                     ('Modeled InSAR', d_mod),
                     ('Residual InSAR', d_obs - d_mod))):
                ax = axes[r, j]
                sc = ax.scatter(lon, lat, c=value, s=4, cmap='jet')
                fig.colorbar(sc, ax=ax, shrink=0.7, label='LOS (mm)')
                ax.set_title(title)
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')
                ax.set_aspect('equal', adjustable='datalim')

        fig.tight_layout()
        prefix = outfile or os.path.join(self.outdir, 'obs_mod')
        return list(self._save_fig(fig, prefix, show=show))

    def plot_sar_ramp(self, show=False, outfile=None):
        """Plot the estimated per-track InSAR ramp contribution."""
        obs = self.result['obs']
        ramp = obs.get('d_sar_ramp')
        flags = obs.get('sar_ramp')
        llh = np.asarray(obs['llh_sar'])
        if ramp is None or llh.size == 0:
            logging.info('No InSAR ramp to plot.')
            return None
        ramp = np.asarray(ramp).reshape(-1)
        n_sar = obs['n_sar']
        if n_sar.size == 0:
            n_sar = np.array([len(llh)])
        flags = ([True] * len(n_sar) if flags is None
                 else np.asarray(flags, dtype=bool))
        starts = np.concatenate([[0], np.cumsum(n_sar)]).astype(int)
        tracks = [(i, starts[i], starts[i + 1]) for i, n in enumerate(n_sar)
                  if n > 0 and i < len(flags) and flags[i]]
        if not tracks:
            logging.info('No track has an estimated tilt ramp.')
            return None

        fig, axes = plt.subplots(1, len(tracks),
                                 figsize=(5 * len(tracks), 4), squeeze=False)
        for ax, (i, start, end) in zip(axes[0], tracks):
            sc = ax.scatter(llh[start:end, 1], llh[start:end, 0],
                            c=ramp[start:end], s=5, cmap='jet')
            fig.colorbar(sc, ax=ax, shrink=0.7, label='Ramp (mm)')
            ax.set_title('InSAR track %d' % (i + 1))
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.set_aspect('equal', adjustable='datalim')
        fig.tight_layout()
        prefix = outfile or os.path.join(self.outdir, 'sar_ramp')
        return self._save_fig(fig, prefix, show=show)

    # ------------------------------------------------------------------ #
    # Trade-off curve
    # ------------------------------------------------------------------ #
    def plot_lcurve(self, show=False, outfile=None):
        """Plot roughness versus misfit for the smoothing scan."""
        sol = self.result['sol']
        ruff = sol.get('ruff')
        misfit = sol.get('misfit')
        smo = sol.get('smo_facts')
        if ruff is None or misfit is None:
            logging.info('No roughness/misfit history to plot.')
            return None
        ruff = np.asarray(ruff).reshape(-1)
        misfit = np.asarray(misfit).reshape(-1)
        if ruff.size != misfit.size:
            raise ValueError('ruff and misfit have different lengths.')

        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(misfit, ruff, '-o')
        if smo is not None and np.size(smo) == ruff.size:
            for x, y, s in zip(misfit, ruff, np.asarray(smo).reshape(-1)):
                ax.annotate('%.3g' % s, (x, y), fontsize=8,
                            textcoords='offset points', xytext=(4, 4))
        ax.set_xlabel('Misfit (WRSS)')
        ax.set_ylabel('Roughness')
        ax.set_title('Smoothing trade-off curve')
        ax.grid(True, alpha=0.3)
        if (misfit > 0).all() and (ruff > 0).all():
            ax.set_xscale('log')
            ax.set_yscale('log')
        fig.tight_layout()
        prefix = outfile or os.path.join(self.outdir, 'lcurve')
        return self._save_fig(fig, prefix, show=show)
