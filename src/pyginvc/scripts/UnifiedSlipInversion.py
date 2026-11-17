#!/usr/bin/env python
"""Command-line entry point for unified rectangular/triangular inversion.

The runner reads a ``schema: unified-v1`` YAML file, loads every enabled
fault (rectangle or triangle), loads the requested GPS/InSAR/leveling files,
and hands everything to
:class:`pyginvc.Inversion.UnifiedSlipInversion.UnifiedSlipInversion`.

Usage::

    pyginvc unified-slip --cfgfile config/unified_cfg.yaml
    pyginvc unified-slip --cfgfile config/unified_cfg.yaml --outpath ./run01

All enabled faults share one ``[lat, lon]`` origin.  Set ``global.origin``
explicitly or let the runner derive it from the enabled fault files.
"""

import os
import logging
import argparse

import numpy as np
import yaml

from pyginvc.GeoData.GPSData import GPSData
from pyginvc.GeoData.SARData import SARData
from pyginvc.GeoData.LEVData import LEVData
from pyginvc.GeoData.GeoDataCollection import GeoDataCollection
from pyginvc.Geometry.Patch import Fault
from pyginvc.Geometry.Triangle import Triangle
from pyginvc.Geometry.UnifiedFaultCollection import UnifiedFaultCollection
from pyginvc.Inversion.UnifiedSlipInversion import UnifiedSlipInversion
from pyginvc.Unified.UnifiedOutput import UnifiedOutput

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(filename)s:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")

__all__ = [
    'build_origin', 'build_faults', 'build_data', 'resolve_gfiletype',
    'solution_file', '_resolve_view', 'run', 'main',
]

SCHEMA = 'unified-v1'


def _as_bool(value):
    """Coerce a YAML scalar to bool."""
    if isinstance(value, str):
        return value.strip().lower() in ('1', 'true', 'yes', 'on')
    return bool(value)


def _has_file(entry, key='file'):
    if not isinstance(entry, dict):
        return False
    return bool(str(entry.get(key) or '').strip())


def _is_enabled(entry, key='file'):
    """Return True when an entry is enabled and names an input file.

    Triangle faults identify their input with ``vertexfile`` instead of
    ``file``; when the generic ``file`` key is probed, accept either.
    """
    if not isinstance(entry, dict):
        return False
    enabled = entry.get('enable', True)
    if enabled is not None and not _as_bool(enabled):
        return False
    if _has_file(entry, key):
        return True
    return key == 'file' and _has_file(entry, 'vertexfile')


def _fault_type(entry, index):
    kind = str(entry.get('type') or '').strip().lower()
    if not kind:
        kind = 'triangle' if _has_file(entry, 'vertexfile') else 'rectangle'
    if kind not in ('rectangle', 'triangle'):
        raise ValueError(
            'faults[%d] has unsupported type %r; use rectangle or triangle.'
            % (index, kind))
    return kind


def _rectangle_origin(entry):
    geom = np.genfromtxt(entry['file'], comments='#')
    if geom.size == 0:
        raise ValueError('Fault geometry file is empty: %s' % entry['file'])
    if geom.ndim == 1:
        geom = geom.reshape(1, -1)
    if geom.shape[1] < 7:
        raise ValueError(
            'Rectangle fault file must have at least 7 columns, got %d: %s'
            % (geom.shape[1], entry['file']))
    lat = np.concatenate([geom[:, 3], geom[:, 5]])
    lon = np.concatenate([geom[:, 4], geom[:, 6]])
    return lat, lon


def _triangle_origin(entry):
    vertex = np.genfromtxt(entry['vertexfile'], comments='#')
    if vertex.size == 0:
        raise ValueError(
            'Triangle vertex file is empty: %s' % entry['vertexfile'])
    if vertex.ndim == 1:
        vertex = vertex.reshape(1, -1)
    if vertex.shape[1] < 3:
        raise ValueError(
            'Triangle vertex file needs at least 3 columns (index lat lon): %s'
            % entry['vertexfile'])
    return vertex[:, 1], vertex[:, 2]


def build_origin(fault_cfg):
    """Derive one ``[lat, lon]`` origin from every enabled fault."""
    lats, lons = [], []
    for i, entry in enumerate(fault_cfg):
        if not _is_enabled(entry):
            continue
        kind = _fault_type(entry, i)
        if kind == 'rectangle':
            lat, lon = _rectangle_origin(entry)
        else:
            lat, lon = _triangle_origin(entry)
        lats.append(np.asarray(lat, dtype=float).ravel())
        lons.append(np.asarray(lon, dtype=float).ravel())
    if not lats:
        raise ValueError('No enabled fault file is configured.')
    lat = np.concatenate(lats)
    lon = np.concatenate(lons)
    if not (np.isfinite(lat).all() and np.isfinite(lon).all()):
        raise ValueError('Fault files contain non-finite coordinates.')
    return [float(np.mean(lat)), float(np.mean(lon))]


def build_faults(fault_cfg, origin):
    """Load enabled faults and return a ``UnifiedFaultCollection``."""
    flts, configs = [], []
    for i, entry in enumerate(fault_cfg):
        if not _is_enabled(entry):
            continue
        kind = _fault_type(entry, i)
        if kind == 'rectangle':
            required = ('file', 'nsegs', 'ndeps')
            missing = [k for k in required if k not in entry]
            if missing:
                raise ValueError(
                    'Rectangle faults[%d] is missing key(s): %s.'
                    % (i, ', '.join(missing)))
            flt = Fault(entry['file'], entry['nsegs'], entry['ndeps'],
                        entry.get('doSubFault', True), origin=origin)
        else:
            required = ('vertexfile', 'elementfile')
            missing = [k for k in required if k not in entry]
            if missing:
                raise ValueError(
                    'Triangle faults[%d] is missing key(s): %s.'
                    % (i, ', '.join(missing)))
            flt = Triangle(entry['vertexfile'], entry['elementfile'],
                           origin=origin)
        flt.load_fault()
        flts.append(flt)
        configs.append(entry)

    if not flts:
        raise ValueError('No enabled fault is configured.')
    return UnifiedFaultCollection(flts, configs)


def resolve_gfiletype(global_cfg, gps_cfg):
    """Resolve the GPS file type from global or the first enabled entry."""
    value = str((global_cfg or {}).get('gfiletype') or '').strip()
    if value:
        return value
    for entry in gps_cfg or []:
        if _is_enabled(entry) and entry.get('type'):
            value = str(entry['type']).strip()
            if value:
                return value
    return 'GMT2D'


def build_data(gps_cfg, sar_cfg, lev_cfg, gfiletype):
    """Load every enabled observation file and aggregate it."""
    files, configs = [], []
    for entry in gps_cfg:
        if not _is_enabled(entry):
            continue
        obj = GPSData(entry['file'], entry.get('type') or gfiletype)
        obj.load()
        files.append(obj)
        configs.append(entry)

    for entry in sar_cfg:
        if not _is_enabled(entry):
            continue
        obj = SARData(entry['file'],
                      scale=entry.get('scale', 1.0),
                      sigma=entry.get('sigma'),
                      lamda=entry.get('lamda'))
        obj.load()
        files.append(obj)
        configs.append(entry)

    for entry in lev_cfg:
        if not _is_enabled(entry):
            continue
        obj = LEVData(entry['file'])
        try:
            obj.load()
        except NotImplementedError as exc:
            raise NotImplementedError(
                'Leveling support is not implemented in the unified runner '
                '(configured file: %s).' % entry['file']) from exc
        files.append(obj)
        configs.append(entry)

    if not files:
        raise ValueError('No GPS/InSAR/leveling data are enabled.')
    return GeoDataCollection(files, configs)


def _check_compatibility(faults, data):
    """Reject combinations that the available Green's functions cannot model."""
    kinds = set(faults.kinds())
    has_lev = len(np.asarray(data.llh_lev)) > 0
    if has_lev and 'triangle' in kinds:
        raise ValueError(
            'Leveling data cannot be combined with triangular faults: '
            'the triangle Green classes do not generate leveling rows.')
    if has_lev:
        raise NotImplementedError(
            'Leveling data are not supported by the unified runner.')


def solution_file(output, outpath, filename='solutions.h5'):
    """Return the working copy of the solution file."""
    return os.path.join(outpath, filename)


_NONINTERACTIVE_BK = ('agg', 'pdf', 'ps', 'svg', 'cairo', 'template',
                      'pgf', 'inline')


def _resolve_view(export_cfg):
    """Return True when figures should be shown interactively."""
    value = (export_cfg or {}).get('view')
    if value is None:
        value = (export_cfg or {}).get('show')
    show = _as_bool(value) if value is not None else False
    if show:
        try:
            import matplotlib
            backend = str(matplotlib.get_backend()).lower()
        except Exception:
            return show
        if any(backend.startswith(k) for k in _NONINTERACTIVE_BK):
            logging.warning(
                'export.view=True but the matplotlib backend is %r; figures '
                'are saved only.', backend)
    return show


def run(cfgfile, outpath=None):
    """Read the YAML file, run the inversion and write the viewer products."""
    with open(cfgfile, 'r') as fid:
        cfg = yaml.safe_load(fid) or {}
    if not isinstance(cfg, dict):
        raise ValueError(
            'Config file %s must contain a YAML mapping.' % cfgfile)
    if cfg.get('schema') != SCHEMA:
        raise ValueError(
            'Config file %s must declare ``schema: %s``; got %r.  Old '
            'multi-fault configs are not supported.'
            % (cfgfile, SCHEMA, cfg.get('schema')))

    G = cfg.get('global') or {}
    FAULTS = cfg.get('faults') or []
    GPSC = cfg.get('gps') or []
    SARC = cfg.get('insar') or []
    LEVC = cfg.get('leveling') or []
    EXPC = dict(cfg.get('export') or {})
    if outpath is not None:
        EXPC['outpath'] = outpath

    for name, value in (('faults', FAULTS), ('gps', GPSC),
                        ('insar', SARC), ('leveling', LEVC)):
        if not isinstance(value, list):
            raise ValueError(
                'Config section %r must be a list, got %s.'
                % (name, type(value).__name__))
    if not FAULTS:
        raise ValueError('No fault is defined in the config.')

    enabled = [f for f in FAULTS if _is_enabled(f)]
    if not enabled:
        raise ValueError('Every fault has enable: false or no input file.')

    origin = G.get('origin')
    if origin:
        origin = np.asarray(origin, dtype=float).reshape(-1)
        if origin.size != 2 or not np.isfinite(origin).all():
            raise ValueError(
                'global.origin must be two finite numbers [lat, lon], got %r.'
                % (origin.tolist(),))
        origin = origin.tolist()
    else:
        origin = build_origin(FAULTS)

    faults = build_faults(FAULTS, origin)
    gfiletype = resolve_gfiletype(G, GPSC)
    data = build_data(GPSC, SARC, LEVC, gfiletype)
    _check_compatibility(faults, data)

    inversion = UnifiedSlipInversion(data, faults, origin=origin,
                                     config=G, export=EXPC)

    outdir = EXPC.get('outpath') or './'
    output = UnifiedOutput(data, faults, sol=inversion,
                           outdir=outdir, config=G, export=EXPC)
    results = inversion.run(output=output)

    print('\n=== Summary (sm, WRSS_geod, WRSS_sar, roughness, Mo, Mw) ===')
    for row in results:
        print('%7.3f  %12.3f  %12.3f  %12.3f  %12.3e  %6.2f' % row)
    print('Slip models written to %s/' % outdir)

    show = _resolve_view(EXPC)
    h5file = solution_file(output, outdir)
    if _as_bool(EXPC.get('view_plots', True)):
        try:
            from pyginvc.Unified.UnifiedView import UnifiedView
            viewer = UnifiedView(h5file)
            viewer.plot_lcurve(show=show)
            viewer.plot_sar_ramp(show=show)
            viewer.plot_obs_mod(scale=EXPC.get('vecscale'), show=show)
            viewer.plot_slip_3d(
                depth_exaggeration=EXPC.get('depth_exaggeration', 2.0),
                show=show)
        except Exception as exc:
            logging.warning('UnifiedView plotting failed: %s', exc)

    embed = _as_bool(EXPC.get('report_embed', False))
    output.WriteReport(embed=embed, cfgfile=cfgfile)

    if _as_bool(EXPC.get('interactive_report', True)):
        try:
            from pyginvc.Export.PlotlyReport import make_report
            make_report(h5file,
                        depth_exaggeration=EXPC.get('depth_exaggeration',
                                                    1.0),
                        arrow_scale=EXPC.get('arrow_scale', 1.0),
                        embed=embed)
        except Exception as exc:
            logging.warning('Interactive report could not be written: %s',
                            exc)
    return inversion


def main():
    parser = argparse.ArgumentParser(
        description='Unified rectangular/triangular slip inversion.')
    parser.add_argument('--cfgfile', type=str, required=True,
                        help='path to the unified YAML config')
    parser.add_argument('--outpath', type=str, default=None,
                        help='output directory; overrides export.outpath')
    args = parser.parse_args()
    run(args.cfgfile, outpath=args.outpath)


if __name__ == '__main__':
    main()
