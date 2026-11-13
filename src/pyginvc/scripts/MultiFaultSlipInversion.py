#!/usr/bin/env python
"""
MultiFaultSlipInversion.py -- 多断层 + 多数据文件滑动反演的命令行入口。

反演类实现位于 pyginvc.Inversion.MultiFaultInversion，这里只负责：
    * 读取 YAML 配置；
    * 逐断层构建 Fault，逐文件加载数据并聚合；
    * 调用 MultiFaultInversion 完成反演并打印汇总。

用法：
    python MultiFaultSlipInversion.py --cfgfile config/slipinv_multi.yaml

限制：
    * 所有 GPS 文件必须是同一种格式（2D 或 3D），不混用；
    * 水准数据尚未实现：配置 leveling 会明确报错，而不是静默忽略。
"""

import os
import sys
import logging
import argparse
import numpy as np
import yaml

from pyginvc.Geometry.Patch import Fault
from pyginvc.Geometry.FaultCollection import FaultCollection
from pyginvc.GeoData.GPSData import GPSData
from pyginvc.GeoData.SARData import SARData
from pyginvc.GeoData.LEVData import LEVData
from pyginvc.GeoData.GeoDataCollection import GeoDataCollection
from pyginvc.Inversion.MultiFaultInversion import (
    MultiFaultInversion,
)
from pyginvc.Export.MultiOutput import MultiOutput
from pyginvc.Export import View

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

__all__ = [
    'MultiFaultInversion',
    'MultiOutput',
    'build_origin',
    'build_faults',
    'build_data',
    'resolve_gfiletype',
    'solution_file',
    'run',
    'main',
]


def _has_file(entry):
    """Return True when a config entry names an input file.

    None (YAML null / empty ``file:``) and blank strings are placeholders
    and are skipped, so a config can keep unused example entries.
    """
    if not isinstance(entry, dict):
        return False
    return bool(str(entry.get('file') or '').strip())


def resolve_gfiletype(global_cfg, gps_cfg):
    """Resolve the default GPS file type used by GPSData.

    Priority: ``global.gfiletype``; else the ``type`` of the first GPS entry
    that actually has a file; else 'GMT2D'.  YAML null / empty strings are
    treated as "not set", so a placeholder entry can never leak ``None``.
    """
    gt = str((global_cfg or {}).get('gfiletype') or '').strip()
    if gt:
        return gt
    for g in gps_cfg or []:
        if _has_file(g) and g.get('type'):
            gt = str(g['type']).strip()
            if gt:
                return gt
    return 'GMT2D'


def build_origin(fault_cfg):
    """
    get origin [lat, lon] from fault files

    Placeholder entries (without a file) are ignored. A clear error is raised
    when no usable fault file exists or when endpoint coordinates are not
    finite.
    """
    geoms = []
    for f in fault_cfg:
        if not _has_file(f):
            continue
        geom = np.genfromtxt(f['file'], comments='#')
        if geom.size == 0:
            raise ValueError('Fault geometry file is empty: %s' % f['file'])
        if geom.ndim == 1:
            geom = geom.reshape(1, -1)
        if geom.shape[1] < 7:
            raise ValueError(
                'Fault geometry file must have at least 7 columns, got %d: %s'
                % (geom.shape[1], f['file']))
        geoms.append(geom)

    if not geoms:
        raise ValueError('No fault file is configured; cannot derive origin.')

    allg = np.vstack(geoms)
    lat  = np.concatenate([allg[:, 3], allg[:, 5]])
    lon  = np.concatenate([allg[:, 4], allg[:, 6]])
    if not (np.isfinite(lat).all() and np.isfinite(lon).all()):
        raise ValueError(
            'Fault geometry contains non-finite endpoint coordinates.')
    return [float(np.mean(lat)), float(np.mean(lon))]


def build_faults(fault_cfg, origin):
    """
    return instance of class FaultCollection

    Placeholder entries (without a file) are skipped; configurations are
    filtered in lockstep so FaultCollection.pairs() stays aligned.
    """
    flts, configs = [], []
    for i, f in enumerate(fault_cfg):
        if not _has_file(f):
            continue
        missing = [k for k in ('nsegs', 'ndeps', 'doSubFault')
                   if k not in f]
        if missing:
            raise ValueError(
                'faults[%d] (%s) is missing required key(s): %s.'
                % (i, f.get('file'), ', '.join(missing)))
        flt = Fault(f['file'], f['nsegs'], f['ndeps'], f['doSubFault'],
                    origin=origin)
        flt.load_fault()
        flts.append(flt)
        configs.append(f)
    if not flts:
        raise ValueError('No fault file is configured.')
    return FaultCollection(flts, configs)


def build_data(gps_cfg, sar_cfg, lev_cfg, gfiletype):
    """
    return instance of class GeoDataCollection

    Entries whose 'file' is empty are placeholders and are skipped, so a
    config can use only one data type (e.g. InSAR only). Only when no data
    file at all is configured does the program abort.
    """
    files, configs = [], []

    for g in gps_cfg:
        if not _has_file(g):
            continue
        obj = GPSData(g['file'], g.get('type') or gfiletype)
        obj.load()
        files.append(obj)
        configs.append(g)

    for s in sar_cfg:
        if not _has_file(s):
            continue
        obj = SARData(s['file'],
                      scale=s.get('scale', 1.0),
                      sigma=s.get('sigma'),
                      lamda=s.get('lamda'))
        obj.load()
        files.append(obj)
        configs.append(s)

    for l in lev_cfg:
        if not _has_file(l):
            continue
        obj = LEVData(l['file'])
        obj.load()
        files.append(obj)
        configs.append(l)

    if not files:
        logging.error('No GPS/InSAR/leveling data.')
        sys.exit(1)

    return GeoDataCollection(files, configs)


def solution_file(output, outpath, filename='solutions.h5'):
    """
    return path of the solution HDF5 file

    MultiOutput archives *.h5 into archdir when it differs from outdir, so
    prefer the archived copy and fall back to the output directory.
    """
    archdir = getattr(output, 'archdir', None)
    if archdir:
        path = os.path.join(archdir, filename)
        if os.path.isfile(path):
            return path
    return os.path.join(outpath, filename)


# --------------------------------------------------------------------------- #
# Main program
# --------------------------------------------------------------------------- #
def run(cfgfile):
    with open(cfgfile, 'r') as fid:
        cfg = yaml.safe_load(fid) or {}
    if not isinstance(cfg, dict):
        raise ValueError(
            'Config file %s must contain a YAML mapping.' % cfgfile)

    # YAML null sections (section present but empty) must behave like an
    # omitted section; otherwise None.get / iteration over None crashes.
    G      = cfg.get('global') or {}
    FAULTS = cfg.get('faults') or []
    GPSC   = cfg.get('gps') or []
    SARC   = cfg.get('insar') or []
    LEVC   = cfg.get('leveling') or []
    EXPC   = cfg.get('export') or {}

    for name, value in (('faults', FAULTS), ('gps', GPSC),
                        ('insar', SARC), ('leveling', LEVC)):
        if not isinstance(value, list):
            raise ValueError(
                'Config section %r must be a list, got %s.'
                % (name, type(value).__name__))

    if not FAULTS:
        logging.error('No faults defined in config.')
        sys.exit(1)

    # 全局 GPS 格式：优先 global.gfiletype，否则取第一条有数据的 GPS 的 type
    gfiletype = resolve_gfiletype(G, GPSC)

    # 1) get origin from YAML file or from fault files
    origin = G.get('origin')
    if origin:
        if not isinstance(origin, (list, tuple)) or len(origin) != 2:
            raise ValueError(
                'global.origin must be [lat, lon], got %r' % (origin,))
        try:
            origin = [float(origin[0]), float(origin[1])]
        except (TypeError, ValueError):
            raise ValueError(
                'global.origin must contain two numbers, got %r' % (origin,))
        if not np.isfinite(origin).all():
            raise ValueError(
                'global.origin must be finite, got %r' % (origin,))
    else:
        origin = build_origin(FAULTS)
    faults = build_faults(FAULTS, origin)

    # 2) get geodetic data
    data = build_data(GPSC, SARC, LEVC, gfiletype)

    # 3) 反演（格林函数、平滑、ramp、权重、边界与扫描求解）
    inversion = MultiFaultInversion(data, faults, origin=origin,
                                    config=G, export=EXPC)

    # 输出全部交给 MultiOutput：逐文件观测拟合、逐断层滑移模型、
    # HDF5/summary/L 曲线/GMT 滑移与归档
    outpath = EXPC.get('outpath', './')
    output  = MultiOutput(data, faults, sol=inversion,
                          archdir=EXPC.get('archdir', outpath),
                          outdir=outpath)
    results = inversion.run(output=output)

    # 4) 汇总
    print('\n=== Summary (sm, WRSS_gps, WRSS_sar, roughness, Mo, Mw) ===')
    for row in results:
        print('%7.3f  %12.3f  %12.3f  %12.3f  %12.3e  %6.2f' % row)
    print('Slip models written to %s/' % outpath)

    # 5) 可视化：观测值与模拟值对比、滑动分布、折衷（L）曲线
    if EXPC.get('view', False):
        h5file = solution_file(output, outpath)
        View.plot_lcurve(h5file)
        View.plot_sar_ramp(h5file)
        View.plot_obs_mod(h5file, scale=EXPC.get('vecscale'))
        View.plot_slip_3d(h5file,
                          depth_exaggeration=EXPC.get('depth_exaggeration', 2.0))


def main():
    parser = argparse.ArgumentParser(description='Multi-fault / multi-file slip inversion.')
    parser.add_argument('--cfgfile', type=str, required=True, help='path to YAML config')
    args = parser.parse_args()
    run(args.cfgfile)


if __name__ == '__main__':
    main()
