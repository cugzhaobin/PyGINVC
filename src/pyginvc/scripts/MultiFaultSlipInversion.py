#!/usr/bin/env python
"""
MultiFaultSlipInversion.py -- 多断层 + 多数据文件滑动反演的命令行入口。

反演类实现位于 pyginvc.Inversion.MultiFaultSlipInversion，这里只负责：
    * 读取 YAML 配置；
    * 逐断层构建 Fault，逐文件加载数据并聚合；
    * 调用 MultiFaultSlipInversion 完成反演并打印汇总。

用法：
    python MultiFaultSlipInversion.py --cfgfile config/slipinv_multi.yaml

限制：
    * 所有 GPS 文件必须是同一种格式（2D 或 3D），不混用；
    * 水准数据当前未真正加载（LEVData 为 stub），仅保留接口。
"""

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
from pyginvc.Inversion.MultiFaultSlipInversion import (
    MultiFaultSlipInversion,
    build_green,
    build_lap,
    build_ramp,
    build_bounds,
    total_moment,
    write_data_fit,
)
from pyginvc.Export.MultiOutput import MultiOutput

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

__all__ = [
    'MultiFaultSlipInversion',
    'MultiOutput',
    'build_origin',
    'build_faults',
    'build_data',
    'build_green',
    'build_lap',
    'build_ramp',
    'build_bounds',
    'total_moment',
    'write_data_fit',
    'run',
    'main',
]


# --------------------------------------------------------------------------- #
# 输入构建
# --------------------------------------------------------------------------- #
def build_origin(fault_cfg):
    """由所有断层几何算公共参考原点 [lat, lon]。"""
    geoms = [np.genfromtxt(f['file'], comments='#') for f in fault_cfg]
    geoms = [g.reshape(1, -1) if g.ndim == 1 else g for g in geoms]
    allg  = np.vstack(geoms)
    return [float(np.mean(allg[:, [3, 5]])), float(np.mean(allg[:, [4, 6]]))]


def build_faults(fault_cfg, origin):
    """每条断层建 Fault 对象并加载，返回聚合后的 FaultCollection。"""
    flts = []
    for f in fault_cfg:
        flt = Fault(f['file'], f['nsegs'], f['ndeps'], f['doSubFault'], origin=origin)
        flt.load_fault()
        flts.append(flt)
    return FaultCollection(flts, fault_cfg)


def build_data(gps_cfg, sar_cfg, lev_cfg, gfiletype):
    """每个文件建一个数据对象，再交给 GeoDataCollection 聚合。"""
    files, configs = [], []

    for g in gps_cfg:
        obj = GPSData(g['file'], g.get('type', gfiletype))
        obj.load()
        files.append(obj)
        configs.append(g)

    for s in sar_cfg:
        obj = SARData(s['file'],
                      scale=s.get('scale', 1.0),
                      sigma=s.get('sigma'),
                      lamda=s.get('lamda'))
        obj.load()
        files.append(obj)
        configs.append(s)

    for l in lev_cfg:
        obj = LEVData(l['file'])
        obj.load()
        files.append(obj)
        configs.append(l)

    return GeoDataCollection(files, configs)


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def run(cfgfile):
    with open(cfgfile, 'r') as fid:
        cfg = yaml.safe_load(fid)

    G      = cfg.get('global', {})
    FAULTS = cfg.get('faults', [])
    GPSC   = cfg.get('gps', [])
    SARC   = cfg.get('insar', [])
    LEVC   = cfg.get('leveling', [])
    EXPC   = cfg.get('export', {})

    if not FAULTS:
        logging.error('No faults defined in config.')
        sys.exit(1)

    # 全局 GPS 格式：优先 global.gfiletype，否则取第一条 GPS 的 type
    gfiletype = G.get('gfiletype') or (GPSC[0].get('type') if GPSC else 'GMT2D')

    # 1) 公共原点 + 断层
    origin = G.get('origin') or build_origin(FAULTS)
    if len(origin) != 2:
        origin = build_origin(FAULTS)
    faults = build_faults(FAULTS, origin)

    # 2) 数据（逐文件加载后聚合，config 与文件一一配对）
    data = build_data(GPSC, SARC, LEVC, gfiletype)

    # 3) 反演（格林函数、平滑、ramp、权重、边界与扫描求解）
    inversion = MultiFaultSlipInversion(data, faults, origin=origin,
                                        config=G, export=EXPC)

    # 输出全部交给 MultiOutput：逐文件观测拟合、逐断层滑移模型、
    # HDF5/summary/L 曲线/GMT 滑移与归档
    outpath = EXPC.get('outpath', './new')
    output  = MultiOutput(data, faults, sol=inversion,
                          archdir=EXPC.get('archdir', outpath),
                          outdir=outpath)
    results = inversion.run(output=output)

    # 4) 汇总
    print('\n=== Summary (sm, WRSS_gps, WRSS_sar, roughness, Mo, Mw) ===')
    for row in results:
        print('%7.3f  %12.3f  %12.3f  %12.3f  %12.3e  %6.2f' % row)
    print('Slip models written to %s/' % outpath)


def main():
    parser = argparse.ArgumentParser(description='Multi-fault / multi-file slip inversion.')
    parser.add_argument('--cfgfile', type=str, required=True, help='path to YAML config')
    args = parser.parse_args()
    run(args.cfgfile)


if __name__ == '__main__':
    main()
