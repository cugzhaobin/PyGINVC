#!/usr/bin/env python
"""
MultiFaultSlipInversion.py -- 多断层 + 多数据文件的滑动反演。

设计原则：
    * 每条断层建一个 Fault 对象，再把 dis_geom_grid / geom_grid 纵向拼接；
    * 每个数据文件建一个数据对象，再由 GeoDataCollection 聚合；
    * 每条断层建一个 Okada，各自算格林函数后横向拼接（支持逐断层 greentype/nu/modulus/rake_beta）；
    * 每条断层建一个 RectPlaneLap，再用 block_diag 拼成块对角平滑；
    * 权重、ramp、边界在类内自拼，求解/残差/边界复用 BaseInversion 的接口。

命令行入口见 pyginvc/scripts/MultiFaultSlipInversion.py。
"""

import os
import logging
import numpy as np
from scipy import linalg

from pyginvc.libs import geotools as gt
from pyginvc.Greens.Okada import Okada
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap
from pyginvc.Inversion.BaseInversion import BaseInversion
from pyginvc.Export.MultiOutput import MultiOutput

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')


# --------------------------------------------------------------------------- #
# 各步骤构建
# --------------------------------------------------------------------------- #
def build_green(faults, data, origin):
    """逐断层算格林函数，横向拼接。GPS/水准复用 make_gps_green，InSAR 自拼以绕开 n_sar bug。"""
    xy_gps = np.array([gt.llh2utm(p, origin) for p in data.llh_gps]).reshape(-1, 2)
    xy_lev = np.array([gt.llh2utm(p, origin) for p in data.llh_lev]).reshape(-1, 2)
    xy_sar = np.array([gt.llh2utm(p, origin) for p in data.llh_sar]).reshape(-1, 2)

    G_blocks, Gsar_blocks = [], []
    for flt, f in faults.pairs():
        g = dict(f['green'])
        g['greenfile'] = 'NONE'            # BaseGreen 必需键，这里不落盘
        ok = Okada(flt, data, g)
        ss, ds, op = g['greentype']
        nu = g['nu']

        # GPS + 水准（每点 ndim 个分量）
        G = np.empty((0, 3 * flt.nf))
        if xy_gps.size:
            G = ok.make_gps_green(flt.dis_geom_grid, xy_gps, nu, ss, ds, op, data.ndim)
        if xy_lev.size:
            Glev = ok.make_gps_green(flt.dis_geom_grid, xy_lev, nu, ss, ds, op, data.ndim)
            G = np.vstack((G, Glev))

        # InSAR（每点 1 个 LOS 值）
        Gsar = np.empty((0, 3 * flt.nf))
        if xy_sar.size:
            Gsar = ok.make_sar_green(data.unit, flt.dis_geom_grid, xy_sar, nu, ss, ds, op)

        # 逐断层 rake_beta 旋转（内联，避免 rotate_greens 写 green_new.h5 的副作用）
        beta = g.get('rake_beta', 0)
        if beta:
            c, s = np.cos(np.deg2rad(beta)), np.sin(np.deg2rad(beta))
            for M in (G, Gsar):
                for i in range(flt.nf):
                    u = M[:, 3 * i].copy()
                    v = M[:, 3 * i + 1].copy()
                    M[:, 3 * i]     = u * c + v * s
                    M[:, 3 * i + 1] = -u * s + v * c

        G_blocks.append(G)
        Gsar_blocks.append(Gsar)

    G     = np.hstack(G_blocks) if G_blocks else np.empty((0, 0))
    G_sar = np.hstack(Gsar_blocks) if Gsar_blocks else np.empty((0, 0))
    return G, G_sar


def build_lap(faults, lapmethod):
    """逐断层建 Laplacian，再拼成块对角。nsegs/ndeps < 2 时退回 Tikhonov(单位阵)。"""
    laps = []
    for flt, f in faults.pairs():
        method = lapmethod if (flt.nsegs >= 2 and flt.ndeps >= 2) else '4'
        lap = RectPlaneLap(flt.nsegs, flt.ndeps, f['bcs'], method)
        lap.build_laps()
        laps.append(lap.G_lap)
    return linalg.block_diag(*laps)


def build_ramp(data, origin):
    """逐文件/逐 track 构建 ramp 设计矩阵，返回 (G_ramp, n_ramp)。
    行序 = [GPS(ndim), 水准(1), InSAR(1)]，未启用 ramp 的文件用 0 列块占位。"""
    ndim = data.ndim
    xy_gps = np.array([gt.llh2utm(p, origin) for p in data.llh_gps]).reshape(-1, 2)
    xy_sar = np.array([gt.llh2utm(p, origin) for p in data.llh_sar]).reshape(-1, 2)

    blocks = []
    n_ramp = 0

    # GPS：每个文件一个块，未启用 ramp 时用 0 列块占位
    sta0 = 0
    for gps, g in data.gps_pairs():
        g = g or {}
        n = len(gps.llh_gps)
        if g.get('ramp', False):
            ok = Okada(None, gps, {'greenfile': 'NONE'})
            M = ok.make_green_gps_ramp(xy_gps[sta0:sta0 + n], ndim,
                                       g.get('ramp_method', 1))
            blocks.append(M)
            n_ramp += M.shape[1]
        else:
            blocks.append(np.zeros((n * ndim, 0)))
        sta0 += n

    # 水准：不建 ramp，只占行
    blocks.append(np.zeros((len(data.d_lev), 0)))

    # InSAR：每个 track 一个块
    sta0 = 0
    for sar, s in data.sar_pairs():
        s = s or {}
        n = len(sar.llh_sar)
        if s.get('ramp', False):
            ok = Okada(None, sar, {'greenfile': 'NONE'})
            M = ok.make_green_sar_ramp(xy_sar[sta0:sta0 + n])
            blocks.append(M)
            n_ramp += M.shape[1]
        else:
            blocks.append(np.zeros((n, 0)))
        sta0 += n

    if n_ramp == 0:
        return np.empty((0, 0)), 0
    return linalg.block_diag(*blocks), n_ramp


def build_bounds(faults, n_ramp):
    """逐断层构造滑移上下界（复用 BaseInversion.set_slip_bounds），ramp 参数不加约束。"""
    bl_parts, bu_parts = [], []
    for flt, f in faults.pairs():
        b = dict((f or {}).get('bound', {}))
        # BaseInversion.set_slip_bounds 缺省为 0 界，这里保持原先的无穷界
        b.setdefault('ss_range', [-np.inf, np.inf])
        b.setdefault('ds_range', [-np.inf, np.inf])
        b.setdefault('op_range', [-np.inf, np.inf])
        bu, bl = BaseInversion.set_slip_bounds(flt.nsegs, flt.ndeps, 0, **b)
        bl_parts.append(np.asarray(bl).reshape(-1))
        bu_parts.append(np.asarray(bu).reshape(-1))

    bl = np.concatenate(bl_parts)
    bu = np.concatenate(bu_parts)
    if n_ramp > 0:
        bl = np.hstack([bl, np.full(n_ramp, -np.inf)])
        bu = np.hstack([bu, np.full(n_ramp,  np.inf)])
    return bl, bu


def total_moment(faults, slip):
    """逐断层累加地震矩（复用 Fault.moment）。"""
    Mo = 0.0
    off = 0
    for flt, f in faults.pairs():
        s = slip[off:off + 3 * flt.nf].reshape(-1, 3)
        Mo += flt.moment(s, float((f or {}).get('green', {}).get('modulus', 3e10)))[0]
        off += 3 * flt.nf
    Mw = (2.0 / 3.0) * np.log10(Mo) - 6.067 if Mo > 0 else 0.0
    return Mo, Mw


def write_data_fit(data, pred, outdir):
    """每个观测文件输出对应的模拟值(*_mod)与残差(*_res)。

    保留此函数作为向后兼容的入口，实际写出由 ``Export.MultiOutput``
    委托给 ``OutputUtils`` 完成。
    """
    return MultiOutput(data, None).WriteData(pred, outdir=outdir)


# --------------------------------------------------------------------------- #
# 反演类
# --------------------------------------------------------------------------- #
class MultiFaultSlipInversion(BaseInversion):
    '''
    MultiFaultSlipInversion -- 多断层滑动反演（继承 BaseInversion）。

    接收已加载并聚合的观测数据与断层，内部完成格林函数、平滑矩阵、ramp、
    权重、边界的构建，并按平滑因子扫描求解与输出。

    Input:
        data    = GeoDataCollection，已加载并聚合的观测数据
        faults  = FaultCollection，已加载并聚合的断层
        origin  = 公共参考原点 [lat, lon]；缺省取 faults.origin
        config  = 全局配置（YAML 的 global 段）
        export  = 输出配置（YAML 的 export 段）

    复用 BaseInversion 提供的 get_smoothing_factors / solve_least_squares /
    compute_residuals / set_slip_bounds；由于 GeoDataCollection 已按文件
    处理权重，这里重写 get_weight 返回其聚合权重的对角矩阵。
    '''

    def __init__(self, data, faults, origin=None, config=None, export=None):
        config = config or {}
        dict_weight = {
            'smoothfactor': config.get('smoothfactor', [0.1, 0.5, 0.1]),
        }
        dict_bound = {
            'bound_switch': config.get('bound_switch', True),
        }
        super().__init__(data, dict_weight, dict_bound)

        self.faults = faults
        self.origin = list(origin) if origin is not None else list(faults.origin)
        self.config = config
        self.export = export or {}

        self.prepared = False
        self.Gmat = self.G_sar = None
        self.G_data = None
        self.G_lap = None
        self.G_ramp = None
        self.n_ramp = 0
        self.bl = self.bu = None
        self.d_geod = None
        self.d_sar = None
        self.n_geod = 0
        self.n_sar = 0
        self.nf3 = 0
        self.geom_grid = None

    def get_weight(self):
        '''
        Return diagonal weight matrices for (GPS + leveling) and InSAR data.

        GeoDataCollection already applies the per-file weights and supplies the
        variances, so reuse its aggregation instead of the BaseInversion
        implementation (which expects GeoData's W_sar attribute).
        '''
        return self.data.get_weights()

    def _ensure_prepared(self):
        if not self.prepared:
            self.prepare()

    def prepare(self):
        '''构建反演所需的全部矩阵与向量，可重复调用。'''
        data, faults = self.data, self.faults

        # 格林函数（逐断层横向拼接）
        self.Gmat, self.G_sar = build_green(faults, data, self.origin)
        self.G_data = np.vstack([self.Gmat, self.G_sar])

        # 平滑矩阵（逐断层块对角）
        self.G_lap = build_lap(faults, self.config.get('lapmethod', '1'))

        # 权重 / ramp
        self.W, self.WSAR = self.get_weight()
        self.G_ramp, self.n_ramp = build_ramp(data, self.origin)

        # 数据向量与维度
        d_gps, d_lev, d_sar = data.d_gps, data.d_lev, data.d_sar
        self.d_geod = np.hstack([x for x in (d_gps, d_lev) if len(x) > 0])
        self.d_sar  = d_sar
        self.n_geod = len(self.d_geod)
        self.n_sar  = len(d_sar)
        self.nf3    = 3 * faults.nf
        self.geom_grid = faults.geom_grid

        # 边界
        self.bound_switch = self.dict_bound.get('bound_switch', True)
        if self.bound_switch:
            self.bl, self.bu = build_bounds(faults, self.n_ramp)

        self.prepared = True
        return self

    def solve(self, smo):
        '''给定平滑因子求解一次，返回解及各项统计量。'''
        self._ensure_prepared()

        G2I_slip = np.vstack([self.W @ self.Gmat,
                              self.WSAR @ self.G_sar,
                              smo * self.G_lap])
        d2I = np.hstack([self.W @ self.d_geod,
                         self.WSAR @ self.d_sar,
                         np.zeros(self.nf3)])

        if self.n_ramp > 0:
            G2I = np.hstack([G2I_slip,
                             np.vstack([self.G_ramp,
                                        np.zeros((self.nf3, self.n_ramp))])])
        else:
            G2I = G2I_slip

        x = self.solve_least_squares(G2I, d2I, self.bl, self.bu)
        slip = x[:self.nf3]

        # 残差（加权空间，复用 BaseInversion.compute_residuals）
        G_pred = (np.hstack([self.G_data, self.G_ramp])
                  if self.n_ramp > 0 else self.G_data)
        d2R_parts = []
        if self.n_geod:
            d2R_parts.append(self.W @ self.d_geod)
        if self.n_sar:
            d2R_parts.append(self.WSAR @ self.d_sar)
        d2R = np.hstack(d2R_parts) if d2R_parts else np.array([])

        r_vec, pred, misfit_gps, misfit_sar, _ = self.compute_residuals(
            d2R, G_pred, x, self.W, self.WSAR, self.n_geod, self.n_sar)

        g_lap = slip @ self.G_lap
        ruff  = float(g_lap @ g_lap)
        Mo, Mw = total_moment(self.faults, slip)

        return {
            'smo'       : smo,
            'x'         : x,
            'slip'      : slip,
            'pred'      : pred,
            'dhat'      : pred,
            'r'         : r_vec,
            'misfit_gps': float(misfit_gps),
            'misfit_sar': float(misfit_sar),
            'roughness' : ruff,
            'Mo'        : Mo,
            'Mw'        : Mw,
        }

    def run(self, smo_facts=None, output=None):
        '''
        按平滑因子扫描求解并输出，返回汇总行列表。

        Input:
            smo_facts = 平滑因子数组；缺省用 get_smoothing_factors()
            output    = 输出对象（通常为 Export.MultiOutput）；缺省按
                        self.export 里的 outpath/archdir 自建一个
        '''
        self._ensure_prepared()

        if smo_facts is None:
            smo_facts = self.get_smoothing_factors()
        smo_facts = np.asarray(smo_facts).reshape(-1)
        if smo_facts.size == 0:
            raise ValueError('No smoothing factors to invert for.')

        outpath = self.export.get('outpath', './new')
        os.makedirs(outpath, exist_ok=True)

        if output is None:
            output = MultiOutput(self.data, self.faults, sol=self,
                                 archdir=self.export.get('archdir', outpath),
                                 outdir=outpath)

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
            logging.info('sm=%7.3f  WRSS_gps=%10.3f  WRSS_sar=%10.3f  roughness=%10.3f  Mw=%6.2f'
                         % (smo, r['misfit_gps'], r['misfit_sar'],
                            r['roughness'], r['Mw']))

        # 与 BaseInversion 的解决方案属性保持一致
        self.smo_facts  = np.asarray(smo_facts)
        self.slip       = np.asarray(slip)
        self.misfit_gps = np.asarray(misfit_gps)
        self.misfit_sar = np.asarray(misfit_sar)
        self.misfit     = self.misfit_gps + self.misfit_sar
        self.ruff       = np.asarray(rough)
        self.smoothness = self.ruff * self.smo_facts ** 2
        self.moment     = np.asarray(moment)
        self.dhat       = np.asarray(last['dhat'])
        self.r          = np.asarray(last['r'])

        # 全局输出：HDF5、汇总、L 曲线、GMT 滑移（取最后一个平滑因子）
        output.WriteHDF5()
        output.WriteSummary()
        output.WriteLCurve()
        output.WriteGMTSlip()
        output.archive_outfile()
        return results
