#!/usr/bin/env python
"""
MultiOutput -- 多断层 / 多数据文件的输出适配器。

``Output`` 现在负责所有公共编排（HDF5/summary/滑移模型/GMT/观测拟合/归档），
``MultiOutput`` 只是它的子类，覆写下面几个钩子：

    _obs_blocks()      逐 GPS 文件 / 水准文件 / InSAR track 生成观测块
    _fault_groups()    逐断层生成 FaultGroup（含名称与 rake_beta）
    _solution_view()   把 MultiFaultInversion 的解映射成 SolutionView
    _origin()          多断层共用的参考原点
    _summary_counts()  逐断层的 nsegs/ndeps

所以单断层与多断层走完全相同的落盘代码，不会再出现两套实现分叉。
"""

import logging
import re

import numpy as np

from pyginvc.Export.Output import Output, ObsBlock, FaultGroup, SolutionView

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")


def _safe_name(name, default):
    '''Sanitise a fault name for use inside a file name.'''
    if not name:
        return default
    return re.sub(r'[^0-9A-Za-z_.-]+', '_', str(name))


class MultiOutput(Output):
    '''
    MultiOutput -- 输出多断层反演/正演结果。

    Input:
        data    = GeoDataCollection，已加载并聚合的观测数据
        faults  = FaultCollection，已加载并聚合的矩形断层
        sol     = 解算结果对象（MultiFaultInversion / Forward 等），
                  需要提供 slip；输出观测拟合时还需要 dhat
        archdir = 归档目录；为 None 或与 outdir 相同时不移动文件
        outdir  = 输出目录，默认当前目录

    常用方法（继承自 Output）：
        OutputSolution()      一次性写出 HDF5 + summary + 模型 + 观测拟合
        WriteHDF5()           solutions.h5
        WriteSummary()        summary
        WriteLCurve()         ruff_misfit
        WriteSlipModel()      逐平滑因子、逐断层的 .faultgeom
        WriteGMTSlip()        ss/ds/ts_slip_llh.gmtlin
        WriteData(dhat)       逐文件的观测/模型/残差
        WriteSmoothing(k, r)  单个平滑因子的模型 + 观测拟合
        archive_outfile()     按默认 patterns 归档
    '''

    def __init__(self, data, faults, sol=None, archdir=None, outdir='.'):
        self.data    = data
        self.faults  = faults
        self.sol     = sol
        self.archdir = archdir
        self.outdir  = outdir or '.'

        # 基类钩子不读取这些属性，保留是为了兼容外部代码
        self.flt     = None
        self.green   = None
        self.flttype = 'rectangle'

    # ------------------------------------------------------------------ #
    # 适配钩子
    # ------------------------------------------------------------------ #
    def _obs_blocks(self):
        '''逐文件生成 ObsBlock，顺序为 GPS -> 水准 -> InSAR。'''
        blocks = []
        ndim = self.data.ndim

        for i, (gps, _) in enumerate(self.data.gps_pairs()):
            n = len(gps.llh_gps)
            if n == 0:
                continue
            blocks.append(ObsBlock(
                'gps', 'gps%d' % (i + 1), gps.llh_gps, gps.d_gps,
                ndim=ndim, sigma=self._gps_sigma(gps, n, ndim),
                station=getattr(gps, 'station_gps', None)))

        levs = list(self.data.lev_pairs())
        for i, (lev, _) in enumerate(levs):
            n = len(np.asarray(lev.llh_lev))
            if n == 0:
                continue
            prefix = 'lev%d' % (i + 1) if len(levs) > 1 else 'lev'
            blocks.append(ObsBlock('lev', prefix, lev.llh_lev, lev.d_lev,
                                   ndim=1))

        if len(self.data.d_sar):
            blocks.append(ObsBlock(
                'sar', 'sar', self.data.llh_sar, self.data.d_sar, ndim=1,
                unit=getattr(self.data, 'unit', None),
                n_sar=self._n_sar()))

        return blocks

    def _fault_groups(self):
        '''逐断层生成 FaultGroup；名称用于逐断层文件名。'''
        groups = []
        for i, (flt, cfg) in enumerate(self.faults.pairs()):
            name = _safe_name((cfg or {}).get('name'),
                              'fault%02d' % (i + 1))
            beta = ((cfg or {}).get('green') or {}).get('rake_beta', 0) or 0
            groups.append(FaultGroup(flt, name, 'rectangle', beta))
        return groups

    def _solution_view(self):
        sol = self.sol
        if sol is None:
            return None
        return SolutionView(
            slip=getattr(sol, 'slip', None),
            dhat=getattr(sol, 'dhat', None),
            r=getattr(sol, 'r', None),
            misfit=getattr(sol, 'misfit', None),
            misfit_gps=getattr(sol, 'misfit_gps', None),
            misfit_sar=getattr(sol, 'misfit_sar', None),
            smo_facts=getattr(sol, 'smo_facts', None),
            moment=getattr(sol, 'moment', None),
            ruff=getattr(sol, 'ruff', None),
            ramp_sar=getattr(sol, 'ramp_sar', None),
            sar_ramp_flag=getattr(sol, 'sar_ramp_flag', None),
            forward=False)

    def _origin(self):
        origin = getattr(self.sol, 'origin', None)
        if origin is None:
            origin = self.faults.origin
        return np.asarray(origin)

    def _flttype(self):
        return 'rectangle'

    def _summary_counts(self, groups):
        return list(self.faults.nsegs), list(self.faults.ndeps)

    # ------------------------------------------------------------------ #
    # 工具
    # ------------------------------------------------------------------ #
    @staticmethod
    def _gps_sigma(gps, n, ndim):
        '''Per-station formal uncertainty from a GPSData covariance.'''
        cov = getattr(gps, 'cov_gps', None)
        if cov is None or np.size(cov) == 0:
            return None
        cov = np.asarray(cov)
        var = cov if cov.ndim == 1 else np.diag(cov)
        return np.sqrt(var).reshape(n, ndim)
