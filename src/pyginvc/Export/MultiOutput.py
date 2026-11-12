#!/usr/bin/env python
"""
MultiOutput -- 多断层 / 多数据文件的输出封装。

与 ``Export.Output`` 的区别：
    * 构造时接收 ``GeoDataCollection``（逐文件聚合的观测）和
      ``FaultCollection``（逐断层聚合的矩形断层），不要求它们是某一个
      ``GeoData`` / ``Fault`` 实例；
    * 所有落盘工作复用 ``Export.OutputUtils`` 的纯函数，本类只负责从
      集合对象里提取数组、按文件/断层切分，并拼出文件名；
    * 观测拟合按输入文件分别输出（GPS 每个文件一组、
      InSAR 每个 track 一组），滑移模型同时给出“全部断层合并”与
      “逐断层”两份。

目前仅支持矩形断层（``FaultCollection`` 的元素为 ``Patch.Fault``）。
三角网格的多断层输出仍需在 ``Output`` 侧扩展。
"""

import os
import re
import logging

import numpy as np

from pyginvc.Export import OutputUtils as ou

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")


def _safe_name(name, default):
    '''Sanitise a fault name for use inside a file name.'''
    if not name:
        return default
    return re.sub(r'[^0-9A-Za-z_.-]+', '_', str(name))


class MultiOutput(object):
    '''
    MultiOutput -- 输出多断层反演/正演结果。

    Input:
        data    = GeoDataCollection，已加载并聚合的观测数据
        faults  = FaultCollection，已加载并聚合的矩形断层
        sol     = 解算结果对象（MultiFaultSlipInversion / Forward 等），
                  需要提供 slip；输出观测拟合时还需要 dhat
        archdir = 归档目录；为 None 或与 outdir 相同时不移动文件
        outdir  = 输出目录，默认当前目录

    常用方法：
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

    # ------------------------------------------------------------------ #
    # 内部：组装纯函数需要的 dict / 数组
    # ------------------------------------------------------------------ #
    def _obs_dict(self):
        data = self.data
        return {
            'llh_gps': data.llh_gps,
            'llh_lev': data.llh_lev,
            'llh_sar': data.llh_sar,
            'unit'   : data.unit,
            'd_gps'  : data.d_gps,
            'd_lev'  : data.d_lev,
            'd_sar'  : data.d_sar,
            'ndim'   : data.ndim,
        }

    def _fault_dict(self):
        faults = self.faults
        origin = getattr(self.sol, 'origin', None)
        return {
            'origin'  : np.asarray(faults.origin if origin is None else origin),
            'flttype' : 'rectangle',
            'dis_geom': faults.dis_geom_grid,
            'geom'    : faults.geom_grid,
        }

    @staticmethod
    def _as_array(value):
        '''None -> empty array, so empty datasets can still be written.'''
        return np.array([]) if value is None else np.asarray(value)

    def _solution_dict(self):
        sol = self.sol
        if sol is None:
            raise ValueError('MultiOutput.sol is None: cannot write a solution.')
        return {
            'slip'  : self._as_array(sol.slip),
            'r'     : self._as_array(getattr(sol, 'r', None)),
            'misfit': self._as_array(getattr(sol, 'misfit', None)),
            'dhat'  : self._as_array(getattr(sol, 'dhat', None)),
        }

    def _slip_rows(self, slip=None):
        '''Return slip as a 2-D array of shape (nmodel, 3*nf).'''
        if slip is None:
            if self.sol is None or getattr(self.sol, 'slip', None) is None:
                raise ValueError('No slip given and the solution has none.')
            slip = self.sol.slip
        slip = np.asarray(slip)
        if slip.ndim == 1:
            slip = slip[None, :]
        return slip

    @staticmethod
    def _gps_sigma(gps, n, ndim):
        '''Per-station formal uncertainty from a GPSData covariance.'''
        cov = getattr(gps, 'cov_gps', None)
        if cov is None or np.size(cov) == 0:
            return None
        cov = np.asarray(cov)
        var = cov if cov.ndim == 1 else np.diag(cov)
        return np.sqrt(var).reshape(n, ndim)

    @staticmethod
    def _fault_vertices(flt):
        '''(nf, 4, 3) [lon, lat, depth] patch corners of one fault.'''
        felem = flt.FaultGeom2AllVertex()
        corners = (
            ('lt_lon', 'lt_lat', 'lt_dep'),
            ('rt_lon', 'rt_lat', 'rt_dep'),
            ('rb_lon', 'rb_lat', 'rb_dep'),
            ('lb_lon', 'lb_lat', 'lb_dep'),
        )
        return np.stack(
            [np.column_stack([np.asarray(felem[lon]),
                              np.asarray(felem[lat]),
                              np.asarray(felem[dep])])
             for lon, lat, dep in corners], axis=1)

    def _path(self, filename):
        if os.path.isabs(filename):
            return filename
        return os.path.join(self.outdir, filename)

    # ------------------------------------------------------------------ #
    # HDF5
    # ------------------------------------------------------------------ #
    def WriteHDF5(self, filename='solutions.h5'):
        '''Write observations, concatenated fault geometry and solution.'''
        path = self._path(filename)
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        ou.write_solution_h5(path, self._obs_dict(),
                             self._fault_dict(), self._solution_dict())
        logging.info('Solution is written to %s.' % path)
        return path

    # ------------------------------------------------------------------ #
    # 观测拟合
    # ------------------------------------------------------------------ #
    def WriteData(self, dhat=None, outdir=None):
        '''
        逐文件写出观测/模型/残差。

        Input:
            dhat  = 模型预测值，行序为 [GPS(逐文件), 水准, InSAR(逐 track)]
            outdir= 输出目录，默认 self.outdir

        Output:
            GPS 第 i 个文件 -> gps{i}_obs/mod/res[_up].gmtvec
            水准             -> lev_mod.xyz, lev_res.xyz
            InSAR 第 i 个 track -> sar_mod_{i}.xyz, sar_res_{i}.xyz
        '''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        if dhat is None:
            if self.sol is None or getattr(self.sol, 'dhat', None) is None:
                raise ValueError('No dhat given and the solution has none.')
            dhat = self.sol.dhat
        dhat = np.asarray(dhat).reshape(-1)

        ndim   = self.data.ndim
        offset = 0
        written = []

        # GPS：每个输入文件一组，行数 = n * ndim
        for i, (gps, _) in enumerate(self.data.gps_pairs()):
            n = len(gps.llh_gps)
            if n == 0:
                continue
            mod = dhat[offset:offset + n * ndim].reshape(n, ndim)
            obs = np.asarray(gps.d_gps).reshape(n, ndim)
            written += ou.write_gps_data(
                gps.llh_gps, obs, mod,
                sigma=self._gps_sigma(gps, n, ndim),
                station=getattr(gps, 'station_gps', None),
                outdir=outdir, prefix='gps%d' % (i + 1))
            offset += n * ndim

        # 水准：聚合输出
        n_lev = len(self.data.d_lev)
        if n_lev:
            mod = dhat[offset:offset + n_lev]
            written += ou.write_lev_data(self.data.llh_lev,
                                         self.data.d_lev, mod,
                                         outdir=outdir)
            offset += n_lev

        # InSAR：按 track 切分在同一函数内完成
        n_sar = len(self.data.d_sar)
        if n_sar:
            mod = dhat[offset:offset + n_sar]
            written += ou.write_sar_data(self.data.llh_sar,
                                         self.data.d_sar, mod,
                                         n_sar=self.data.n_sar,
                                         outdir=outdir)
            offset += n_sar

        logging.info('%d data-fit file(s) written to %s.' % (len(written), outdir))
        return written

    # ------------------------------------------------------------------ #
    # 滑移模型
    # ------------------------------------------------------------------ #
    def WriteSlipModel(self, slip=None, k=None, outdir=None,
                       prefix='slipmodel'):
        '''
        写矩形滑移模型。

        Input:
            slip   = (3*nf,) 或 (nsmooth, 3*nf)；缺省取 self.sol.slip
            k      = 只写第 k 个平滑因子（0 基）；缺省写全部
            outdir = 输出目录
            prefix = 文件名前缀

        Output:
            {prefix}_{kkk}.faultgeom             全部断层合并且
            {prefix}_{kkk}_{fault}.faultgeom     逐断层
            {prefix}_{kkk}_{fault}_rot.faultgeom 配置了 rake_beta 时的旋转结果
        '''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        single = slip is not None and np.asarray(slip).ndim == 1
        rows = self._slip_rows(slip)
        if single:
            # 已经是单个模型，k 只用于文件编号
            indices = [0 if k is None else k]
        elif k is None:
            indices = list(range(rows.shape[0]))
        else:
            rows = rows[k:k + 1]
            indices = [k]

        nf3  = 3 * self.faults.nf
        legs = list(self.faults.pairs())
        written = []

        for row, idx in zip(rows, indices):
            row  = np.asarray(row).reshape(-1)[:nf3]
            stem = '%s_%03d' % (prefix, idx + 1)

            # 全部断层合并（与 Export.Output 的 slipmodel_*.faultgeom 对应）
            fname = os.path.join(outdir, stem + '.faultgeom')
            ou.write_faultgeom(fname, self.faults.geom_grid, row)
            written.append(fname)

            # 逐断层
            off = 0
            for i, (flt, cfg) in enumerate(legs):
                islip = row[off:off + 3 * flt.nf].reshape(flt.nf, 3)
                fname = os.path.join(
                    outdir, '%s_%s.faultgeom'
                    % (stem, _safe_name((cfg or {}).get('name'),
                                        'fault%02d' % (i + 1))))
                ou.write_faultgeom(fname, flt.geom_grid, islip)
                written.append(fname)

                beta = ((cfg or {}).get('green') or {}).get('rake_beta', 0)
                if beta:
                    rname = fname[:-len('.faultgeom')] + '_rot.faultgeom'
                    ou.write_faultgeom(rname, flt.geom_grid,
                                       ou.rotate_slip(islip, beta))
                    written.append(rname)
                off += 3 * flt.nf

        return written

    def WriteGMTSlip(self, slip=None, outdir=None):
        '''写 GMT 多段面片形式的滑移分布（取最后一个平滑因子）。'''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        row = np.asarray(self._slip_rows(slip)[-1]).reshape(-1)

        vertices, values, off = [], [], 0
        for flt, _ in self.faults.pairs():
            vertices.append(self._fault_vertices(flt))
            values.append(row[off:off + 3 * flt.nf].reshape(flt.nf, 3))
            off += 3 * flt.nf

        vertices = np.vstack(vertices)
        slip3    = np.vstack(values)
        slip_mag = np.sqrt(slip3[:, 0] ** 2 + slip3[:, 1] ** 2)

        header = '# Rectangle fault slip model'
        fmt    = '%10.4f %10.4f %10.4f'
        written = []
        for fname, value in (('ss_slip_llh.gmtlin', slip3[:, 0]),
                             ('ds_slip_llh.gmtlin', slip3[:, 1]),
                             ('ts_slip_llh.gmtlin', slip_mag)):
            path = os.path.join(outdir, fname)
            ou.write_gmt_slip(path, vertices, value, header=header, fmt=fmt)
            written.append(path)

        logging.info('Fault slip distribution is output in geography system.')
        return written

    # ------------------------------------------------------------------ #
    # 统计
    # ------------------------------------------------------------------ #
    def WriteSummary(self, filename='summary'):
        '''写数据/断层/反演汇总（断层参数逐条列出）。'''
        sol = self.sol
        if sol is None:
            raise ValueError('MultiOutput.sol is None: cannot write summary.')

        path = self._path(filename)
        ou.write_summary(
            path,
            ndim       = self.data.ndim,
            n_gps      = len(self.data.llh_gps),
            n_sar      = len(self.data.llh_sar),
            n_lev      = len(self.data.llh_lev),
            flttype    = 'rectangle',
            nstrike    = list(self.faults.nsegs),
            ndip       = list(self.faults.ndeps),
            smo_facts  = self._as_array(sol.smo_facts),
            misfit_gps = self._as_array(sol.misfit_gps),
            misfit_sar = self._as_array(sol.misfit_sar),
            n_d_gps    = len(self.data.d_gps),
            n_d_sar    = len(self.data.d_sar),
            moment     = self._as_array(sol.moment),
        )
        return path

    def WriteLCurve(self, filename='ruff_misfit'):
        '''写 roughness / misfit / 平滑因子 / 矩 的权衡表。'''
        sol = self.sol
        if sol is None:
            raise ValueError('MultiOutput.sol is None: cannot write L-curve.')

        path = self._path(filename)
        smo_facts = self._as_array(sol.smo_facts)
        ok = ou.write_lcurve(path,
                             self._as_array(sol.ruff),
                             self._as_array(sol.misfit),
                             smo_facts,
                             self._as_array(sol.moment))
        if not ok:
            logging.warning('Only %d smoothing factors.' % len(smo_facts))
            return None
        logging.info('Roughness vs misfit is printed into %s.' % path)
        return path

    # ------------------------------------------------------------------ #
    # 单个平滑因子
    # ------------------------------------------------------------------ #
    def WriteSmoothing(self, k, result, outdir=None):
        '''
        输出第 k 个平滑因子的结果。

        Input:
            k      = 平滑因子序号（0 基）
            result = MultiFaultSlipInversion.solve() 返回的 dict
                     （至少含 slip 与 pred）
        '''
        outdir = outdir or self.outdir
        self.WriteSlipModel(result['slip'], k=k, outdir=outdir)

        subdir = os.path.join(outdir, 'sm%03d' % k)
        os.makedirs(subdir, exist_ok=True)
        self.WriteData(result['pred'], outdir=subdir)
        return subdir

    # ------------------------------------------------------------------ #
    # 总入口 / 归档
    # ------------------------------------------------------------------ #
    def OutputSolution(self):
        '''写出 HDF5、汇总、模型与观测拟合，并归档。'''
        os.makedirs(self.outdir, exist_ok=True)
        self.WriteHDF5()
        self.WriteSummary()
        self.WriteSlipModel()
        self.WriteGMTSlip()
        self.WriteLCurve()
        self.WriteData()
        self.archive_outfile()

    def archive_outfile(self):
        '''将 outdir 中匹配默认 patterns 的文件移动到 archdir。'''
        if not self.archdir:
            logging.error('The input archive directory is null')
            return []
        return ou.archive_files(self.archdir,
                                ou.DEFAULT_ARCHIVE_PATTERNS,
                                src_dir=self.outdir)
