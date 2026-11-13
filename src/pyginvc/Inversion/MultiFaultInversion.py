#!/usr/bin/env python
import os
import logging
import numpy as np
from scipy import linalg

from pyginvc.libs import geotools as gt
from pyginvc.Greens.BaseGreen import BaseGreen
from pyginvc.Greens.Okada import Okada
from pyginvc.Greens.EDCMP import EDCMP
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap
from pyginvc.Inversion.BaseInversion import BaseInversion
from pyginvc.Export.MultiOutput import MultiOutput

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')


class MultiFaultInversion(BaseInversion):
    '''
    多断层、多数据源的滑动反演，是 GeoInversion 的升级版：
    接收已加载并聚合的观测数据与断层，内部完成格林函数、平滑矩阵、ramp、
    权重、边界的构建，并按平滑因子扫描求解与输出。

    Input:
        data    = GeoDataCollection，已加载并聚合的观测数据
        faults  = FaultCollection，已加载并聚合的断层
        origin  = 公共参考原点 [lat, lon]；缺省取 faults.origin
        config  = 全局配置（YAML 的 global 段）
        export  = 输出配置（YAML 的 export 段）

    复用 BaseInversion 提供的 get_smoothing_factors / solve_least_squares /
    compute_residuals / set_slip_bounds。权重由
    GeoDataCollection.get_weights() 按文件聚合，在 prepare() 里一次性取得，
    不再使用 BaseInversion.get_weight()（它面向单断层 GeoData）。
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
        self.outpath = self.export.get('outpath', './')

        self.prepared = False
        self.Gmat = self.G_sar = None
        self.G_lap  = None
        self.G_data = None
        self.G_ramp = None
        self.n_ramp = 0
        # 每个 InSAR 文件是否估计了 ramp，以及最终模型的 ramp 系数与
        # 其对观测的影响（长度 = InSAR 总点数），供绘图使用
        self.sar_ramp_flag = None
        self.ramp = None
        self.ramp_sar = None
        self.bl = self.bu = None
        self.d_geod = None
        self.d_sar = None
        self.n_geod = 0
        self.n_sar = 0
        self.nf3 = 0
        self.geom_grid = None

    def _ensure_prepared(self):
        if not self.prepared:
            self.prepare()


    def prepare(self):
        '''构建反演所需的全部矩阵与向量，可重复调用。'''
        faults = self.faults

        # 数据向量与维度必须先于 build_green 确定：格林函数的行数直接取自
        # n_geod/n_sar，不再按 llh 推算第二遍。
        d_gps, d_lev, d_sar = self.data.d_gps, self.data.d_lev, self.data.d_sar
        geod_parts = [x for x in (d_gps, d_lev) if len(x) > 0]
        self.d_geod = (np.hstack(geod_parts) if geod_parts
                       else np.array([]))
        self.d_sar  = d_sar
        self.n_geod = len(self.d_geod)
        self.n_sar  = len(d_sar)
        self.nf3    = 3 * faults.nf

        # 格林函数（逐断层横向拼接）
        self.Gmat, self.G_sar = self.build_green()
        self.G_data = np.vstack([self.Gmat, self.G_sar])

        # 平滑矩阵（逐断层块对角）
        self.G_lap = self.build_laps()

        # 权重 / ramp
        self.W, self.WSAR = self.data.get_weights()
        self.G_ramp, self.n_ramp = self.build_ramp()

        # 加权矩阵/向量只依赖数据与权重，按平滑因子循环求解时直接复用。
        self._Gw      = (self.W @ self.Gmat if self.n_geod
                         else np.zeros((0, self.nf3)))
        self._Gsarw   = (self.WSAR @ self.G_sar if self.n_sar
                         else np.zeros((0, self.nf3)))
        self._d_geodw = (self.W @ self.d_geod if self.n_geod
                         else np.array([]))
        self._d_sarw  = (self.WSAR @ self.d_sar if self.n_sar
                         else np.array([]))
        self._d2R     = np.hstack([self._d_geodw, self._d_sarw])
        self._d2I     = np.hstack([self._d2R, np.zeros(self.nf3)])

        if self.n_ramp > 0:
            # Ramp 列必须与它们描述的数据使用相同的行权重；否则 ramp 会在
            # 未加权尺度上吸收噪声并偏置滑动估计。
            ramp_parts = []
            if self.n_geod:
                ramp_parts.append(self.W @ self.G_ramp[:self.n_geod])
            if self.n_sar:
                ramp_parts.append(self.WSAR @ self.G_ramp[self.n_geod:])
            self._G_ramp_w = (np.vstack(ramp_parts) if ramp_parts
                              else np.zeros((0, self.n_ramp)))
        else:
            self._G_ramp_w = np.zeros((0, 0))

        self.geom_grid = faults.geom_grid

        # 边界
        self.bound_switch = self.dict_bound.get('bound_switch', True)
        if self.bound_switch:
            self.bl, self.bu = self.build_bounds()

        self.prepared = True
        return self

    def solve(self, smo):
        '''给定平滑因子求解一次，返回解及各项统计量。'''
        self._ensure_prepared()

        G2I_slip = np.vstack([self._Gw,
                              self._Gsarw,
                              smo * self.G_lap])
        d2I = self._d2I

        if self.n_ramp > 0:
            G2I = np.hstack([G2I_slip,
                             np.vstack([self._G_ramp_w,
                                        np.zeros((self.nf3, self.n_ramp))])])
        else:
            G2I = G2I_slip

        x = self.solve_least_squares(G2I, d2I, self.bl, self.bu)
        slip = x[:self.nf3]

        # ramp 系数及其对 InSAR 观测的贡献（无 ramp 时为全 0）
        ramp_hat = x[self.nf3:] if self.n_ramp > 0 else np.array([])
        if self.n_ramp > 0 and self.n_sar > 0:
            ramp_sar = self.G_ramp[self.n_geod:] @ ramp_hat
        else:
            ramp_sar = np.zeros(self.n_sar)

        # 预测与残差（加权空间，复用 BaseInversion.compute_residuals）。
        # pred 由未加权的 G_data/G_ramp 直接给出，传入 dhat 可省掉
        # compute_residuals 里重复的 G2R @ x。
        pred = self.G_data @ slip
        if self.n_ramp > 0:
            pred = pred + self.G_ramp @ ramp_hat
        r_vec, _, misfit_gps, misfit_sar, _ = self.compute_residuals(
            self._d2R, None, x, self.W, self.WSAR,
            self.n_geod, self.n_sar, dhat=pred)

        # The penalty in G2I is smo*G_lap @ slip, so the roughness is
        # ||G_lap @ slip||^2 (not ||slip @ G_lap||^2 = s^T L L^T s, which
        # differs when L is not symmetric, e.g. for mirrored bcs=1 edges).
        lap_slip = self.G_lap @ slip
        ruff     = float(lap_slip @ lap_slip)
        Mo, Mw = self.faults.total_moment(
            slip, default_modulus=self._default_modulus())

        return {
            'smo'       : smo,
            'x'         : x,
            'slip'      : slip,
            'dhat'      : pred,
            'r'         : r_vec,
            'misfit_gps': float(misfit_gps),
            'misfit_sar': float(misfit_sar),
            'roughness' : ruff,
            'Mo'        : Mo,
            'Mw'        : Mw,
            'ramp'      : ramp_hat,
            'ramp_sar'  : ramp_sar,
        }

    def _default_modulus(self):
        '''全局默认剪切模量：顶层 modulus 优先，其次 green.modulus。'''
        mod = self.config.get('modulus')
        if mod is None:
            mod = (self.config.get('green') or {}).get('modulus', 3e10)
        return mod

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

        outpath = self.export.get('outpath', './')
        os.makedirs(outpath, exist_ok=True)

        if output is None:
            output = MultiOutput(self.data,
                                 self.faults,
                                 sol = self,
                                 archdir = self.export.get('archdir', outpath),
                                 outdir = outpath)

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
        self.ramp       = np.asarray(last['ramp'])
        self.ramp_sar   = np.asarray(last['ramp_sar'])

        # 全局输出：HDF5、汇总、L 曲线、GMT 滑移（取最后一个平滑因子）
        output.WriteHDF5()
        output.WriteSummary()
        output.WriteLCurve()
        output.WriteGMTSlip()
        output.archive_outfile()
        return results


    def build_bounds(self):
        """
        set bounds for parameters to be estimated.
        """

        bl_parts, bu_parts = [], []
        for flt, cfg in self.faults.pairs():
            # ``bound:`` without a mapping (YAML null) means "use defaults".
            b = dict((cfg or {}).get('bound') or {})
            b.setdefault('ss_range', [-np.inf, np.inf])
            b.setdefault('ds_range', [-np.inf, np.inf])
            b.setdefault('op_range', [-np.inf, np.inf])
            bu, bl = BaseInversion.set_slip_bounds(flt.nsegs, flt.ndeps, 0, **b)
            bl_parts.append(np.asarray(bl).reshape(-1))
            bu_parts.append(np.asarray(bu).reshape(-1))

        bl = np.concatenate(bl_parts)
        bu = np.concatenate(bu_parts)
        if self.n_ramp > 0:
            bl = np.hstack([bl, np.full(self.n_ramp, -np.inf)])
            bu = np.hstack([bu, np.full(self.n_ramp,  np.inf)])
        return bl, bu


    def build_laps(self):
        """逐断层建 Laplacian，再拼成块对角。nsegs/ndeps < 2 时退回 Tikhonov(单位阵)。"""
        laps = []
        for flt, cfg in self.faults.pairs():
            cfg       = cfg or {}
            lapmethod = str(self.config.get('lapmethod', '1')).strip()
            if lapmethod not in ('1', '4'):
                raise ValueError(
                    "Unsupported lapmethod %r; use '1' (Laplacian) or "
                    "'4' (Tikhonov)." % (lapmethod,))
            method    = lapmethod if (flt.nsegs >= 2 and flt.ndeps >= 2) else '4'
            lap       = RectPlaneLap(flt.nsegs, flt.ndeps,
                                     cfg.get('bcs', [0, 0, 0, 0]), method)
            lap.build_laps()
            laps.append(lap.G_lap)
        return linalg.block_diag(*laps)


    def build_green(self):
        """逐断层构建（或复用缓存）格林函数，并按断层横向拼接。

        全局 global.green 段与逐断层 green 段合并，逐断层优先；未给出的
        greentype/nu/method/modulus 用安全默认值补齐，避免 KeyError。
        method 支持 'okada'（均匀半空间）与 'wang'（R. J. Wang 层状模型，
        通过 EDCMP 计算，需要 edgrn 生成的格林函数文件）。
        每个断层使用独立的缓存文件，防止多条断层共用同一文件而读到彼此
        错位的矩阵。
        """
        G_blocks, Gsar_blocks = [], []
        os.makedirs(self.outpath, exist_ok=True)

        # n_geod/n_sar 由 prepare() 按实际数据向量长度给出。
        n_geod = self.n_geod
        n_sar  = self.n_sar

        for i, (flt, cfg) in enumerate(self.faults.pairs()):
            g   = self._green_config(i, cfg)
            nf3 = 3 * flt.nf
            # Use the origin stored on the Fault itself.  Its patch geometry
            # (dis_geom_grid) was computed relative to that origin when it was
            # loaded; overwriting flt.origin here would shift the station
            # coordinates used by the Green's functions without moving the
            # patches.
            if g['method'] == 'wang':
                green = EDCMP(flt, self.data, g)
            else:
                green = Okada(flt, self.data, g)
            green.build_greens()       # load 或 generate+save，然后 rotate
            G_blocks.append(self._green_block(green.G, n_geod, nf3, i, 'GPS'))
            Gsar_blocks.append(self._green_block(green.G_sar, n_sar, nf3, i, 'InSAR'))

        self.Gmat  = np.hstack(G_blocks)
        self.G_sar = np.hstack(Gsar_blocks)
        return self.Gmat, self.G_sar

    def _green_config(self, index, fault_cfg):
        """合并全局/逐断层 Green 配置并补默认值。"""
        cfg = fault_cfg or {}
        g   = dict(self.config.get('green') or {})
        # rake_beta is per-fault only.  A value in the global green section
        # would rotate the Green's functions for every fault while the output
        # layer (which only reads per-fault configs) would write unrotated
        # slips.  Ignoring it here keeps both layers consistent.
        g.pop('rake_beta', None)
        g.update(cfg.get('green') or {})

        g.setdefault('method', 'okada')
        g.setdefault('greentype', [1, 1, 0])
        g.setdefault('nu', 0.25)
        g.setdefault('modulus', 3e10)
        g.setdefault('gps_ramp', False)
        g.setdefault('sar_ramp', False)
        g['greenfile'] = self._green_file(index, g.get('greenfile'))

        method = str(g['method']).lower()
        if method in ('okada', 'rectokada'):
            g['method'] = 'okada'
        elif method in ('wang', 'edcmp', 'layerwang'):
            g['method'] = 'wang'
            # 层状模型 Green 函数：先用 edgrn 生成，再交给 edcmp 计算形变
            g.setdefault('grndir', './')
            g.setdefault('grn_ssfile', 'edgrnhs.ss')
            g.setdefault('grn_dsfile', 'edgrnhs.ds')
            g.setdefault('grn_clfile', 'edgrnhs.cl')
            g.setdefault('edcmp_exe', 'edcmp')
        else:
            raise NotImplementedError(
                "Green's function method %r is not implemented; "
                "use 'okada' or 'wang'." % g['method'])
        return g

    def _green_file(self, index, explicit):
        """返回逐断层唯一的 Green 缓存路径或 'NONE'。

        逐断层 greenfile 优先，其次 global.greenfile：
            'NONE' -> 不缓存（不读也不写）
            'SAVE' 或空 -> outpath/green_XX.h5
            其它路径 -> {path}_XX.h5
        统一使用 .h5，因为 BaseGreen.save_greens 只写 HDF5。
        """
        globalfile = self.config.get('greenfile')
        requested = explicit
        if requested is None or str(requested).strip() == '':
            requested = globalfile
        if requested is None or str(requested).strip() == '':
            return os.path.join(self.outpath, 'green_%02d.h5' % index)

        key = str(requested).strip().upper()
        if key == 'NONE':
            return 'NONE'
        if key == 'SAVE':
            return os.path.join(self.outpath, 'green_%02d.h5' % index)

        base = str(requested)
        root, ext = os.path.splitext(base)
        if ext.lower() != '.h5':
            if ext:
                logging.warning(
                    'Green cache %s is not HDF5; using .h5 instead.', base)
            ext = '.h5'
        return '%s_%02d%s' % (root, index, ext)

    @staticmethod
    def _green_block(M, nrows, ncols, index, label):
        """校验一个断层的 Green 块形状；对应数据为空时返回 (0, 3*nf)。"""
        M = np.asarray(M)
        if M.size == 0:
            if nrows == 0:
                return np.zeros((0, ncols))
            raise ValueError(
                "Green's function block %d (%s) is empty but %d observations "
                "are expected; delete the cached Green file or change "
                "greenfile." % (index, label, nrows))
        if M.shape != (nrows, ncols):
            raise ValueError(
                "Green's function block %d (%s) has shape %s, expected %s; "
                "delete the cached Green file or change greenfile."
                % (index, label, M.shape, (nrows, ncols)))
        return M

    def build_ramp(self):
        """逐文件/逐 track 构建 ramp 设计矩阵，返回 (G_ramp, n_ramp)。
        行序 = [GPS(ndim), 水准(1), InSAR(1)]，未启用 GPS ramp 的文件用 0 列块占位。

        InSAR 的 LOS 绝对基准未知，因此每个有数据的 track 至少估计一个
        常数（参考电平）；ramp=True 时再额外估计东西/南北两个坡度项。
        只有确实需要坐标时才做投影，且使用单次构造投影的批量转换；
        否则全部 InSAR 点逐个新建 pyproj.Proj 会非常浪费。
        """
        data   = self.data
        origin = self.origin
        ndim   = data.ndim

        gps_entries = list(data.gps_pairs())
        sar_entries = list(data.sar_pairs())

        gps_need = any(bool((g or {}).get('ramp', False)) and
                       len(gps.llh_gps) > 0 for gps, g in gps_entries)
        # 每个非空 InSAR track 都需要常数项，所以只要有数据就要投影坐标。
        sar_need = any(len(sar.llh_sar) > 0 for sar, _ in sar_entries)

        xy_gps = (gt.llh2utm_batch(data.llh_gps, origin)
                  if gps_need else None)
        xy_sar = (gt.llh2utm_batch(data.llh_sar, origin)
                  if sar_need else None)

        blocks = []
        n_ramp = 0

        # GPS：每个文件一个块，未启用 ramp 时用 0 列块占位
        sta0 = 0
        for gps, g in gps_entries:
            g = g or {}
            n = len(gps.llh_gps)
            if g.get('ramp', False) and n > 0:
                helper = BaseGreen(None, gps, {'greenfile': 'NONE'})
                M = helper.make_green_gps_ramp(xy_gps[sta0:sta0 + n], ndim,
                                               g.get('ramp_method', 1))
                blocks.append(M)
                n_ramp += M.shape[1]
            else:
                blocks.append(np.zeros((n * ndim, 0)))
            sta0 += n

        # 水准：不建 ramp，只占行
        blocks.append(np.zeros((len(data.d_lev), 0)))

        # InSAR：每个有数据的 track 至少估计一个参考电平（常数项）；
        # ramp=True 时再加上 [1, east, north] 三列。
        sta0 = 0
        sar_ramp_flag = []
        for sar, s in sar_entries:
            s = s or {}
            n = len(sar.llh_sar)
            if n > 0:
                M  = BaseGreen.make_green_sar_ramp(xy_sar[sta0:sta0 + n])
                if not s.get('ramp', False):
                    M = M[:, :1]          # 只保留常数（参考电平）
                blocks.append(M)
                n_ramp += M.shape[1]
                sar_ramp_flag.append(bool(s.get('ramp', False)))
            else:
                blocks.append(np.zeros((0, 0)))
                sar_ramp_flag.append(False)
            sta0 += n
        self.sar_ramp_flag = np.asarray(sar_ramp_flag, dtype=bool)

        if n_ramp == 0:
            return np.empty((0, 0)), 0
        return linalg.block_diag(*blocks), n_ramp
