#!/usr/bin/env python
"""PlotlyReport -- interactive HTML report (zoom / pan / rotate 3D).

Static PNG figures in ``report.html`` cannot be zoomed and the 3D fault view
cannot be rotated.  This module renders the same content with Plotly:

    * 3D slip distribution as an interactive Mesh3d (drag to rotate,
      scroll to zoom, hover shows per-patch slip / depth / strike / dip);
    * GPS observation vs prediction vs residual as vector segments;
    * InSAR observation / prediction / residual per track as WebGL scatters
      with a shared RdBu_r scale (zero pinned to white by ``cmid=0``).

Output is a single HTML file; pass ``embed=True`` to inline plotly.js so it
opens offline.

Usage::

    from pyginvc.Export.PlotlyReport import make_report
    make_report('GPS_WX_T070A/solutions.h5')
"""

import logging
import os

import numpy as np

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from pyginvc.Export import View

logger = logging.getLogger(__name__)

KM_PER_DEG = 111.19


def _diverge_colorscale(cmin, cmax, name='RdBu_r'):
    """构造 diverging 色标，使 v=0 对应白色（色标中点）。

    颜色位置按 cmin..0 和 0..cmax 两段分别拉伸，正负不对称时
    仍保证零值落在白色。
    """
    from plotly.colors import sample_colorscale

    p = np.linspace(0.0, 1.0, 257)
    vals = cmin + p * (cmax - cmin)
    # v<0: cmin -> 0.0, 0 -> 0.5；v>0: 0 -> 0.5, cmax -> 1.0
    pos = np.where(vals <= 0.0,
                   0.5 * (vals - cmin) / (-cmin) if cmin < 0 else 0.0,
                   0.5 + 0.5 * vals / cmax if cmax > 0 else 0.5)
    pos = np.clip(pos, 0.0, 1.0)
    cols = sample_colorscale(name, pos)
    return [[float(t), c] for t, c in zip(p, cols)]


def _fmt(v, nd=2):
    try:
        return ('%.' + str(nd) + 'f') % float(v)
    except (TypeError, ValueError):
        return str(v)


def _slip_stats(felem):
    """总滑动、矩震级等摘要（length/width 为 km，slip 为 mm）。"""
    slot = np.asarray(felem['total_slip'], float)
    length = np.asarray(felem['length'], float)
    width = np.asarray(felem['width'], float)
    area = length * width * 1e6                    # km^2 -> m^2
    slip_m = slot / 1000.0                         # mm -> m
    mu = 3.0e10
    mo = float(np.sum(mu * slip_m * area))
    mw = (2.0 / 3.0) * (np.log10(mo) - 9.1) if mo > 0 else float('nan')
    return {
        'max_slip_mm': float(slot.max()) if slot.size else 0.0,
        'mean_slip_mm': float(slot.mean()) if slot.size else 0.0,
        'n_patch': int(slot.size),
        'moment': mo,
        'mw': mw,
    }


def load_fault_felem(outdir):
    """按 faultgeom2geojson.py 的方式从 slipmodel 的 FaultGeom 文件构造
    带滑移的 patch 顶点表（DataFrame，FaultGeom2AllVertex 输出格式）。

    优先使用逐断层文件 slipmodel_*_seg*.faultgeom；没有时退回到合并文件。
    """
    import glob

    from pyginvc.Geometry.Patch import Fault

    files = sorted(glob.glob(os.path.join(outdir, 'slipmodel_*_seg*.faultgeom')))
    if not files:
        files = sorted(glob.glob(os.path.join(outdir, 'slipmodel_*.faultgeom')))
    if not files:
        return None

    frames = []
    for fname in files:
        geom = np.atleast_2d(np.genfromtxt(fname, comments='#'))
        flt = Fault(fname, len(geom), 1, False)
        flt.load_fault()
        frames.append(flt.FaultGeom2AllVertex())
    if len(frames) == 1:
        return frames[0]
    import pandas as pd
    return pd.concat(frames, ignore_index=True)


def figure_slip_3d(felem, depth_exaggeration=1.0, show_arrows=True,
                   slip_scale=0.6):
    """按 Fault.plot_faultgeom 的构造绘制 3D 滑动分布。

    每个子断层是独立的闭合四边形：四角点 (lb, rb, rt, lt) 取自
    FaultGeom2AllVertex，顶点与三角面索引均按 patch 分组（patch p 占
    顶点 4p..4p+3），绝不跨 patch 连线。关闭光照使四边形整块同色，
    不出现对角三角着色；逐块纯色（jet）对应 plot_faultgeom 的 facecolors。
    坐标为经纬度 + 深度（深度向下为正，z 轴反向显示）。
    """
    # patch p 的四角点: lb, rb, rt, lt（与 plot_faultgeom 的顺序一致）
    corners = np.stack([
        np.column_stack([felem['lb_lon'], felem['lb_lat'], -felem['lb_dep']]),
        np.column_stack([felem['rb_lon'], felem['rb_lat'], -felem['rb_dep']]),
        np.column_stack([felem['rt_lon'], felem['rt_lat'], -felem['rt_dep']]),
        np.column_stack([felem['lt_lon'], felem['lt_lat'], -felem['lt_dep']]),
    ], axis=1).astype(float)                       # (n, 4, 3)
    corners[:, :, 2] *= depth_exaggeration

    v = corners.reshape(-1, 3)                     # patch p: 4p .. 4p+3
    n = len(felem)
    base = 4 * np.arange(n)
    # 每块两个三角形 (lb,rb,rt) 与 (lb,rt,lt)，索引只在本 patch 内
    i = np.concatenate([base, base])
    j = np.concatenate([base + 1, base + 2])
    k = np.concatenate([base + 2, base + 3])

    slot = np.nan_to_num(np.asarray(felem['total_slip'], float))
    intensity = np.repeat(slot, 4)                 # 每块四顶点同色
    hover = np.repeat(
        np.array(['ss %.0f mm  ds %.0f mm<br>total %.0f mm<br>'
                  'dip %.1f&deg;  strike %.1f&deg;<br>shallow edge depth %.2f km'
                  % (a, b, s, d, st, td) for a, b, s, d, st, td in zip(
                      np.nan_to_num(np.asarray(felem['strike_slip'], float)),
                      np.nan_to_num(np.asarray(felem['dip_slip'], float)),
                      slot,
                      np.asarray(felem['dip'], float) % 180.0,
                      np.asarray(felem['strike'], float) % 360.0,
                      -np.asarray(felem['lb_dep'], float))]),
        4)

    smax = float(slot.max()) if slot.size else 1.0
    fig = go.Figure(go.Mesh3d(
        x=v[:, 0], y=v[:, 1], z=v[:, 2], i=i, j=j, k=k,
        intensity=intensity, colorscale='Jet',
        cmin=0.0, cmax=smax,
        intensitymode='vertex',
        text=hover, hovertemplate='%{text}<extra></extra>',
        colorbar=dict(title='Total slip (mm)', len=0.7),
        # 纯 ambient：非平面四边形的两半不会被光照染成不同颜色
        lighting=dict(ambient=1.0, diffuse=0.0, specular=0.0),
        flatshading=False,
    ))

    # 每个 patch 的闭合边界（lb->rb->rt->lt->lb），无对角线
    xs, ys, zs = [], [], []
    for p in range(n):
        q = corners[p]
        xs += [q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[0, 0], None]
        ys += [q[0, 1], q[1, 1], q[2, 1], q[3, 1], q[0, 1], None]
        zs += [q[0, 2], q[1, 2], q[2, 2], q[3, 2], q[0, 2], None]
    fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                               line=dict(color='rgba(30,30,30,0.45)', width=1),
                               hoverinfo='skip', showlegend=False))

    # patch 平面内的滑动矢量箭头（与 plot_faultgeom 相同构造）
    if show_arrows:
        segments = View._slip_arrow_segments(felem, 'llh', slip_scale)
        if segments:
            ax_, ay_, az_ = [], [], []
            for p0, p1 in segments:
                ax_ += [p0[0], p1[0], None]
                ay_ += [p0[1], p1[1], None]
                az_ += [-p0[2] * depth_exaggeration,
                        -p1[2] * depth_exaggeration, None]
            fig.add_trace(go.Scatter3d(x=ax_, y=ay_, z=az_, mode='lines',
                                       line=dict(color='rgba(20,20,20,0.8)',
                                                 width=3),
                                       hoverinfo='skip', showlegend=False))

    fig.update_layout(
        scene=dict(
            xaxis_title='Longitude', yaxis_title='Latitude',
            zaxis_title='Depth (km, x%.1f)' % depth_exaggeration,
            zaxis=dict(autorange='reversed'),
            aspectmode='manual',
            aspectratio=dict(x=1.0, y=1.0, z=0.55),
            camera=dict(eye=dict(x=-1.5, y=-1.5, z=0.9)),
        ),
        margin=dict(l=0, r=0, t=10, b=0),
        height=650,
    )
    return fig


def _split_fault_blocks(felem):
    """把合并的 patch 表按断层拆开。

    FaultGeom 文件的 patch 按「倾角向为主序」排列，块内深度逐行加深，
    新断层块从最浅一行重新开始，因此以深度回升作为分界。
    """
    dep = np.round(np.asarray(felem['lb_dep'], float), 3)
    cuts = [0] + list(np.where(np.diff(dep) > 1e-6)[0] + 1) + [len(dep)]
    return [felem.iloc[a:b].reset_index(drop=True)
            for a, b in zip(cuts[:-1], cuts[1:])]


def _slip_grid_2d(block):
    """把一个断层块投影到 (沿走向, 倾角向) 平面网格，参考 plot_slip_2d.py。

    返回 (x_centers, y_centers, slip_grid)，失败返回 None。
    """
    n = len(block)
    dep = np.round(np.asarray(block['lb_dep'], float), 3)
    ndip = len(np.unique(dep))
    if ndip < 2 or n % ndip:
        return None
    nstrike = n // ndip

    # patch 中心（浅边中点）
    lat_c = 0.5 * (np.asarray(block['lb_lat'], float) +
                   np.asarray(block['rb_lat'], float))
    lon_c = 0.5 * (np.asarray(block['lb_lon'], float) +
                   np.asarray(block['rb_lon'], float))
    lat_ref, lon_ref = lat_c.mean(), lon_c.mean()
    coslat = np.cos(np.radians(lat_ref))
    x_e = (lon_c - lon_ref) * KM_PER_DEG * coslat
    x_n = (lat_c - lat_ref) * KM_PER_DEG

    # 最浅一层 patch 中心的主方向作为走向
    shallow = dep == dep.min()
    pts = np.vstack([x_e[shallow], x_n[shallow]]).T
    pts = pts - pts.mean(axis=0)
    _, _, vh = np.linalg.svd(pts, full_matrices=False)
    strike_vec = vh[0]

    strike_pos = (np.vstack([x_e, x_n]).T @ strike_vec).reshape(ndip, nstrike)
    order = np.argsort(strike_pos[0])
    strike_pos = strike_pos[:, order]

    width = np.asarray(block['width'], float).reshape(ndip, nstrike)[:, order]
    slip = np.asarray(block['total_slip'], float).reshape(ndip, nstrike)[:, order]

    p = strike_pos[0]
    mid = 0.5 * (p[:-1] + p[1:])
    left = p[0] - 0.5 * (p[1] - p[0])
    right = p[-1] + 0.5 * (p[-1] - p[-2])
    x_edges = np.concatenate([[left], mid, [right]])
    x_edges -= x_edges[0]
    y_edges = np.concatenate([[0.0], np.cumsum(width[:, 0])])

    x_c = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_c = 0.5 * (y_edges[:-1] + y_edges[1:])
    return x_c, y_c, slip


def figure_slip_2d(felem):
    """断层面上的 2D 滑动分布：横轴沿走向长度，纵轴倾角向宽度。

    参考 test/plot_slip_2d.py：按断层块分别建网格，go.Heatmap 渲染
    （pcolormesh 的 plotly 等价物），共享 jet 色标、0 起。
    """
    blocks = _split_fault_blocks(felem)
    grids = [(k, _slip_grid_2d(b)) for k, b in enumerate(blocks)]
    grids = [(k, g) for k, g in grids if g is not None]
    if not grids:
        return None

    smax = max(float(g[2].max()) for _, g in grids)
    fig = make_subplots(
        rows=1, cols=len(grids), shared_yaxes=False,
        subplot_titles=['Segment %d' % (k + 1) for k, _ in grids])
    for m, (k, (x_c, y_c, slip)) in enumerate(grids, start=1):
        fig.add_trace(go.Heatmap(
            x=x_c, y=y_c, z=slip,
            colorscale='Jet', zmin=0.0, zmax=smax,
            showscale=(m == 1),
            colorbar=dict(title='Slip (mm)', len=0.75),
            hovertemplate='along-strike %{x:.1f} km<br>'
                          'down-dip %{y:.1f} km<br>'
                          'slip %{z:.0f} mm<extra></extra>',
        ), row=1, col=m)
        fig.update_xaxes(title_text='Along-strike (km)', row=1, col=m)
        fig.update_yaxes(title_text='Down-dip (km)', autorange='reversed',
                         row=1, col=m)
        # 等比例，便于看长宽
        fig.update_yaxes(scaleanchor='x' if m == 1 else 'x%d' % m,
                         scaleratio=1.0, row=1, col=m)
    fig.update_layout(height=430, margin=dict(l=10, r=10, t=50, b=10),
                      showlegend=False)
    return fig


def _arrow_scale(lon, lat, *uv, fraction=0.2):
    """返回 (sx, sy)：mm -> 经纬度偏移的换算系数。

    使最长箭头约占图幅短边的 fraction；经度按 cos(lat) 修正，
    两轴在真实距离上长度一致。
    """
    span_lon = max(np.ptp(lon), 1e-3)
    span_lat = max(np.ptp(lat), 1e-3)
    coslat = np.cos(np.radians(np.nanmean(lat)))
    span_x_km = span_lon * KM_PER_DEG * coslat
    span_y_km = span_lat * KM_PER_DEG
    amp = max([float(np.nanmax(np.abs(np.asarray(u)))) for u in uv] + [1e-6])
    target_km = fraction * min(span_x_km, span_y_km)   # 最长箭头的物理长度
    km_per_mm = target_km / amp                       # 1 mm 对应的箭头长度 (km)
    sx = km_per_mm / (KM_PER_DEG * coslat)            # deg per mm (经度)
    sy = km_per_mm / KM_PER_DEG                       # deg per mm (纬度)
    return sx, sy


def figure_gps(llh_gps, d_gps, m_gps, felem=None, arrow_scale=1.0):
    """GPS 观测/预测（线段）+ 残差；三分量时叠加竖直分量散点。

    arrow_scale : 箭头长度倍数，1.0 为默认（最长箭头约占图幅短边 20%），
                  大于 1 更长，小于 1 更短。
    """
    lat, lon = llh_gps[:, 0], llh_gps[:, 1]
    ndim = d_gps.shape[1]
    sx, sy = _arrow_scale(lon, lat, d_gps[:, 0], d_gps[:, 1],
                          m_gps[:, 0], m_gps[:, 1],
                          fraction=0.2 * float(arrow_scale))

    def segments(u, v):
        xs, ys = [], []
        for x0, y0, du, dv in zip(lon, lat, u, v):
            xs += [x0, x0 + du * sx, None]
            ys += [y0, y0 + dv * sy, None]
        return xs, ys

    fig = go.Figure()
    for name, (u, v), color in (
            ('Observation', (d_gps[:, 0], d_gps[:, 1]), '#1f4fd8'),
            ('Prediction', (m_gps[:, 0], m_gps[:, 1]), '#d62728')):
        xs, ys = segments(u, v)
        fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines',
                                 line=dict(color=color, width=2),
                                 name=name, hoverinfo='skip'))
    fig.add_trace(go.Scatter(x=lon, y=lat, mode='markers',
                             marker=dict(size=6, color='black'),
                             name='GPS stations', customdata=np.column_stack([lon, lat]),
                             hovertemplate='lon %{customdata[0]:.4f}<br>'
                                           'lat %{customdata[1]:.4f}<extra></extra>'))

    if ndim == 3:
        up_min = float(min(d_gps[:, 2].min(), m_gps[:, 2].min()))
        up_max = float(max(d_gps[:, 2].max(), m_gps[:, 2].max()))
        fig.add_trace(go.Scatter(
            x=lon, y=lat, mode='markers',
            marker=dict(size=14, color=d_gps[:, 2],
                        colorscale=_diverge_colorscale(up_min, up_max),
                        cmin=up_min, cmax=up_max,
                        line=dict(color='black', width=0.5),
                        colorbar=dict(title='Up obs (mm)'), showscale=True),
            name='Up obs', customdata=d_gps[:, 2],
            hovertemplate='Up obs %{customdata:.1f} mm<extra></extra>'))

    if felem is not None:
        fig.add_trace(_fault_outline(felem))

    fig.update_layout(
        xaxis_title='Longitude', yaxis_title='Latitude',
        yaxis=dict(scaleanchor='x', scaleratio=1.0 /
                   np.cos(np.radians(np.mean(lat)))),
        height=560, margin=dict(l=10, r=10, t=30, b=10),
        legend=dict(orientation='h', y=1.08),
        title='GPS observation (blue) vs prediction (red)',
    )
    return fig


def _fault_outline(felem):
    """断层 patch 网格：每个子断层画闭合四边形（用于 GPS 平面图）。"""
    x = np.column_stack([felem[c + '_lon'] for c in ('lb', 'rb', 'rt', 'lt')])
    y = np.column_stack([felem[c + '_lat'] for c in ('lb', 'rb', 'rt', 'lt')])
    xs, ys = [], []
    for p in range(len(felem)):
        xs += [x[p, 0], x[p, 1], x[p, 2], x[p, 3], x[p, 0], None]
        ys += [y[p, 0], y[p, 1], y[p, 2], y[p, 3], y[p, 0], None]
    return go.Scatter(x=xs, y=ys, mode='lines',
                      line=dict(color='rgba(60,60,60,0.45)', width=1),
                      name='Fault patches', hoverinfo='skip')


def figure_insar(llh_sar, counts, ob_all, mb_all, label='InSAR',
                 felem=None, slip_overlay=True):
    """每个 track 一组 观测/预测/残差；三者共享 RdBu_r、零值白。

    可选在残差图上叠加断层滑动分布（逐 patch 闭合多边形填充，jet 配色）。
    叠加图层带 meta='slip_overlay'，页面上的开关通过 Plotly.restyle 显隐。
    """
    from plotly.colors import sample_colorscale

    overlay_traces = []            # 用于页面上加开关
    if felem is not None and slip_overlay:
        slot = np.nan_to_num(np.asarray(felem['total_slip'], float))
        smax = float(slot.max()) if slot.size else 1.0
        xs = np.column_stack([felem[c + '_lon'] for c in ('lb', 'rb', 'rt', 'lt')])
        ys = np.column_stack([felem[c + '_lat'] for c in ('lb', 'rb', 'rt', 'lt')])
        norm = slot / smax if smax > 0 else np.zeros_like(slot)
        rgba = []
        for rgb in sample_colorscale('Jet', norm):
            r, g, b = [int(float(v)) for v in rgb[4:-1].split(',')]
            rgba.append('rgba(%d,%d,%d,0.72)' % (r, g, b))
        overlay_traces = list(zip(xs, ys, rgba))

    start = 0
    figs = []
    for k, count in enumerate(counts):
        end = start + count
        if count > 0:
            ob, mb = ob_all[start:end], mb_all[start:end]
            rb = ob - mb
            lon = llh_sar[start:end, 1]
            lat = llh_sar[start:end, 0]
            vmin = min(ob.min(), mb.min(), rb.min())
            vmax = max(ob.max(), mb.max(), rb.max())
            if vmax <= vmin:
                vmin, vmax = vmin - 1.0, vmax + 1.0

            fig = make_subplots(
                rows=1, cols=3, shared_yaxes=True,
                subplot_titles=('%s%d Observation' % (label, k + 1),
                                '%s%d Prediction' % (label, k + 1),
                                '%s%d Residual' % (label, k + 1)))
            for col, (val, name) in enumerate(
                    ((ob, 'obs'), (mb, 'pred'), (rb, 'res')), start=1):
                marker = dict(size=3, color=val,
                              colorscale=_diverge_colorscale(vmin, vmax),
                              cmin=vmin, cmax=vmax,
                              showscale=(col == 3))
                if col == 3:
                    marker['colorbar'] = dict(title='LOS (mm)', len=0.45, y=0.78)
                fig.add_trace(go.Scattergl(
                    x=lon, y=lat, mode='markers',
                    marker=marker, name=name,
                    hovertemplate='lon %{x:.4f}<br>lat %{y:.4f}<br>'
                                  'LOS %{marker.color:.1f} mm<extra></extra>'
                ), row=1, col=col)

            # 残差图上叠加断层滑动分布（可开关）
            if overlay_traces:
                for qx, qy, fc in overlay_traces:
                    fig.add_trace(go.Scatter(
                        x=list(qx), y=list(qy), fill='toself', mode='lines',
                        line=dict(color='rgba(0,0,0,0.3)', width=0.3),
                        fillcolor=fc, hoverinfo='skip', showlegend=False,
                        meta='slip_overlay'), row=1, col=3)
                # 滑动色标（不可见的散点只为生成 colorbar）
                fig.add_trace(go.Scattergl(
                    x=[None], y=[None], mode='markers',
                    marker=dict(color=[0.0], colorscale='Jet', cmin=0.0,
                                cmax=float(np.nan_to_num(np.asarray(
                                    felem['total_slip'], float)).max()),
                                showscale=True,
                                colorbar=dict(title='Slip (mm)', len=0.45,
                                              y=0.22)),
                    hoverinfo='skip', showlegend=False, meta='slip_overlay'),
                    row=1, col=3)

            fig.update_yaxes(scaleanchor='x2', scaleratio=1.0 /
                             np.cos(np.radians(np.mean(lat))), row=1, col=1)
            fig.update_layout(height=420, margin=dict(l=10, r=10, t=40, b=10),
                              showlegend=False,
                              title='%s track %d  (%d pts, shared scale %s..%s mm)'
                              % (label, k + 1, count, _fmt(vmin, 1), _fmt(vmax, 1)))
            figs.append(fig)
        start = end
    return figs


def make_report(solutionfile, outfile=None, embed=True,
                depth_exaggeration=1.0, title=None, arrow_scale=1.0):
    """从 solutions.h5 生成交互式 HTML 报告，返回文件路径。

    arrow_scale : GPS 箭头长度倍数（1.0 为默认，越大箭头越长）。
    """
    import h5py

    outdir = os.path.dirname(os.path.abspath(solutionfile))
    outfile = outfile or os.path.join(outdir, 'report_interactive.html')
    title = title or ('PyGINVC interactive report - %s'
                      % os.path.basename(outdir))

    with h5py.File(solutionfile, 'r') as h5:
        llh_gps = h5['obs/llh_gps'][()]
        ndim = int(h5['obs/ndim'][()])
        d_gps = h5['obs/d_gps'][()].reshape(len(llh_gps), ndim)
        m_gps = h5['sol/dhat'][0:ndim * len(llh_gps)].reshape(len(llh_gps), ndim)
        llh_sar = h5['obs/llh_sar'][()]
        n_lev = len(h5['obs/llh_lev'][()]) if 'llh_lev' in h5['obs'] else 0
        offset = ndim * len(llh_gps) + n_lev
        d_sar = h5['obs/d_sar'][()]
        m_sar = h5['sol/dhat'][offset:offset + len(d_sar)]
        counts = View._sar_counts(h5)
        # 与 View.plot_slip_3d 一致，先从 h5 构造 patch；失败再读 faultgeom 文件
        try:
            felem = View._fault_felem(h5)
        except Exception as exc:
            logger.warning('View._fault_felem failed: %s', exc)
            felem = None
        if felem is None:
            try:
                felem = load_fault_felem(outdir)
            except Exception:
                felem = None

    sections = []
    stats = _slip_stats(felem) if felem is not None else None
    _js_done = {'done': False}

    def add_fig(fig, sid):
        # plotly.js 只内联一次（第一个图），后续图复用
        include = (not _js_done['done'])
        _js_done['done'] = True
        div = pio.to_html(
            fig, full_html=False,
            include_plotlyjs=(('inline' if embed else 'cdn') if include else False),
            config={'displaylogo': False,
                    'modeBarButtonsToRemove': ['lasso2d', 'select2d']},
            div_id=sid)
        sections.append('<div class="fig">%s</div>' % div)

    header = ['<h1>%s</h1>' % title]
    if stats:
        header.append(
            '<div class="cards">'
            '<div class="card"><b>%d</b><span>fault patches</span></div>'
            '<div class="card"><b>%.0f mm</b><span>max slip</span></div>'
            '<div class="card"><b>%.0f mm</b><span>mean slip</span></div>'
            '<div class="card"><b>%.2e</b><span>Mo (N m)</span></div>'
            '<div class="card"><b>%.2f</b><span>Mw</span></div>'
            '</div>' % (stats['n_patch'], stats['max_slip_mm'],
                        stats['mean_slip_mm'], stats['moment'], stats['mw']))

    if felem is not None:
        sections.append('<h2>3D slip distribution</h2>'
                        '<p class="hint">Drag to rotate, scroll to zoom, '
                        'hover for per-patch values.</p>')
        add_fig(figure_slip_3d(felem, depth_exaggeration), 'slip3d')

        fig2d = figure_slip_2d(felem)
        if fig2d is not None:
            sections.append('<h2>2D slip distribution</h2>'
                            '<p class="hint">Slip on the fault plane: x = '
                            'along-strike length, y = down-dip width '
                            '(surface at top). Hover for values.</p>')
            add_fig(fig2d, 'slip2d')

    if len(llh_gps):
        sections.append('<h2>GPS observation vs prediction</h2>')
        add_fig(figure_gps(llh_gps, d_gps, m_gps, felem,
                           arrow_scale=arrow_scale), 'gps')

    if len(llh_sar):
        sections.append(
            '<h2>InSAR observation / prediction / residual</h2>'
            '<p class="hint">Shared RdBu_r scale, zero at white. '
            'Drag a box or scroll to zoom into a track.</p>'
            '<div class="controls"><label>'
            '<input type="checkbox" id="slipOverlayToggle" checked> '
            '在残差图上叠加断层滑动分布</label></div>')
        for k, fig in enumerate(figure_insar(llh_sar, counts, d_sar, m_sar,
                                             felem=felem)):
            add_fig(fig, 'sar%d' % (k + 1))

    downloads = []
    for name in sorted(os.listdir(outdir)):
        if name.endswith(('.h5', '.faultgeom', '.gmtlin', '.txt', '.summary')) \
                or name == 'summary':
            downloads.append('<li><a href="%s" download>%s</a></li>'
                             % (name, name))
    if downloads:
        sections.append('<h2>Downloads</h2><ul class="dl">%s</ul>'
                        % ''.join(downloads))

    html = """<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s</title>
<style>
  :root{--bg:#f6f7f9;--card:#fff;--line:#e3e6ea;--fg:#20242b;--muted:#667085}
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--fg);
       font-family:-apple-system,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif}
  main{max-width:1400px;margin:0 auto;padding:24px}
  h1{font-size:22px;margin:8px 0 18px}
  h2{font-size:17px;margin:28px 0 8px}
  .hint{color:var(--muted);font-size:13px;margin:0 0 8px}
  .cards{display:flex;flex-wrap:wrap;gap:12px;margin:10px 0 4px}
  .card{background:var(--card);border:1px solid var(--line);border-radius:8px;
        padding:10px 16px;min-width:130px}
  .card b{display:block;font-size:20px}
  .card span{color:var(--muted);font-size:12px}
  .fig{background:var(--card);border:1px solid var(--line);border-radius:8px;
       padding:8px;margin:10px 0 18px;overflow:hidden}
  ul.dl{columns:2;font-size:13px}
  ul.dl a{color:#1f4fd8;text-decoration:none}
  .controls{display:inline-block;margin:2px 0 10px;padding:8px 12px;
            background:var(--card);border:1px solid var(--line);border-radius:8px;
            font-size:14px;cursor:pointer}
  .controls input{vertical-align:-2px;margin-right:4px}
</style></head><body><main>
%s
</main>
<script>
(function () {
  var cb = document.getElementById('slipOverlayToggle');
  if (!cb) return;
  function apply() {
    var on = cb.checked;
    document.querySelectorAll('.js-plotly-plot').forEach(function (gd) {
      if (!gd.data) return;
      var idx = [];
      gd.data.forEach(function (t, i) { if (t && t.meta === 'slip_overlay') idx.push(i); });
      if (idx.length && window.Plotly) {
        Plotly.restyle(gd, {visible: on}, idx);
      }
    });
  }
  cb.addEventListener('change', apply);
  window.addEventListener('load', apply);
})();
</script>
</body></html>""" % (title, '\n'.join(header + sections))

    with open(outfile, 'w', encoding='utf-8') as fid:
        fid.write(html)
    logger.info('Interactive report written to %s', outfile)
    return outfile


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Interactive plotly report')
    parser.add_argument('solutionfile', help='path to solutions.h5')
    parser.add_argument('-o', '--outfile', default=None)
    parser.add_argument('--no-embed', action='store_true',
                        help='use CDN plotly.js instead of inlining it')
    parser.add_argument('--depth-exaggeration', type=float, default=1.0)
    parser.add_argument('--arrow-scale', type=float, default=1.0,
                        help='GPS 箭头长度倍数（默认 1.0；2=加倍，0.5=减半）')
    args = parser.parse_args()
    print(make_report(args.solutionfile, args.outfile,
                      embed=not args.no_embed,
                      depth_exaggeration=args.depth_exaggeration,
                      arrow_scale=args.arrow_scale))
