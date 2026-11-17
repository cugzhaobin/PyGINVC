#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build a non-crossing quad fault mesh from a curved surface trace.

Pipeline
--------
1. Extrude the surface trace down dip exactly like trace2faultgeom.py does
   (straight rays along the local surface normal).  This is the "preliminary"
   geometry; on strongly curved traces its columns cross at depth.
2. Extract the iso-depth contours (one polyline per depth level).
3. Smooth along strike and resample: the default ``direction`` method smooths
   the down-dip direction field, so the contours keep the along-strike shape
   of the trace while the columns stop converging; ``nodes`` smooths the
   contour coordinates directly (flattens bends more strongly).  The window
   is increased until every quader is non-degenerate.
4. Split the smoothed contours into quadrilateral subfaults.

Outputs
-------
``<base>_initial.faultgeom``   preliminary extrusion (PyGINV format)
``<base>_contours.dat``        contours: level depth_km lon lat
``<base>_quads.dat``           exact quad corners (+ strike/dip/width)
``<base>_smoothed.faultgeom``  smoothed mesh in PyGINV faultgeom format
``<base>_mesh.png/.pdf``       QC figure
``<base>_profile.png/.pdf``    depth versus horizontal distance cross-section
                               (true scale: 1 km is the same on both axes)
``<base>_vertex.dat``/``<base>_element.dat``
                               triangular mesh for TriSlipInversion.py
                               (written with ``--tri``; every quad is split
                               into two triangles)

Note: PyGINV faultgeom encodes each patch as the shallow edge plus a
perpendicular width/dip offset, so a general quad is written as its
best-fit parallelogram.  The script reports the resulting encoding error
(checked by reloading the file with pyginvc.Geometry.Patch.Fault).

All geometry arguments are required; meshing details (surface smoothing,
contour smoothing, strike mode, dip profile) have built-in sensible
defaults.  ``--strike-mode segment`` extrudes every patch along its own
segment normal (the iso-depth contours may cross on a strongly curved
trace).  ``--strike-mode mean`` (default) keeps the exact surface trace
shape at every depth: each iso-depth contour is a translated copy, shifted
along the normal to the mean strike by the cumulative ``dz/tan(dip)`` from
the dip profile (a pure vertical shift when the dip is 90 deg), so no
contours can cross and the deep fault has the same curved shape as the
surface.  Dip profiles: ``linear``, ``exp`` / ``log`` (steep at the top,
progressively gentler down dip; ``--profile-decay``), ``power``
(f=t**--power: p<1 fast at top, p>1 fast at depth), and ``logistic``
(sharp turnover at ``--profile-t0`` with width ``--profile-s``), plus
``sine`` / ``cosine`` / ``smoothstep``.

Usage
-----
    python smooth_fault_mesh.py --trace GaoFen_Fault7.txt --nsegs 23 \
        --ndeps 10 --top-depth 0 --bottom-depth 18 \
        --top-dip 60 --bottom-dip 40 --profile exp \
        --dip-azimuth 250 --output fault_mesh

Optional: ``--strike-mode {mean,segment}`` (default mean), ``--profile
{linear,exp,power,log,logistic,sine,cosine,smoothstep}`` (default linear),
``--power`` / ``--profile-t0`` / ``--profile-s``, ``--smooth`` (trace smoothing window,
default 10), ``--contour-method {direction,nodes}``, ``--contour-smooth-km``,
``--tri`` (also write TriSlipInversion.py vertex/element files).  Sizing
options for an independently re-meshed triangular grid are described below.

Triangular mesh
---------------
With ``--tri`` the smoothed quad grid is triangulated (two triangles per
quad, split along the shorter diagonal) and written in the format read by
``pyginvc.Geometry.Triangle.load_fault``:

    <base>_vertex.dat    idx lat lon depth_km   (depth positive downward)
    <base>_element.dat   idx1 idx2 idx3 ss ds op (mm; initial slip = 0)

For a mesh that is *not* tied to the quad grid, pass a target edge length
(km).  The fault surface is then re-sampled on staggered depth rows: alternate
rows are shifted by half a node and consecutive rows are ``L*sqrt(3)/2`` apart
down dip, so the Delaunay triangulation in the (along-strike, down-dip)
parameter plane consists of near-equilateral triangles with edge length ~L:

    --tri-edge 1.5                     uniform target edge length
    --tri-edge-top 1.0 --tri-edge-bottom 3.0
                                       finer shallow, coarser deep
    --tri-edge 2.0 --tri-edge-bottom 4.0
                                       --tri-edge is the default for an
                                       unspecified top/bottom value

``--tri-grading {geometric,linear,smoothstep}`` controls how the size
transitions from the shallow to the deep edge (default geometric, i.e. a
constant ratio between neighbouring rows).  Sizing options imply ``--tri``.
Graded meshes also write ``<base>_tri_mesh.png/.pdf`` (map view plus edge
length versus depth) for QC.

Point ``dict_fault.vertexfile``/``elementfile`` in the TriSlipInversion
configuration at these two files.
"""

import argparse
import os
import sys

import numpy as np

try:
    from pyproj import Proj
except ImportError:                                    # pragma: no cover
    Proj = None


# ---------------------------------------------------------------------------
# trace input / preparation
# ---------------------------------------------------------------------------
def angular_diff(a, b):
    """Signed smallest difference a - b wrapped into [-180, 180)."""
    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0


def mean_azimuth(east, north):
    """Circular mean of azimuths (degrees clockwise from north)."""
    az = np.degrees(np.arctan2(east, north))
    z = np.mean(np.exp(1j * np.radians(az)))
    return np.degrees(np.angle(z)) % 360.0


class LocalProjection(object):
    """Local transverse-Mercator projection, same style as PyGINV geotools."""

    def __init__(self, lat0, lon0):
        if Proj is None:
            raise ImportError("pyproj is required (pip install pyproj)")
        self.lat0 = float(lat0)
        self.lon0 = float(lon0)
        self.proj = Proj(ellps='WGS84', proj='tmerc', lon_0=self.lon0)
        self.e0, self.n0 = self.proj(self.lon0, self.lat0)

    def to_xy(self, lon, lat):
        """lon/lat (deg) -> local east/north (km)."""
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        e, n = self.proj(lon, lat)
        return (np.asarray(e) - self.e0) / 1000.0, \
               (np.asarray(n) - self.n0) / 1000.0

    def to_ll(self, east, north):
        """local east/north (km) -> lon/lat (deg)."""
        east = np.asarray(east, dtype=float)
        north = np.asarray(north, dtype=float)
        lon, lat = self.proj(east * 1000.0 + self.e0,
                             north * 1000.0 + self.n0, inverse=True)
        return np.asarray(lon), np.asarray(lat)


def read_trace(filename, order='auto'):
    """Read a trace file.  Returns lat, lon, resolved column order."""
    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(1, data.size)
    if data.shape[1] < 2:
        raise ValueError("trace file needs at least two columns (lon lat)")
    c0 = data[:, 0]
    c1 = data[:, 1]

    if order == 'auto':
        if np.any(np.abs(c0) > 90.0) and np.all(np.abs(c1) <= 90.0):
            order = 'lonlat'
        elif np.any(np.abs(c1) > 90.0) and np.all(np.abs(c0) <= 90.0):
            order = 'latlon'
        else:
            order = 'lonlat'
            print("WARNING: cannot tell lon/lat column order from the values; "
                  "assuming 'lon lat'.  Use --order latlon if wrong.",
                  file=sys.stderr)

    if order == 'lonlat':
        lon, lat = c0.copy(), c1.copy()
    else:
        lat, lon = c0.copy(), c1.copy()

    if np.any(np.abs(lat) > 90.0):
        raise ValueError("latitude values outside [-90, 90]; check --order")
    if np.any(np.abs(lon) > 360.0):
        raise ValueError("longitude values outside [-360, 360]; check --order")
    lon = (lon + 180.0) % 360.0 - 180.0             # normalise to [-180,180)

    # drop exactly duplicated consecutive vertices
    keep = np.ones(len(lat), dtype=bool)
    if len(lat) > 1:
        keep[1:] = (np.abs(np.diff(lat)) > 0.0) | (np.abs(np.diff(lon)) > 0.0)
    lat, lon = lat[keep], lon[keep]
    if len(lat) < 2:
        raise ValueError("trace needs at least two distinct points")
    return lat, lon, order


def smooth_xy(x, y, window):
    """Moving-average smoothing of a polyline, endpoints preserved."""
    if window is None or window <= 1:
        return x, y
    if window % 2 == 0:
        window += 1
    half = window // 2
    xp = np.r_[np.repeat(x[0], half), x, np.repeat(x[-1], half)]
    yp = np.r_[np.repeat(y[0], half), y, np.repeat(y[-1], half)]
    kernel = np.ones(window) / float(window)
    xs = np.convolve(xp, kernel, mode='valid')
    ys = np.convolve(yp, kernel, mode='valid')
    xs[0], ys[0] = x[0], y[0]
    xs[-1], ys[-1] = x[-1], y[-1]
    return xs, ys


def resample_nodes(x, y, nsegs):
    """Uniform arc-length resampling of a polyline into nsegs segments."""
    if nsegs < 1:
        raise ValueError("need at least one segment along strike")
    seg = np.hypot(np.diff(x), np.diff(y))
    s = np.r_[0.0, np.cumsum(seg)]
    total = s[-1]
    if total <= 0.0:
        raise ValueError("trace has zero length")
    sn = np.linspace(0.0, total, nsegs + 1)
    return np.interp(sn, s, x), np.interp(sn, s, y), total


# ---------------------------------------------------------------------------
# dip profiles
# ---------------------------------------------------------------------------
def profile_shape(t, mode, power=1.0, decay=1.0, t0=0.5, s=0.1):
    """Shape f(t) running from 0 (top) to 1 (bottom)."""
    t = np.clip(np.asarray(t, dtype=float), 0.0, 1.0)
    if mode == 'linear':
        f = t
    elif mode == 'power':
        if power <= 0.0:
            raise ValueError("--power must be > 0")
        f = t ** power
    elif mode == 'exp':
        if decay <= 0.0:
            raise ValueError("--decay must be > 0")
        if decay < 1.0e-8:
            f = t
        else:
            f = (1.0 - np.exp(-decay * t)) / (1.0 - np.exp(-decay))
    elif mode == 'log':
        if decay <= 0.0:
            raise ValueError("--decay must be > 0")
        if decay < 1.0e-8:
            f = t
        else:
            f = np.log1p(decay * t) / np.log1p(decay)
    elif mode == 'logistic':
        if s <= 0.0:
            raise ValueError("--s must be > 0")
        if not (0.0 < t0 < 1.0):
            raise ValueError("--t0 must be between 0 and 1")
        h = 1.0 / (1.0 + np.exp((t - t0) / s))
        h0 = 1.0 / (1.0 + np.exp(-t0 / s))
        h1 = 1.0 / (1.0 + np.exp((1.0 - t0) / s))
        f = (h - h0) / (h1 - h0)
    elif mode == 'sine':
        f = np.sin(0.5 * np.pi * t)
    elif mode == 'cosine':
        f = 0.5 * (1.0 - np.cos(np.pi * t))
    elif mode == 'smoothstep':
        f = t * t * (3.0 - 2.0 * t)
    else:
        raise ValueError("unknown dip profile: %s" % mode)
    return np.clip(f, 0.0, 1.0)


def dip_profile(t, top_dip, bottom_dip, mode, power, decay, t0, s):
    f = profile_shape(t, mode, power, decay, t0, s)
    return top_dip + (bottom_dip - top_dip) * f


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------
def build_surface(args):
    """Smoothed, resampled surface trace in local (east, north) km.

    ``--strike-mode`` 只影响深部几何（见 ``blend_mean_strike``），
    地表迹线始终保留真实弯曲形状。
    """
    lat, lon, _ = read_trace(args.trace, 'auto')
    proj = LocalProjection(np.mean(lat), np.mean(lon))
    x, y = proj.to_xy(lon, lat)
    x, y = smooth_xy(x, y, args.smooth)
    xn, yn, total = resample_nodes(x, y, args.nsegs)

    az0 = mean_azimuth(np.diff(xn), np.diff(yn))
    d_right = abs(angular_diff((az0 + 90.0) % 360.0, args.dip_azimuth))
    d_left = abs(angular_diff((az0 - 90.0) % 360.0, args.dip_azimuth))
    reverse = d_left < d_right
    if reverse:
        xn, yn = xn[::-1].copy(), yn[::-1].copy()
        az0 = (az0 + 180.0) % 360.0

    A = np.column_stack([xn[:-1], yn[:-1]])
    B = np.column_stack([xn[1:], yn[1:]])
    seg = B - A
    tvec = seg / np.hypot(seg[:, 0], seg[:, 1])[:, None]
    uvec = np.column_stack([tvec[:, 1], -tvec[:, 0]])   # right normal
    return proj, A, B, uvec, total, az0, reverse


def build_mean_strike_mesh(A, B, z, dips, az_mean, dip_azimuth):
    """按平均走向把地表迹线整体拷贝到各深度层。

    返回 (C, cum, u)：
        C   = (nlev, nnode, 2)，每层都是地表迹线的纯平移副本——形状和
              总长度与地表完全相同，不会随深度扩散；
        cum = 各层的累计水平位移 (km)；
        u   = 平移方向（平均走向的法向，指向下倾侧）。

    水平位移 cum = cumsum(dz / tan(dip))，垂向就是深度 z 本身；
    倾角 90 deg 时水平位移为 0，即纯垂直拷贝。
    """
    C0 = np.vstack([A, B[-1][None, :]])            # 地表迹线节点
    t = np.array([np.sin(np.radians(az_mean)), np.cos(np.radians(az_mean))])
    u = np.array([t[1], -t[0]])                    # 平均走向的右法向
    az_u = np.degrees(np.arctan2(u[0], u[1])) % 360.0
    if abs(angular_diff(az_u, dip_azimuth)) > 90.0:
        u = -u
    cum = np.r_[0.0, np.cumsum(np.diff(z) / np.tan(np.radians(dips)))]
    C = C0[None, :, :] + cum[:, None, None] * u[None, None, :]
    return C, cum, u


def depth_levels(args):
    """Depth nodes and the dip of every row (same convention as trace2faultgeom)."""
    z = np.linspace(args.top_depth, args.bottom_depth, args.ndeps + 1)
    t = (0.5 * (z[:-1] + z[1:]) - args.top_depth) / (args.bottom_depth -
                                                     args.top_depth)
    dips = dip_profile(t, args.top_dip, args.bottom_dip, args.profile,
                       args.power, args.profile_decay,
                       args.profile_t0, args.profile_s)
    return z, dips


def initial_contours(A, B, uvec, z, dips):
    """Straight-ray extrusion: contour i is the top edge of row i."""
    dz = np.diff(z)
    hstep = dz / np.tan(np.radians(dips))
    cum = np.r_[0.0, np.cumsum(hstep)]
    n = len(A)
    C = np.empty((len(cum), n + 1, 2))
    for i, c in enumerate(cum):
        C[i, :n] = A + c * uvec
        C[i, n] = B[-1] + c * uvec[-1]
    return C, cum


def cell_area(C):
    """Plan-view signed area of every quad (rows = depth, cols = strike)."""
    A = C[:-1, :-1]
    B = C[:-1, 1:]
    E = C[1:, 1:]
    D = C[1:, :-1]

    def cross(p, q):
        return p[..., 0] * q[..., 1] - p[..., 1] * q[..., 0]

    return 0.5 * (cross(A, B) + cross(B, E) + cross(E, D) + cross(D, A))


def cell_validity(C, rel_area=0.02):
    area = cell_area(C)
    ref = np.median(np.abs(area[0]))
    sign = np.sign(np.median(area))
    bad = (area * sign <= rel_area * ref)
    return area, bad


def smooth_levels(C, window):
    """Smooth contours 1..end along strike and resample them evenly."""
    Cs = C.copy()
    for i in range(1, C.shape[0]):
        xs, ys = smooth_xy(C[i, :, 0], C[i, :, 1], window)
        xs, ys, _ = resample_nodes(xs, ys, C.shape[1] - 1)
        Cs[i, :, 0] = xs
        Cs[i, :, 1] = ys
    return Cs


def smooth_vectors(v, window):
    """Moving average of a unit-vector field along strike, renormalised."""
    if window <= 1:
        return v
    if window % 2 == 0:
        window += 1
    half = window // 2
    vp = np.vstack([np.repeat(v[0][None, :], half, 0), v,
                    np.repeat(v[-1][None, :], half, 0)])
    k = np.ones(window) / float(window)
    vs = np.column_stack([np.convolve(vp[:, 0], k, 'valid'),
                          np.convolve(vp[:, 1], k, 'valid')])
    return vs / np.hypot(vs[:, 0], vs[:, 1])[:, None]


def direction_contours(nodes, uvec, cum, window):
    """Contours from a smoothed down-dip direction field.

    Every column is a straight ray, but the ray direction is the *smoothed*
    surface normal.  Compared with smoothing the contour coordinates this
    keeps the along-strike shape (bends are not flattened) while still
    preventing the columns from crossing.
    """
    un = np.vstack([uvec[0], 0.5 * (uvec[:-1] + uvec[1:]), uvec[-1]])
    un /= np.hypot(un[:, 0], un[:, 1])[:, None]
    un = smooth_vectors(un, window)
    C = np.empty((len(cum), len(nodes), 2))
    for i, c in enumerate(cum):
        C[i] = nodes + c * un
    return C


def choose_smoothing(C, w0, wmax):
    """Window giving the fewest degenerate cells (ties -> smaller window).

    Strong along-strike smoothing can fold a contour, so a larger window is
    not always better; all windows from ``w0`` to ``wmax`` are tested and the
    best one is used.
    """
    best = None
    window = max(1, int(w0))
    while window <= max(1, int(wmax)):
        Cs = smooth_levels(C, window)
        area, bad = cell_validity(Cs)
        score = mesh_score(Cs, bad)
        if best is None or score < best[0]:
            best = (score, window, Cs, area, bad)
        if score == 0:
            break
        window += 2
    score, window, Cs, area, bad = best
    if score:
        print('WARNING: even the best smoothing window (%d) leaves %d '
              'problem cells.' % (window, score))
    return Cs, window, area, bad


def mesh_score(C, bad=None, ratio=0.5):
    """Degenerate quads plus columns squeezed below ``ratio`` of surface spacing."""
    if bad is None:
        _, bad = cell_validity(C)
    d = np.hypot(np.diff(C[:, :, 0], axis=1), np.diff(C[:, :, 1], axis=1))
    ref = np.median(d[0])
    squeezed = int(np.sum(d[1:] < ratio * ref))
    return int(bad.sum()) + squeezed


def choose_direction(nodes, uvec, cum, w0, wmax):
    """Choose the direction-field smoothing window (direction method)."""
    best = None
    window = max(1, int(w0))
    while window <= max(1, int(wmax)):
        C = direction_contours(nodes, uvec, cum, window)
        area, bad = cell_validity(C)
        score = mesh_score(C, bad)
        if best is None or score < best[0]:
            best = (score, window, C, area, bad)
        if score == 0:
            break
        window += 2
    score, window, C, area, bad = best
    if score:
        print('WARNING: even the best direction window (%d) leaves %d '
              'problem cells.' % (window, score))
    return C, window, area, bad


def patches(C, z):
    """Quad corner table and per-patch best-fit (strike, dip, width)."""
    nlev, nnode = C.shape[0], C.shape[1]
    dz = np.diff(z)
    rows = []
    for i in range(nlev - 1):
        for j in range(nnode - 1):
            p00 = C[i, j]
            p01 = C[i, j + 1]
            p10 = C[i + 1, j]
            p11 = C[i + 1, j + 1]
            top = p01 - p00
            slen = np.hypot(*top)
            nrm = np.array([top[1], -top[0]]) / max(slen, 1e-12)  # right normal
            w1 = np.dot(p10 - p00, nrm)
            w2 = np.dot(p11 - p01, nrm)
            wperp = 0.5 * (w1 + w2)
            width = np.hypot(wperp, dz[i])
            dip = np.degrees(np.arctan2(dz[i], wperp))
            strike = np.degrees(np.arctan2(top[0], top[1])) % 360.0
            rows.append(dict(i=i, j=j, p00=p00, p01=p01, p11=p11, p10=p10,
                             strike=strike, dip=dip, width=width,
                             wperp=wperp))
    return rows


# ---------------------------------------------------------------------------
# writers
# ---------------------------------------------------------------------------
GEOM_FMT = ' '.join(['%12.4f'] * 3 + ['%15.6f'] * 4 + ['%15.3f'] * 3)


def write_faultgeom(path, rows, proj, header):
    nlev = max(r['i'] for r in rows) + 1
    nseg = max(r['j'] for r in rows) + 1
    with open(path, 'w') as fid:
        for line in header:
            fid.write('# ' + line + '\n')
        for i in range(nlev):
            for j in range(nseg):
                r = rows[i * nseg + j]
                lo0, la0 = proj.to_ll(r['p00'][0], r['p00'][1])
                lo1, la1 = proj.to_ll(r['p01'][0], r['p01'][1])
                vals = [r['width'], r['depth'], 180.0 + r['dip'],
                        la0, lo0, la1, lo1, 0.0, 0.0, 0.0]
                fid.write((GEOM_FMT % tuple(vals)) + '\n')


def write_contours(path, C, z, proj, nsegs):
    with open(path, 'w') as fid:
        fid.write('# level depth_km lon lat\n')
        for i in range(C.shape[0]):
            for j in range(C.shape[1]):
                lo, la = proj.to_ll(C[i, j, 0], C[i, j, 1])
                fid.write('%3d %8.3f %12.7f %11.7f\n' % (i, z[i], lo, la))


def write_quads(path, rows, z, proj):
    with open(path, 'w') as fid:
        fid.write('# i j depth_top_km '
                  'lon00 lat00 lon01 lat01 lon11 lat11 lon10 lat10 '
                  'strike_deg dip_deg width_km\n')
        for r in rows:
            lo00, la00 = proj.to_ll(r['p00'][0], r['p00'][1])
            lo01, la01 = proj.to_ll(r['p01'][0], r['p01'][1])
            lo11, la11 = proj.to_ll(r['p11'][0], r['p11'][1])
            lo10, la10 = proj.to_ll(r['p10'][0], r['p10'][1])
            lines = ('%3d %3d %8.3f  %12.7f %11.7f  %12.7f %11.7f  '
                     '%12.7f %11.7f  %12.7f %11.7f  %8.2f %7.2f %7.3f\n')
            fid.write(lines % (r['i'], r['j'], z[r['i']],
                               lo00, la00, lo01, la01, lo11, la11,
                               lo10, la10, r['strike'], r['dip'],
                               r['width']))


def write_triangles(base, C, z, proj):
    """将四边形网格剖分为三角元，输出 TriSlipInversion.py 可读取的文件。

    C = (nlev, nnode, 2) 的局部 (east, north) km 网格节点（即平滑后的
    iso-depth 等高线）。每个四边形沿较短的横对角线分成两个三角形。

    输出：
        <base>_vertex.dat    idx lat lon depth_km（深度向下为正）
        <base>_element.dat   idx1 idx2 idx3 ss ds op（mm，初始滑动为 0）

    三角形顶点按 (v2-v1) x (v3-v1) 的 z 分量为正排列，与
    pyginvc.Geometry.Triangle.load_fault 强制采用的环流方向一致，
    载入时不会被翻转，滑动分量的符号约定因此保持一致。
    """
    nlev, nnode = C.shape[0], C.shape[1]

    def vid(i, j):
        return 1 + i * nnode + j

    tris = []
    for i in range(nlev - 1):
        for j in range(nnode - 1):
            v00, v01, v10, v11 = (i, j), (i, j + 1), (i + 1, j), (i + 1, j + 1)
            # split along the shorter diagonal to keep triangles well shaped
            if (np.hypot(*(C[v01] - C[v10]))
                    < np.hypot(*(C[v00] - C[v11]))):
                cells = ((v01, v11, v10), (v01, v10, v00))
            else:
                cells = ((v00, v01, v11), (v00, v11, v10))
            for a, b, c in cells:
                p1, p2, p3 = C[a], C[b], C[c]
                cross_z = ((p2[0] - p1[0]) * (p3[1] - p1[1])
                           - (p2[1] - p1[1]) * (p3[0] - p1[0]))
                if cross_z < 0.0:              # normal must point downward
                    b, c = c, b
                tris.append((vid(*a), vid(*b), vid(*c)))

    with open(base + '_vertex.dat', 'w') as fid:
        fid.write('# idx lat lon depth_km (depth positive downward)\n')
        for i in range(nlev):
            for j in range(nnode):
                lo, la = proj.to_ll(C[i, j, 0], C[i, j, 1])
                fid.write('%6d %13.7f %13.7f %10.4f\n'
                          % (vid(i, j), la, lo, z[i]))

    with open(base + '_element.dat', 'w') as fid:
        fid.write('# idx1 idx2 idx3 slip_ss_mm slip_ds_mm slip_ts_mm\n')
        for v1, v2, v3 in tris:
            fid.write('%6d %6d %6d %14.4f %14.4f %14.4f\n'
                      % (v1, v2, v3, 0.0, 0.0, 0.0))
    return len(tris)


# ---------------------------------------------------------------------------
# graded triangular mesh (independent of the quad grid)
# ---------------------------------------------------------------------------
def contour_lengths(C):
    """Along-strike length of every depth contour (km)."""
    d = np.hypot(np.diff(C[:, :, 0], axis=1), np.diff(C[:, :, 1], axis=1))
    return d.sum(axis=1)


def contour_point(Ci, s):
    """Position on one contour at fractional arc length ``s`` in [0, 1]."""
    d = np.hypot(np.diff(Ci[:, 0]), np.diff(Ci[:, 1]))
    cum = np.r_[0.0, np.cumsum(d)]
    total = cum[-1]
    if total <= 0.0:
        raise ValueError('zero-length contour; cannot build a mesh')
    sn = np.clip(np.asarray(s, dtype=float), 0.0, 1.0) * total
    return np.interp(sn, cum, Ci[:, 0]), np.interp(sn, cum, Ci[:, 1])


def surface_point(C, zlev, zq, s):
    """Local (east, north) on the fault surface at depth ``zq``.

    ``s`` is the fractional arc length along the nearest depth contours.
    Horizontal positions are interpolated linearly between the two bounding
    contours at the same strike fraction.
    """
    nlev = C.shape[0]
    if nlev < 2:
        raise ValueError('need at least two depth contours')
    # fractional level index at this depth (works for non-uniform levels too)
    f = float(np.interp(zq, zlev, np.arange(nlev)))
    i = min(int(np.floor(f)), nlev - 2)
    a = f - i
    x0, y0 = contour_point(C[i], s)
    x1, y1 = contour_point(C[i + 1], s)
    return (1.0 - a) * x0 + a * x1, (1.0 - a) * y0 + a * y1


def edge_size_profile(zq, z0, z1, l_top, l_bottom, grading):
    """Target 3-D edge length (km) at depth ``zq``."""
    u = float(np.clip((zq - z0) / max(z1 - z0, 1e-12), 0.0, 1.0))
    if grading == 'linear':
        f = u
    elif grading == 'smoothstep':
        f = u * u * (3.0 - 2.0 * u)
    elif grading == 'geometric':
        return float(np.exp(np.log(l_top) * (1.0 - u)
                            + np.log(l_bottom) * u))
    else:
        raise ValueError('unknown tri grading: %s' % grading)
    return l_top + (l_bottom - l_top) * f


def graded_depth_rows(z0, z1, step_at, dip_at, max_rows=200000):
    """Depth rows whose down-dip 3-D spacing follows ``step_at``.

    ``step_at(z)`` is the desired slant distance between neighbouring rows
    (km).  A dip ``theta`` gives ds = dz / sin(theta), so dz = step sin(theta).
    """
    rows = [float(z0)]
    while rows[-1] < z1 - 1e-9:
        zc = rows[-1]
        if z1 - zc < 1e-6:              # avoid a numerically zero last row
            rows[-1] = float(z1)
            break
        # fixed-point iteration for the vertical step; using the dip at the
        # midpoint keeps the 3-D spacing close to the target when dip varies
        dz = float(step_at(zc)) * np.sin(np.radians(float(dip_at(zc))))
        for _ in range(2):
            zm = zc + 0.5 * dz
            dz = float(step_at(zm)) * np.sin(
                np.radians(float(dip_at(zm))))
        if not np.isfinite(dz) or dz <= 1e-9:
            raise ValueError('cannot advance down dip at depth %.3f km' % zc)

        remaining = z1 - zc
        if remaining <= 1.25 * dz:
            rows.append(float(z1))               # accept up to 25% oversize
        elif remaining < 2.0 * dz:
            # split the last two rows evenly instead of leaving a tiny sliver
            rows.append(zc + 0.5 * remaining)
            rows.append(float(z1))
        else:
            rows.append(min(zc + dz, z1))
        if len(rows) > max_rows:
            raise ValueError(
                'graded mesh needs more than %d depth rows; use a larger '
                'edge length' % max_rows)
    if len(rows) == 1:
        rows.append(float(z1))
    return np.array(rows)


def build_graded_triangles(C, z, dips, l_top, l_bottom, grading,
                           max_vertices=2000000):
    """Delaunay triangulation of the fault surface with a target edge size.

    Points are placed on staggered depth rows, i.e. a triangular lattice whose
    nodal spacing is the target edge length ``L``.  Consecutive rows are
    separated down dip by ``L * sqrt(3) / 2`` and alternate rows are shifted
    by half a node, so the Delaunay triangles are close to equilateral.  The
    triangulation is computed in the (along-strike, down-dip) parameter plane,
    whose boundary is a rectangle, and mapped back to local coordinates.  This
    gives an isotropic mesh that is not constrained by the quad grid.

    Returns ``(xy, z, tris, rows, row_counts, edge_lengths, edge_mid_depth)``.
    """
    if not (l_top > 0.0 and l_bottom > 0.0):
        raise ValueError('triangular edge lengths must be positive')
    z0, z1 = float(z[0]), float(z[-1])
    if z1 <= z0:
        raise ValueError('bottom depth must be larger than top depth')

    zmid = 0.5 * (z[:-1] + z[1:])
    dips = np.asarray(dips, dtype=float)
    dip_at = lambda zq: np.interp(zq, zmid, dips)   # scalar or array
    size_at = lambda zq: edge_size_profile(zq, z0, z1, l_top, l_bottom,
                                           grading)
    # triangular lattice: row spacing = L sin(60 deg) gives equilateral cells
    step_at = lambda zq: size_at(zq) * np.sqrt(3.0) / 2.0

    rows = graded_depth_rows(z0, z1, step_at, dip_at)
    lvl_len = contour_lengths(C)
    len_at = lambda zq: float(np.interp(zq, z, lvl_len))
    l_ref = float(np.median(lvl_len))          # constant chart width

    # chart Y: true slant distance from the shallow edge
    slant = np.r_[0.0, np.cumsum(
        np.diff(rows) / np.sin(np.radians(
            0.5 * (dip_at(rows[:-1]) + dip_at(rows[1:])))))]

    chart, xy, zz, row_counts = [], [], [], []
    for k, zq in enumerate(rows):
        nk = int(max(1, round(len_at(zq) / size_at(zq))))
        dx = l_ref / nk
        if k % 2:
            # half-cell offset: interior nodes plus both boundary nodes
            xchart = np.r_[0.0, (np.arange(nk) + 0.5) * dx, l_ref]
        else:
            xchart = np.linspace(0.0, l_ref, nk + 1)
        s = xchart / l_ref
        n = len(xchart)
        row_counts.append(n)
        x, y = surface_point(C, z, zq, s)
        chart.append(np.column_stack([xchart, np.full(n, slant[k])]))
        xy.append(np.column_stack([x, y]))
        zz.append(np.full(n, zq))
    chart = np.vstack(chart)
    xy = np.vstack(xy)
    zz = np.concatenate(zz)
    if len(zz) > max_vertices:
        raise ValueError(
            'graded mesh would have %d vertices; increase --tri-edge '
            '(limit %d)' % (len(zz), max_vertices))

    try:
        from scipy.spatial import Delaunay
    except ImportError:
        raise ImportError('scipy is required for --tri-edge')
    tris = Delaunay(chart).simplices.copy()

    # enforce the same plan-view CCW convention as write_triangles()
    p = xy[tris]
    cross = ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
             - (p[:, 1, 1] - p[:, 0, 1]) * (p[:, 2, 0] - p[:, 0, 0]))
    npos, nneg = int(np.sum(cross > 0.0)), int(np.sum(cross < 0.0))
    if npos and nneg:
        print('WARNING: graded triangulation changes orientation in plan '
              'view (%d CCW, %d CW); the fault surface may be folded.'
              % (npos, nneg))
    tris[cross < 0.0] = tris[cross < 0.0][:, [0, 2, 1]]
    ndeg = int(np.sum(np.abs(cross) < 1e-12))
    if ndeg:
        print('WARNING: %d zero-area triangles in graded mesh.' % ndeg)

    edges = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    P = np.column_stack([xy, zz])
    d = P[edges[:, 0]] - P[edges[:, 1]]
    elen = np.sqrt((d ** 2).sum(axis=1))
    emid = 0.5 * (zz[edges[:, 0]] + zz[edges[:, 1]])
    return xy, zz, tris, rows, row_counts, elen, emid


def triangle_min_angles(xy, zz, tris):
    """Minimum interior angle (deg) of every triangle in 3-D."""
    P = np.column_stack([xy, zz])
    pts = P[tris]
    out = []
    for a, b, c in ((1, 2, 0), (0, 2, 1), (0, 1, 2)):
        u = pts[:, a] - pts[:, c]
        w = pts[:, b] - pts[:, c]
        den = np.linalg.norm(u, axis=1) * np.linalg.norm(w, axis=1)
        cosang = (u * w).sum(axis=1) / np.maximum(den, 1e-30)
        out.append(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))
    return np.min(np.column_stack(out), axis=1)


def write_triangles_graded(base, C, z, dips, proj, l_top, l_bottom, grading):
    """Write the graded Delaunay mesh and return its statistics."""
    info = {}
    (xy, zz, tris, rows, row_counts, elen,
     emid) = build_graded_triangles(C, z, dips, l_top, l_bottom, grading)

    with open(base + '_vertex.dat', 'w') as fid:
        fid.write('# idx lat lon depth_km (depth positive downward)\n')
        for k in range(len(zz)):
            lo, la = proj.to_ll(xy[k, 0], xy[k, 1])
            fid.write('%6d %13.7f %13.7f %10.4f\n'
                      % (k + 1, la, lo, zz[k]))

    with open(base + '_element.dat', 'w') as fid:
        fid.write('# idx1 idx2 idx3 slip_ss_mm slip_ds_mm slip_ts_mm\n')
        for a, b, c in tris:
            fid.write('%6d %6d %6d %14.4f %14.4f %14.4f\n'
                      % (a + 1, b + 1, c + 1, 0.0, 0.0, 0.0))

    minang = triangle_min_angles(xy, zz, tris)
    info.update(xy=xy, z=zz, tris=tris, rows=rows, row_counts=row_counts,
                elen=elen, emid=emid, min_angle=float(minang.min()),
                nvert=len(zz), ntri=len(tris),
                emin=float(elen.min()), emax=float(elen.max()),
                emean=float(elen.mean()),
                nsliver=int(np.sum(minang < 15.0)))
    return info


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------
def encoding_misfit(path, nsegs, ndeps, C, proj):
    """Reload the written faultgeom and compare its deep edge with the mesh."""
    try:
        from pyginvc.Geometry.Patch import Fault
    except Exception as exc:                       # pragma: no cover
        print('WARNING: pyginvc not importable, skip encoding check: %s' % exc)
        return None
    flt = Fault(path, nsegs, ndeps, False)
    flt.load_fault()
    fe = flt.FaultGeom2AllVertex()
    d = []
    for k, (_, r) in enumerate(fe.iterrows()):
        i, j = divmod(k, nsegs)
        seg = (C[i + 1, j], C[i + 1, j + 1])
        for lon, lat in ((r['lt_lon'], r['lt_lat']),
                         (r['rt_lon'], r['rt_lat'])):
            e, n = proj.to_xy(np.array([lon]), np.array([lat]))
            d.append(point_segment_distance(np.array([e[0], n[0]]), seg))
    return np.array(d)


def point_segment_distance(p, seg):
    a, b = seg
    ab = b - a
    t = np.clip(np.dot(p - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0)
    return float(np.hypot(*(p - (a + t * ab))))


def make_figure(path, Cin, Csm, z, rows, proj, args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 8.0),
                             gridspec_kw={'width_ratios': [1.35, 1.0]})
    ax = axes[0]
    cmap = plt.cm.viridis

    for i in range(Cin.shape[0]):
        lo, la = proj.to_ll(Cin[i, :, 0], Cin[i, :, 1])
        ax.plot(lo, la, color='0.72', lw=0.7, ls='--',
                label='initial contours' if i == 1 else None)
    for i in range(Csm.shape[0]):
        lo, la = proj.to_ll(Csm[i, :, 0], Csm[i, :, 1])
        ax.plot(lo, la, color=cmap(i / max(1, Csm.shape[0] - 1)), lw=1.4,
                label='smoothed contours' if i == 1 else None)

    for r in rows:
        quad = np.vstack([r['p00'], r['p01'], r['p11'], r['p10'],
                          r['p00']])
        lo, la = proj.to_ll(quad[:, 0], quad[:, 1])
        ax.plot(lo, la, color='#d62728', lw=0.35, alpha=0.7)

    lo, la = proj.to_ll(Csm[0, :, 0], Csm[0, :, 1])
    ax.plot(lo, la, color='k', lw=2.2, label='surface trace')

    ax.set_aspect(1.0 / np.cos(np.radians(np.mean(la))))
    ax.set_xlabel('longitude (deg)')
    ax.set_ylabel('latitude (deg)')
    ax.set_title('Contours and quad subfaults (colour = smoothed contour)',
                 fontsize=11, loc='left')
    ax.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax.legend(loc='best', fontsize=8, frameon=True)

    # --- right: along-strike spacing per level ---------------------------
    ax2 = axes[1]

    def min_spacing(C):
        d = np.hypot(np.diff(C[:, :, 0], axis=1),
                     np.diff(C[:, :, 1], axis=1))
        return d.min(axis=1), np.median(d, axis=1)

    dmin_i, dmed_i = min_spacing(Cin)
    dmin_s, dmed_s = min_spacing(Csm)
    ax2.plot(z, dmin_i, 'o--', color='0.55', ms=4,
             label='initial: min spacing')
    ax2.plot(z, dmin_s, 'o-', color='#1f77b4', ms=4,
             label='smoothed: min spacing')
    ax2.plot(z, dmed_s, 's-', color='#1f77b4', ms=3, alpha=0.45,
             label='smoothed: median spacing')
    ax2.axhline(0.0, color='#d62728', lw=0.8, ls=':')
    ax2.set_xlabel('depth (km)')
    ax2.set_ylabel('along-strike spacing (km)')
    ax2.set_title('Column spacing versus depth', fontsize=11, loc='left')
    ax2.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax2.legend(loc='best', fontsize=8, frameon=True)
    ax2.invert_yaxis()

    fig.suptitle('%s: nsegs=%d ndeps=%d, %s smoothing=%d'
                 % (os.path.basename(args.trace), args.nsegs, args.ndeps,
                    args.contour_method, args._window), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ('png', 'pdf'):
        fig.savefig('%s.%s' % (path, ext), dpi=170, bbox_inches='tight')
    plt.close(fig)


def make_profile_figure(path, z, dips, args):
    """断层剖面图：深度 vs 距离。

    左：深度 vs 沿倾向的水平距离 cum = sum(dz/tan(dip))；
    右：深度 vs 沿断层面的斜距 s = sum(dz/sin(dip))。
    两图都标注每层倾角，地表在上，listric 形态一目了然。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    zmid = 0.5 * (z[:-1] + z[1:])
    dz = np.diff(z)
    cum = np.r_[0.0, np.cumsum(dz / np.tan(np.radians(dips)))]
    slant = np.r_[0.0, np.cumsum(dz / np.sin(np.radians(dips)))]

    # Give each panel a width proportional to its distance span so both can
    # use the same km-per-inch scale without wasting horizontal space.  A
    # (near-)vertical fault has no horizontal distance: its panel has no
    # meaningful x scale, so reserve a normal-width cell for it below.
    zspan = float(max(abs(z[-1] - z[0]), 1e-9))
    xspans = [abs(cum[-1] - cum[0]), abs(slant[-1] - slant[0])]
    good = [s if s > 1e-6 * zspan else None for s in xspans]
    ref = max([s for s in good if s] + [zspan])
    width_ratios = [s if s else 0.6 * ref for s in good]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True,
                             gridspec_kw={'width_ratios': width_ratios})
    for ax, dist, xlabel in (
            (axes[0], cum, 'horizontal distance along dip (km)'),
            (axes[1], slant, 'distance along fault plane (km)')):
        ax.plot(dist, z, 'o-', color='#1f77b4', lw=1.8, ms=4)
        for i, dip in enumerate(dips):
            ax.annotate('%.0f°' % dip,
                        xy=(0.5 * (dist[i] + dist[i + 1]), zmid[i]),
                        textcoords='offset points', xytext=(5, 0),
                        fontsize=8, va='center', color='0.25')
        ax.set_xlabel(xlabel)
        if abs(dist[-1] - dist[0]) <= 1e-6 * zspan:
            # Vertical fault: every point is at x = 0, so an equal aspect
            # would collapse the axes to a sliver.  Any x scale is truthful
            # here; use a full-width panel with the fault line at zero.
            ax.set_xlim(-1.0, 1.0)
        else:
            # True scale: one kilometre has the same length on both axes,
            # so a 40 km deep fault with 20 km of horizontal offset is
            # drawn 2:1.
            ax.set_aspect('equal', adjustable='box')
        ax.grid(True, ls=':', lw=0.5, alpha=0.6)
    # sharey=True: inverting both axes would cancel out, so invert once
    axes[0].invert_yaxis()
    axes[0].set_ylabel('depth (km)')
    extra = ''
    if args.profile in ('exp', 'log'):
        extra = ', decay=%.1f' % args.profile_decay
    elif args.profile == 'power':
        extra = ', p=%.2f' % args.power
    elif args.profile == 'logistic':
        extra = ', t0=%.2f, s=%.3f' % (args.profile_t0, args.profile_s)
    axes[0].set_title('Fault cross-section (dip profile: %s%s)'
                      % (args.profile, extra),
                      fontsize=11, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig('%s.%s' % (path, ext), dpi=170, bbox_inches='tight')
    plt.close(fig)


def make_tri_figure(path, info, z, l_top, l_bottom, grading, proj):
    """QC figure for the graded triangular mesh.

    Left: map view of the triangles coloured by depth (the fault surface is
    viewed from above, so deeper rows appear next to the trace).  Right: every
    edge length versus its mean depth, with the requested target profile.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    xy, zz, tris = info['xy'], info['z'], info['tris']
    edges = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.4),
                             gridspec_kw={'width_ratios': [1.35, 1.0]})
    ax = axes[0]

    segs = []
    for a, b in edges:
        lo, la = proj.to_ll(xy[[a, b], 0], xy[[a, b], 1])
        segs.append(np.column_stack([lo, la]))
    lc = LineCollection(segs, cmap='viridis', linewidths=0.45, alpha=0.9)
    lc.set_array(0.5 * (zz[edges[:, 0]] + zz[edges[:, 1]]))
    ax.add_collection(lc)
    ax.autoscale()

    n0 = info['row_counts'][0]                    # shallow edge nodes first
    lo, la = proj.to_ll(xy[:n0, 0], xy[:n0, 1])
    ax.plot(lo, la, color='k', lw=2.0, label='surface trace')

    lo_all, la_all = proj.to_ll(xy[:, 0], xy[:, 1])
    ax.set_aspect(1.0 / np.cos(np.radians(np.mean(la_all))))
    ax.set_xlabel('longitude (deg)')
    ax.set_ylabel('latitude (deg)')
    ax.set_title('Graded triangular mesh in map view', fontsize=11, loc='left')
    ax.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax.legend(loc='best', fontsize=8, frameon=True)
    cb = fig.colorbar(lc, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label('depth (km)')

    ax2 = axes[1]
    ax2.scatter(info['elen'], info['emid'], s=4, c=info['emid'],
                cmap='viridis', alpha=0.5, edgecolors='none')
    zzq = np.linspace(z[0], z[-1], 400)
    target = [edge_size_profile(q, z[0], z[-1], l_top, l_bottom, grading)
              for q in zzq]
    ax2.plot(target, zzq, color='#d62728', lw=2.0, label='target')
    for zr in info['rows']:
        ax2.axhline(zr, color='0.85', lw=0.3, zorder=0)
    ax2.set_xlabel('3-D edge length (km)')
    ax2.set_ylabel('depth (km)')
    ax2.set_title('Edge length versus depth', fontsize=11, loc='left')
    ax2.grid(True, ls=':', lw=0.5, alpha=0.6)
    ax2.legend(loc='best', fontsize=8, frameon=True)
    ax2.invert_yaxis()

    fig.suptitle('Graded triangular mesh: %d vertices, %d triangles; '
                 'target %.2f-%.2f km (%s grading)'
                 % (info['nvert'], info['ntri'], l_top, l_bottom, grading),
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ('png', 'pdf'):
        fig.savefig('%s.%s' % (path, ext), dpi=170, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description='Smooth iso-depth contours and build a quad fault mesh. '
                    'Only the fault geometry (--trace/--nsegs/--ndeps/'
                    '--top-depth/--bottom-depth/--top-dip/--bottom-dip/'
                    '--dip-azimuth/--output) is required; meshing details '
                    'have built-in defaults.')
    p.add_argument('--trace', required=True, help='surface trace file')
    p.add_argument('--nsegs', type=int, required=True,
                   help='patches along strike')
    p.add_argument('--ndeps', type=int, required=True,
                   help='patches down dip')
    p.add_argument('--top-depth', type=float, required=True,
                   help='depth of the shallow edge (km)')
    p.add_argument('--bottom-depth', type=float, required=True,
                   help='depth of the deep edge (km)')
    p.add_argument('--top-dip', type=float, required=True,
                   help='dip at the top (deg)')
    p.add_argument('--bottom-dip', type=float, required=True,
                   help='dip at the bottom (deg)')
    p.add_argument('--dip-azimuth', type=float, required=True,
                   help='absolute down-dip azimuth (deg)')
    p.add_argument('--output', required=True,
                   help='output base path (no extension)')
    # --- optional meshing details (sensible built-in values) ------------
    p.add_argument('--strike-mode', default='mean',
                   choices=['segment', 'mean'],
                   help='default mean: every depth contour is a translated '
                        'copy of the surface trace, shifted along the normal '
                        'to the mean strike by dz/tan(dip); segment: extrude '
                        'each patch along its own segment normal')
    p.add_argument('--profile', default='linear',
                   choices=['linear', 'exp', 'power', 'log', 'logistic',
                            'sine', 'cosine', 'smoothstep'],
                   help='dip profile from top to bottom (default linear). '
                        'exp = steep at the top and progressively gentler '
                        'down dip (top-dip must be larger than bottom-dip). '
                        'power = t**p, p<1 changes fast shallow and keeps '
                        'changing down dip, p>1 changes fast deep. '
                        'logistic = turns over sharply around t0')
    p.add_argument('--profile-decay', type=float, default=3.0,
                   help='exp/log profile: how fast the dip approaches the '
                        'bottom value (default 3.0; larger = flattens sooner)')
    p.add_argument('--power', type=float, default=1.0,
                   help='power profile exponent p in f=t**p (default 1.0). '
                        'p<1: rapid change at the top, sustained change with '
                        'depth; p>1: flat shallow, rapid change at depth')
    p.add_argument('--profile-t0', type=float, default=0.5,
                   help='logistic profile: normalized depth (0-1) of the '
                        'steepest dip change (default 0.5)')
    p.add_argument('--profile-s', type=float, default=0.1,
                   help='logistic profile: transition width in normalized '
                        'depth (default 0.1; smaller = sharper bend)')
    p.add_argument('--smooth', type=int, default=10,
                   help='surface trace smoothing window (default 10; 1 = none)')
    p.add_argument('--contour-method', default='direction',
                   choices=['direction', 'nodes'],
                   help='contour smoothing method (default direction)')
    p.add_argument('--contour-smooth-km', type=float, default=None,
                   help='contour smoothing length in km (default: auto)')
    p.add_argument('--tri', action='store_true',
                   help='also write a triangular mesh (base_vertex.dat, '
                        'base_element.dat) readable by TriSlipInversion.py')
    p.add_argument('--tri-edge', type=float, default=None, metavar='KM',
                   help='target 3-D edge length (km) of an independently '
                        're-meshed triangular grid; implies --tri.  Serves '
                        'as the default for --tri-edge-top/-bottom')
    p.add_argument('--tri-edge-top', type=float, default=None, metavar='KM',
                   help='target edge length (km) at the shallow edge')
    p.add_argument('--tri-edge-bottom', type=float, default=None, metavar='KM',
                   help='target edge length (km) at the deep edge')
    p.add_argument('--tri-grading', default='geometric',
                   choices=['geometric', 'linear', 'smoothstep'],
                   help='how the target edge length transitions from the '
                        'shallow to the deep edge (default geometric: '
                        'constant ratio between neighbouring rows)')
    args = p.parse_args()

    # --tri-edge is the default for both ends; specifying either end alone
    # makes the mesh uniform at that size.
    l_top = args.tri_edge_top if args.tri_edge_top is not None else args.tri_edge
    l_bot = (args.tri_edge_bottom if args.tri_edge_bottom is not None
             else args.tri_edge)
    if l_top is not None and l_bot is None:
        l_bot = l_top
    if l_bot is not None and l_top is None:
        l_top = l_bot
    if l_top is not None:
        if l_top <= 0.0 or l_bot <= 0.0:
            p.error('triangular edge lengths must be positive')
        args.tri = True

    proj, A, B, uvec, total, az0, reverse = build_surface(args)
    contour_smooth, contour_smooth_max = 5, 41
    if args.contour_smooth_km is not None:
        node_km = total / max(1, args.nsegs)
        w = int(round(args.contour_smooth_km / node_km))
        contour_smooth = max(1, w + (1 - w % 2))
        print('contour smoothing length %.1f km -> window %d (node spacing '
              '%.2f km)' % (args.contour_smooth_km, contour_smooth, node_km))
    z, dips = depth_levels(args)
    if args.strike_mode == 'mean':
        # 深部各层 = 地表迹线的纯平移（沿平均走向法向 + 垂向深度），
        # 形状与总长度保持不变，不会随深度扩散
        Csm, cum, u_mean = build_mean_strike_mesh(
            A, B, z, dips, az0, args.dip_azimuth)
        Cin = Csm
        window = 1
        area_s, bad_s = cell_validity(Csm)
        lengths = [np.sum(np.hypot(np.diff(Cin[i, :, 0]),
                                   np.diff(Cin[i, :, 1])))
                   for i in range(len(z))]
        print('mean-strike mesh: %d contours, each a translation of the '
              'surface trace (length %.2f-%.2f km, identical shape)'
              % (len(z), min(lengths), max(lengths)))
        print('mean-strike shift: horizontal offset %.2f km at %.1f km depth, '
              'direction azimuth %.1f deg'
              % (cum[-1], z[-1], np.degrees(np.arctan2(u_mean[0], u_mean[1]))
                 % 360.0))
    else:
        Cin, cum = initial_contours(A, B, uvec, z, dips)
    print('trace length %.2f km, %d segments, %d depth levels'
          % (total, args.nsegs, len(z)))
    print('surface azimuth %.1f deg (%s strike mode), dip azimuth %.1f deg%s'
          % (az0, args.strike_mode, args.dip_azimuth,
             ' (trace reversed)' if reverse else ''))

    area_i, bad_i = cell_validity(Cin)
    nbad = int(bad_i.sum())
    if nbad:
        i0 = int(np.where(bad_i.any(axis=1))[0][0])
        print('initial mesh: %d degenerate quads, first at level %d (%.1f km)'
              % (nbad, i0, z[i0]))
    else:
        print('initial mesh: no degenerate quads')

    if args.strike_mode == 'mean':
        print('mean-strike mode: depth contours are translated copies of the '
              'surface trace (contour smoothing skipped)')
    elif args.contour_method == 'direction':
        nodes = np.vstack([A, B[-1][None, :]])
        Csm, window, area_s, bad_s = choose_direction(
            nodes, uvec, cum, contour_smooth, contour_smooth_max)
    else:
        Csm, window, area_s, bad_s = choose_smoothing(
            Cin, contour_smooth, contour_smooth_max)
    args._window = window
    print('%s smoothing window %d: %d degenerate quads'
          % ('mean' if args.strike_mode == 'mean' else args.contour_method,
             window, int(bad_s.sum())))

    dmin = np.min(np.hypot(np.diff(Csm[:, :, 0], axis=1),
                           np.diff(Csm[:, :, 1], axis=1)))
    print('smoothed mesh: min along-strike spacing %.3f km' % dmin)

    # ---- quad mesh rows -------------------------------------------------
    rows = patches(Csm, z)
    for r in rows:
        r['depth'] = z[r['i']]

    # ---- initial faultgeom ---------------------------------------------
    if args.strike_mode == 'mean':
        # mean 模式下"初始"几何就是平移拷贝，与最终网格完全一致
        rows_i = rows
        header_i = ['quad fault mesh built by smooth_fault_mesh.py',
                    'trace      : %s (%d segments along strike)' % (
                        args.trace, args.nsegs),
                    'grid       : nsegs=%d ndeps=%d' % (args.nsegs,
                                                       args.ndeps),
                    'geometry   : mean-strike translation of the surface '
                    'trace (same as smoothed mesh)',
                    'depth      : %.3f - %.3f km' % (z[0], z[-1])]
    else:
        # Exactly reproduce trace2faultgeom.py: every patch edge is shifted
        # along its own segment normal, so adjacent patches may gap/overlap.
        dz_all = np.diff(z)
        widths_i = dz_all / np.sin(np.radians(dips))
        rows_i = []
        for i in range(args.ndeps):
            for j in range(args.nsegs):
                p00 = A[j] + cum[i] * uvec[j]
                p01 = B[j] + cum[i] * uvec[j]
                top = p01 - p00
                rows_i.append(dict(
                    i=i, j=j, p00=p00, p01=p01, p10=None, p11=None,
                    width=widths_i[i], dip=dips[i],
                    strike=np.degrees(np.arctan2(top[0], top[1])) % 360.0,
                    wperp=dz_all[i] / np.tan(np.radians(dips[i])),
                    depth=z[i]))
        header_i = ['quad fault mesh built by smooth_fault_mesh.py',
                    'trace      : %s (%d segments along strike)' % (
                        args.trace, args.nsegs),
                    'grid       : nsegs=%d ndeps=%d' % (args.nsegs,
                                                       args.ndeps),
                    'geometry   : initial straight-ray extrusion (may cross)',
                    'depth      : %.3f - %.3f km' % (z[0], z[-1])]
    write_faultgeom(args.output + '_initial.faultgeom', rows_i, proj,
                    header_i)

    # ---- smoothed quad mesh --------------------------------------------
    header_s = ['quad fault mesh built by smooth_fault_mesh.py',
                'trace      : %s (%d segments along strike)' % (args.trace,
                                                               args.nsegs),
                'grid       : nsegs=%d ndeps=%d' % (args.nsegs, args.ndeps),
                'geometry   : %s' % (
                    'mean-strike translation of the surface trace'
                    if args.strike_mode == 'mean'
                    else 'smoothed iso-depth contours -> quad subfaults'),
                'smoothing  : %s, window=%d (surface window=%d)'
                % ('translation copy' if args.strike_mode == 'mean'
                   else args.contour_method + ' method', window, args.smooth),
                'depth      : %.3f - %.3f km' % (z[0], z[-1]),
                'dip        : %.2f - %.2f deg'
                % (np.min(dips), np.max(dips)),
                'convention : shallow edge, dip column = 180 + angle']
    geom_s = args.output + '_smoothed.faultgeom'
    write_faultgeom(geom_s, rows, proj, header_s)

    write_contours(args.output + '_contours.dat', Csm, z, proj, args.nsegs)
    write_quads(args.output + '_quads.dat', rows, z, proj)
    tri_info = None
    if args.tri:
        if l_top is not None:
            tri_info = write_triangles_graded(
                args.output, Csm, z, dips, proj, l_top, l_bot,
                args.tri_grading)
            print('graded triangular mesh: %d vertices, %d triangles; '
                  'target edge %.2f-%.2f km (%s), actual %.2f-%.2f km '
                  '(mean %.2f); min angle %.1f deg'
                  % (tri_info['nvert'], tri_info['ntri'], l_top, l_bot,
                     args.tri_grading, tri_info['emin'], tri_info['emax'],
                     tri_info['emean'], tri_info['min_angle']))
            if tri_info['nsliver']:
                print('WARNING: %d triangles have a minimum angle below '
                      '15 deg; try a larger --tri-edge or smoother grading.'
                      % tri_info['nsliver'])
        else:
            ntris = write_triangles(args.output, Csm, z, proj)
            print('triangular mesh: %d triangles, %d vertices'
                  % (ntris, Csm.shape[0] * Csm.shape[1]))

    # ---- encoding check -------------------------------------------------
    d = encoding_misfit(geom_s, args.nsegs, args.ndeps, Csm, proj)
    if d is not None:
        print('faultgeom parallelogram vs exact quads: '
              'median %.1f m, max %.1f m' % (1000 * np.median(d),
                                             1000 * d.max()))

    make_figure(args.output + '_mesh', Cin, Csm, z, rows, proj, args)
    make_profile_figure(args.output + '_profile', z, dips, args)
    if tri_info is not None:
        make_tri_figure(args.output + '_tri_mesh', tri_info, z, l_top, l_bot,
                        args.tri_grading, proj)
    written = [args.output + '_initial.faultgeom',
               args.output + '_contours.dat',
               args.output + '_quads.dat',
               args.output + '_smoothed.faultgeom',
               args.output + '_mesh.png',
               args.output + '_profile.png']
    if args.tri:
        written += [args.output + '_vertex.dat',
                    args.output + '_element.dat']
    if tri_info is not None:
        written += [args.output + '_tri_mesh.png']
    print('written ' + ', '.join(written))


if __name__ == '__main__':
    main()
