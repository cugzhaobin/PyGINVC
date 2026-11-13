#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
trace2faultgeom.py
==================

Build a rectangular-patch faultgeom file from a **surface fault trace**.

The trace (e.g. ``fault_trace.dat``, two columns of lon/lat) is discretised
along strike into ``nsegs`` patches.  Every patch is extended down dip with a
constant or depth-dependent dip, producing the PyGINV faultgeom columns::

    width_km  depth_km  dip_deg  lat0  lon0  lat1  lon1  ss  ds  op

Construction
------------
For row ``i`` (from the top of the fault downwards)::

    dz_i    = z[i+1] - z[i]                      vertical extent (km)
    theta_i = dip of row i                       (degrees from horizontal)
    width_i = dz_i / sin(theta_i)                down-dip patch width (km)
    dx_i    = dz_i / tan(theta_i)                horizontal step (km)

Row ``i`` is attached to the previous one (the bottom edge of row ``i`` is the
top edge of row ``i+1``), so the fault surface is continuous.  The down-dip
direction is perpendicular to the local trace (to the *right* of the along
strike direction A -> B by default; use ``--dip-side left`` to flip it).

Depth range and patch width
---------------------------
Two mutually exclusive ways to size the fault:

* **Depth-controlled** (``--bottom-depth`` / ``--depth-nodes``): the top and
  bottom depths are fixed, the dip gives the widths.
* **Width-controlled** (``--mean-width`` and no bottom depth): every patch has
  the given down-dip width and the bottom depth follows from the dips.

If both ``--bottom-depth`` and ``--mean-width`` are given, the depth range is
kept and the *absolute dip level* is adjusted (``--solve-dip``) so that the
mean patch width equals ``--mean-width`` exactly.

Dip profiles
------------
``t`` is the normalised position from the top (0) to the bottom (1) of the
fault.  The dip varies from ``--top-dip`` at t=0 to ``--bottom-dip`` at t=1.
Available shapes f(t) (all run from 0 at t=0 to 1 at t=1)::

    linear      f = t
    power       f = t**p                        (--power p)
    exp         f = (1 - exp(-k t)) / (1 - exp(-k))       (--decay k)
    log         f = log(1 + k t) / log(1 + k)             (--decay k)
    logistic    f = normalised 1/(1+exp((t-t0)/s))        (--t0, --s)
    sine        f = sin(pi t / 2)
    cosine      f = (1 - cos(pi t)) / 2
    smoothstep  f = t**2 (3 - 2 t)

    dip(t) = top_dip + (bottom_dip - top_dip) * f(t)

The profile is evaluated at the mid-depth of every row.  A single dip is
given with ``--dip`` (or ``--top-dip = --bottom-dip``).  An
explicit per-row list overrides everything: ``--dips 60,52,46,41,37,34``.

Dip convention
--------------
By default (``--dip-mode reverse``) the given ``(lat0, lon0)-(lat1, lon1)``
edge is the **shallow** edge (the surface trace for row 0) and the dip column
is written as ``180 + dip``.  This is the convention used by the existing
Dingri faultgeom files.  ``--dip-mode angle`` writes the true angle (0-90) and
the given edge is the deep edge instead.

Important: PyGINV always puts the down-dip side to the **right** of the A -> B
strike direction.  Therefore a west-dipping fault whose trace is stored from
south to north must be written with the trace reversed; ``--dip-side left``
does this automatically (``--dip-azimuth`` picks the side from an absolute
azimuth).

Examples
--------
Single 49-degree dip, trace segmented into 33 patches, 13 rows to 25 km::

    ./trace2faultgeom.py --trace fault_trace.dat --nsegs 33 --ndeps 13 \\
        --dip 49 --top-depth 0.9 --bottom-depth 25.4 \\
        --output dingri_49.faultgeom

Listric fault: 60 deg at the surface, 32 deg at 18 km, exponential::

    ./trace2faultgeom.py --trace fault_trace.dat --nsegs 33 --ndeps 13 \\
        --top-dip 60 --bottom-dip 32 --profile exp --decay 4 \\
        --top-depth 0 --bottom-depth 18 --dip-side left \\
        --output dingri_listric.faultgeom

Uniform 2.5 km patches, dip solved from the mean width::

    ./trace2faultgeom.py --trace fault_trace.dat --nsegs 33 --ndeps 13 \\
        --mean-width 2.5 --top-depth 0.9 --bottom-depth 25.4 \\
        --solve-dip --output dingri_w25.faultgeom

Plotting: ``--plot`` saves the 2D QC figure; ``--plot3d`` saves a 3D patch
view rendered by pyginvc's own ``Fault.plot_faultgeom`` (same look as
``plot_faultgeom.py``), optionally with aftershocks and patch indices.

Dependencies: numpy, pyproj (only for the local projection); matplotlib for
``--plot``; pyginvc (pandas) for ``--plot3d``.
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
# small helpers
# ---------------------------------------------------------------------------
def angular_diff(a, b):
    """Signed smallest difference a - b wrapped into [-180, 180)."""
    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0


def mean_azimuth(east, north):
    """Circular mean of azimuths (degrees clockwise from north)."""
    az = np.degrees(np.arctan2(east, north))
    z = np.mean(np.exp(1j * np.radians(az)))
    return np.degrees(np.angle(z)) % 360.0


def parse_floats(text):
    return np.array([float(v) for v in text.replace(' ', '').split(',')
                     if v != ''], dtype=float)


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


# ---------------------------------------------------------------------------
# trace input / preparation
# ---------------------------------------------------------------------------
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


def solve_dip_scale(base_dips, dz, target_mean):
    """Scale dips in (90 - dip) space until the mean width is target_mean.

    dip_i(scale) = 90 - scale * (90 - base_dips_i).  scale=1 reproduces the
    input profile; larger scale => gentler dips => wider patches.
    """
    base = np.asarray(base_dips, dtype=float)
    dz = np.asarray(dz, dtype=float)

    def mean_width(scale):
        theta = np.radians(90.0 - scale * (90.0 - base))
        if np.any(theta <= 0.0) or np.any(theta > np.pi / 2.0):
            return np.inf
        return np.mean(dz / np.sin(theta))

    lo = 1.0e-8
    if np.min(base) >= 90.0:
        s_max = np.inf
    else:
        s_max = (90.0 - 1.0e-3) / (90.0 - np.min(base))
    hi = 1.0
    if np.isfinite(s_max):
        hi = min(hi, s_max)
    w_lo = mean_width(lo)
    w_hi = mean_width(hi)
    while w_hi < target_mean and np.isfinite(s_max) and hi < s_max:
        hi = min(hi * 2.0, s_max)
        w_hi = mean_width(hi)
    if w_lo > target_mean:
        raise ValueError(
            "cannot reach mean width %.3f km: with this depth range the "
            "smallest possible mean width is %.3f km"
            % (target_mean, w_lo))
    if w_hi < target_mean:
        raise ValueError(
            "cannot reach mean width %.3f km: with these dips the largest "
            "reachable mean width is %.3f km" % (target_mean, w_hi))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mean_width(mid) < target_mean:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-12 * max(1.0, hi):
            break
    scale = 0.5 * (lo + hi)
    return 90.0 - scale * (90.0 - base), scale


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="Build a faultgeom file from a surface fault trace.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="See the module docstring for the construction and examples.")

    p.add_argument('--trace', required=True,
                   help='surface trace file, columns lon lat (or lat lon)')
    p.add_argument('--output', default='trace.faultgeom',
                   help='output faultgeom file (default: trace.faultgeom)')
    p.add_argument('--order', choices=['auto', 'lonlat', 'latlon'],
                   default='auto', help='trace column order (default: auto)')

    # along strike
    p.add_argument('--nsegs', type=int, default=0,
                   help='number of patches along strike (default: trace '
                        'vertices, i.e. one patch per input segment)')
    p.add_argument('--seg-length', type=float, default=None,
                   help='target patch length along strike (km); ignored when '
                        '--nsegs is given')
    p.add_argument('--smooth', type=int, default=0,
                   help='moving-average window (vertices) for the trace; '
                        '0 = no smoothing (default: 0)')

    # along dip
    p.add_argument('--ndeps', type=int, default=0,
                   help='number of rows along dip (default: 10, or the length '
                        'of --dips / --depth-nodes)')
    p.add_argument('--top-depth', type=float, default=0.0,
                   help='depth of the top of the fault (km, default: 0)')
    p.add_argument('--bottom-depth', type=float, default=None,
                   help='depth of the bottom of the fault (km); omit it to use '
                        'width-controlled mode via --mean-width')
    p.add_argument('--depth-nodes', type=str, default=None,
                   help='explicit depth nodes (km, comma separated, strictly '
                        'increasing); overrides --top-depth/--bottom-depth')
    p.add_argument('--mean-width', type=float, default=None,
                   help='mean down-dip width of a patch (km).  With a fixed '
                        'depth range it is a target (see --solve-dip); without '
                        'a bottom depth every patch gets exactly this width')
    p.add_argument('--solve-dip', action='store_true',
                   help='with a fixed depth range and --mean-width, rescale the '
                        'dip profile so the mean width matches exactly')

    # dip
    p.add_argument('--dip', type=float, default=None,
                   help='constant dip (degrees from horizontal)')
    p.add_argument('--top-dip', type=float, default=None,
                   help='dip at the top of the fault (degrees)')
    p.add_argument('--bottom-dip', type=float, default=None,
                   help='dip at the bottom of the fault (degrees)')
    p.add_argument('--dips', type=str, default=None,
                   help='explicit dip for every row, top to bottom, comma '
                        'separated; overrides the profile')
    p.add_argument('--profile', default='linear',
                   choices=['linear', 'power', 'exp', 'log', 'logistic',
                            'sine', 'cosine', 'smoothstep'],
                   help='dip profile shape (default: linear)')
    p.add_argument('--power', type=float, default=1.0,
                   help='exponent of the power profile (default: 1.0)')
    p.add_argument('--decay', '--k', dest='decay', type=float, default=1.0,
                   help='decay parameter k of the exp/log profiles '
                        '(default: 1.0)')
    p.add_argument('--t0', type=float, default=0.5,
                   help='logistic transition centre, normalised depth '
                        '(default: 0.5)')
    p.add_argument('--s', type=float, default=0.1,
                   help='logistic transition width, normalised depth '
                        '(default: 0.1)')

    # orientation / convention
    p.add_argument('--dip-side', choices=['right', 'left'], default='right',
                   help='down-dip side relative to the A -> B trace order: '
                        'right (default) or left.  PyGINV always dips to the '
                        'right, so "left" reverses the written trace order.')
    p.add_argument('--dip-azimuth', type=float, default=None,
                   help='desired down-dip azimuth (degrees, clockwise from '
                        'north); chooses --dip-side automatically')
    p.add_argument('--dip-mode', choices=['reverse', 'angle'],
                   default='reverse',
                   help='reverse: given edge is the shallow edge and the dip '
                        'column is 180+angle (default); angle: given edge is '
                        'the deep edge and the dip column is the true angle')
    p.add_argument('--flip-dip-direction', action='store_true',
                   help='alias for --dip-side left')

    # slip / output
    p.add_argument('--ss', type=float, default=0.0,
                   help='constant strike-slip value written to every row (mm)')
    p.add_argument('--ds', type=float, default=0.0,
                   help='constant dip-slip value written to every row (mm)')
    p.add_argument('--op', type=float, default=0.0,
                   help='constant opening value written to every row (mm)')
    p.add_argument('--no-header', action='store_true',
                   help='do not write the leading comment header')
    p.add_argument('--plot', nargs='?', const='', default=None,
                   help='save a 2D QC figure (default name: '
                        'trace2faultgeom.png)')
    p.add_argument('--plot-column', type=int, default=0,
                   help='along-strike column used for the cross-section '
                        '(default: 0)')

    # ---- 3D patch view (in the spirit of plot_faultgeom.py) -------------
    p.add_argument('--plot3d', nargs='?', const='', default=None,
                   help='save a 3D patch view like plot_faultgeom.py '
                        '(default name: <output>_3d.png)')
    p.add_argument('--color-by', choices=['depth', 'dip', 'width', 'slip'],
                   default='depth',
                   help='quantity used to colour the patches (default: depth)')
    p.add_argument('--view-azimuth', type=float, default=-40.0,
                   help='3D view azimuth (degrees, default: -40)')
    p.add_argument('--view-elevation', type=float, default=30.0,
                   help='3D view elevation (degrees, default: 30)')
    p.add_argument('--show-index', action='store_true',
                   help='label every patch with its index in the 3D view')
    p.add_argument('--aftershock', type=str, default='',
                   help='optional aftershock catalog (lon lat depth_km) '
                        'overlaid on the 3D view')
    p.add_argument('--plot3d-show', action='store_true',
                   help='open an interactive window in addition to saving')
    return p


def resolve_arguments(args):
    """Apply defaults, validate combinations and return a settings dict."""
    # ---- dip levels -----------------------------------------------------
    explicit_dips = None
    if args.dips:
        explicit_dips = parse_floats(args.dips)
        if np.any(explicit_dips <= 0.0) or np.any(explicit_dips > 90.0):
            raise ValueError("--dips must satisfy 0 < dip <= 90")
        if (args.dip is not None or args.top_dip is not None or
                args.bottom_dip is not None):
            print("NOTE: --dips overrides --dip/--top-dip/--bottom-dip.",
                  file=sys.stderr)

    if args.dip is not None:
        if args.top_dip is not None or args.bottom_dip is not None:
            raise ValueError("use either --dip or --top-dip/--bottom-dip, "
                             "not both")
        if explicit_dips is None:
            if not (0.0 < args.dip <= 90.0):
                raise ValueError("--dip must satisfy 0 < dip <= 90")

    elif explicit_dips is None:
        if (args.top_dip is None) != (args.bottom_dip is None):
            raise ValueError("--top-dip and --bottom-dip must be given "
                             "together")
        if args.top_dip is not None:
            if not (0.0 < args.top_dip <= 90.0 and
                    0.0 < args.bottom_dip <= 90.0):
                raise ValueError("dips must satisfy 0 < dip <= 90")

    # ---- number of rows -------------------------------------------------
    candidates = []
    if args.ndeps:
        candidates.append(('--ndeps', args.ndeps))
    if explicit_dips is not None:
        candidates.append(('--dips', len(explicit_dips)))
    if args.depth_nodes:
        nodes_tmp = parse_floats(args.depth_nodes)
        if len(nodes_tmp) < 2:
            raise ValueError("--depth-nodes needs at least two values")
        if np.any(np.diff(nodes_tmp) <= 0.0):
            raise ValueError("--depth-nodes must be strictly increasing")
        candidates.append(('--depth-nodes', len(nodes_tmp) - 1))
    if len({n for _, n in candidates}) > 1:
        raise ValueError("conflicting row counts: %s"
                         % ', '.join('%s=%d' % c for c in candidates))
    n_out = candidates[0][1] if candidates else 10
    if n_out < 1:
        raise ValueError("need at least one row along dip")
    if n_out > 1000:
        raise ValueError("more than 1000 rows along dip is not supported")

    # ---- depth mode vs width mode --------------------------------------
    explicit_depth = args.depth_nodes is not None or \
        args.bottom_depth is not None
    width_mode = (not explicit_depth) and (args.mean_width is not None)
    if not explicit_depth and args.mean_width is None:
        raise ValueError(
            "give a depth range (--bottom-depth or --depth-nodes) or a patch "
            "width (--mean-width)")
    if width_mode and args.solve_dip:
        raise ValueError("--solve-dip needs a fixed depth range; give "
                         "--bottom-depth")
    if args.solve_dip and args.mean_width is None:
        raise ValueError("--solve-dip needs --mean-width")
    if args.mean_width is not None and args.mean_width <= 0.0:
        raise ValueError("--mean-width must be > 0")

    # ---- depth nodes ----------------------------------------------------
    if args.depth_nodes:
        z_nodes = parse_floats(args.depth_nodes)
        z_top, z_bot = float(z_nodes[0]), float(z_nodes[-1])
        if args.top_depth != 0.0 and abs(args.top_depth - z_top) > 1.0e-6:
            print("NOTE: --top-depth overridden by --depth-nodes.",
                  file=sys.stderr)
    elif explicit_depth:
        z_top = float(args.top_depth)
        z_bot = float(args.bottom_depth)
        if z_bot <= z_top:
            raise ValueError("--bottom-depth (%.3f) must be larger than the "
                             "top depth (%.3f)" % (z_bot, z_top))
        z_nodes = np.linspace(z_top, z_bot, n_out + 1)
    else:
        z_top = float(args.top_depth)
        z_nodes = None
        z_bot = None

    # ---- nominal dips (before any width solving) ------------------------
    dips = None
    if explicit_dips is not None:
        dips = explicit_dips.copy()
    elif args.dip is not None:
        dips = np.full(n_out, float(args.dip))
    elif args.top_dip is not None:
        if width_mode:
            t = (np.arange(n_out) + 0.5) / n_out
        else:
            t = (0.5 * (z_nodes[:-1] + z_nodes[1:]) - z_top) / (z_bot - z_top)
        dips = dip_profile(t, float(args.top_dip), float(args.bottom_dip),
                           args.profile, args.power, args.decay, args.t0,
                           args.s)
    elif not width_mode and args.mean_width is not None:
        # constant dip solved from the depth range and the mean width
        if args.profile != 'linear':
            raise ValueError(
                "a dip profile needs --top-dip and --bottom-dip; with only "
                "--mean-width a constant dip is solved")
        sin_dip = (z_bot - z_top) / (n_out * args.mean_width)
        if sin_dip > 1.0:
            raise ValueError(
                "mean width %.3f km is too small for a depth range of %.3f km "
                "and %d rows: the minimum mean width is %.3f km"
                % (args.mean_width, z_bot - z_top, n_out,
                   (z_bot - z_top) / n_out))
        theta = np.degrees(np.arcsin(sin_dip))
        dips = np.full(n_out, theta)
        print("NOTE: constant dip solved from depth range and mean width: "
              "%.4f deg" % theta, file=sys.stderr)
    else:
        raise ValueError(
            "no dip given: use --dip, --top-dip/--bottom-dip, --dips, or "
            "give --mean-width with a depth range for a constant dip")

    # ---- widths / depths ------------------------------------------------
    solve_scale = 1.0
    if width_mode:
        if dips is None:
            raise ValueError("width-controlled mode needs --dip or "
                             "--top-dip/--bottom-dip")
        widths = np.full(n_out, float(args.mean_width))
        z_nodes = z_top + np.r_[0.0, np.cumsum(
            widths * np.sin(np.radians(dips)))]
        z_bot = float(z_nodes[-1])
        dz = np.diff(z_nodes)
    else:
        dz = np.diff(z_nodes)
        if args.mean_width is not None and args.solve_dip:
            if explicit_dips is not None or args.dip is not None:
                raise ValueError("--solve-dip needs a dip profile; it cannot "
                                 "rescale an explicit --dip / --dips")
            dips, solve_scale = solve_dip_scale(
                dips, dz, float(args.mean_width))
            print("NOTE: dip profile rescaled by %.6f in (90-dip) space to hit "
                  "the requested mean width." % solve_scale,
                  file=sys.stderr)
        widths = dz / np.sin(np.radians(dips))
        if args.mean_width is not None:
            achieved = float(np.mean(widths))
            rel = abs(achieved - args.mean_width) / args.mean_width
            if rel > 1.0e-3:
                print("WARNING: mean patch width is %.4f km, requested %.4f km "
                      "(%.2f%% off).  Use --solve-dip to scale the dip "
                      "profile." % (achieved, args.mean_width, 100.0 * rel),
                      file=sys.stderr)

    if np.any(dips <= 0.0) or np.any(dips > 90.0):
        raise ValueError("dip values must satisfy 0 < dip <= 90 degrees")
    if np.any(widths <= 0.0):
        raise ValueError("computed non-positive patch widths")

    return dict(n_out=n_out, z_nodes=np.asarray(z_nodes, dtype=float),
                z_top=z_top, z_bot=z_bot, dz=np.asarray(dz, dtype=float),
                dips=np.asarray(dips, dtype=float),
                widths=np.asarray(widths, dtype=float),
                width_mode=width_mode, solve_scale=solve_scale,
                explicit_dips=explicit_dips is not None)


def build_geometry(args, cfg):
    """Create the output faultgeom array together with QC information."""
    lat, lon, _ = read_trace(args.trace, args.order)
    proj = LocalProjection(np.mean(lat), np.mean(lon))
    x, y = proj.to_xy(lon, lat)
    x, y = smooth_xy(x, y, args.smooth)

    # ---- along-strike segmentation -------------------------------------
    if args.nsegs:
        nsegs = int(args.nsegs)
        if nsegs < 1:
            raise ValueError("--nsegs must be >= 1")
        xn, yn, total_len = resample_nodes(x, y, nsegs)
    elif args.seg_length is not None:
        total_len = float(np.sum(np.hypot(np.diff(x), np.diff(y))))
        nsegs = max(1, int(round(total_len / args.seg_length)))
        xn, yn, total_len = resample_nodes(x, y, nsegs)
        print("NOTE: --seg-length %.3f km -> %d segments (%.3f km each)."
              % (args.seg_length, nsegs, total_len / nsegs), file=sys.stderr)
    else:
        xn, yn = x, y
        nsegs = len(xn) - 1
        total_len = float(np.sum(np.hypot(np.diff(xn), np.diff(yn))))

    # ---- choose the down-dip side --------------------------------------
    az0 = mean_azimuth(np.diff(xn), np.diff(yn))
    reverse = (args.dip_side == 'left') or args.flip_dip_direction
    chosen_az = None
    if args.dip_azimuth is not None:
        d_right = abs(angular_diff((az0 + 90.0) % 360.0, args.dip_azimuth))
        d_left = abs(angular_diff((az0 - 90.0) % 360.0, args.dip_azimuth))
        reverse = d_left < d_right
        chosen_az = (az0 - 90.0) % 360.0 if reverse else (az0 + 90.0) % 360.0
        if min(d_right, d_left) > 45.0:
            print("WARNING: the trace is not perpendicular to the requested "
                  "dip azimuth (best mismatch %.1f deg); the down-dip "
                  "direction follows the local trace." % min(d_right, d_left),
                  file=sys.stderr)
    if reverse:
        xn = xn[::-1].copy()
        yn = yn[::-1].copy()

    A = np.column_stack([xn[:-1], yn[:-1]])
    B = np.column_stack([xn[1:], yn[1:]])
    seg = B - A
    slen = np.hypot(seg[:, 0], seg[:, 1])
    if np.any(slen <= 1.0e-9):
        raise ValueError("degenerate segment along strike after resampling")
    tvec = seg / slen[:, None]
    # down-dip unit vector: rotate the along-strike vector clockwise by 90 deg
    uvec = np.column_stack([tvec[:, 1], -tvec[:, 0]])
    strikes = np.degrees(np.arctan2(tvec[:, 0], tvec[:, 1])) % 360.0
    dip_az = np.degrees(np.arctan2(uvec[:, 0], uvec[:, 1])) % 360.0

    # ---- rows along dip ------------------------------------------------
    n_out = cfg['n_out']
    dips = cfg['dips']
    widths = cfg['widths']
    z_nodes = cfg['z_nodes']
    dz = cfg['dz']
    reverse_dip = (args.dip_mode == 'reverse')

    out = np.zeros((n_out * nsegs, 10))
    top_xy = A.copy()
    rows_xy = []                      # top edge of every row, for QC
    for i in range(n_out):
        hstep = dz[i] / np.tan(np.radians(dips[i]))
        if reverse_dip:
            write_xy = top_xy.copy()
            depth_col = z_nodes[i]
            dip_col = 180.0 + dips[i]
        else:
            write_xy = top_xy + uvec * hstep
            depth_col = z_nodes[i + 1]
            dip_col = dips[i]

        pts = np.vstack([write_xy, write_xy + B - A])
        plon, plat = proj.to_ll(pts[:, 0], pts[:, 1])
        la = plat[:nsegs]
        lb = plat[nsegs:]
        lo = plon[:nsegs]
        lob = plon[nsegs:]

        block = out[i * nsegs:(i + 1) * nsegs]
        block[:, 0] = widths[i]
        block[:, 1] = depth_col
        block[:, 2] = dip_col
        block[:, 3] = la
        block[:, 4] = lo
        block[:, 5] = lb
        block[:, 6] = lob
        block[:, 7] = args.ss
        block[:, 8] = args.ds
        block[:, 9] = args.op

        rows_xy.append(write_xy.copy())
        top_xy = top_xy + uvec * hstep

    qc = dict(proj=proj, A=A, B=B, uvec=uvec, strikes=strikes, dip_az=dip_az,
              nsegs=nsegs, total_len=total_len, rows_xy=rows_xy,
              reverse_dip=reverse_dip, az0=az0, reverse=reverse,
              chosen_az=chosen_az)
    return out, cfg, qc


def write_output(args, cfg, out, qc):
    n_out = cfg['n_out']
    nsegs = qc['nsegs']
    dips = cfg['dips']
    widths = cfg['widths']
    z_nodes = cfg['z_nodes']

    fmt = ' '.join(['%12.4f'] * 3 + ['%15.6f'] * 4 + ['%15.3f'] * 3)
    if args.dip_mode == 'reverse':
        mode_txt = 'trace is the shallow edge; dip column = 180 + angle'
    else:
        mode_txt = 'trace is the deep edge; dip column = true angle'
    header = [
        "faultgeom built from a surface trace by trace2faultgeom.py",
        "trace      : %s (%d patches along strike)" % (args.trace, nsegs),
        "grid       : nsegs=%d ndeps=%d" % (nsegs, n_out),
        "sizing     : %s" % ('width-controlled (uniform %.4f km patches)'
                             % np.mean(widths) if cfg['width_mode'] else
                             'depth-controlled'),
        "depth      : top=%.4f km bottom=%.4f km"
        % (z_nodes[0], z_nodes[-1]),
        "dip        : %.4f - %.4f deg (mean %.4f deg)"
        % (np.min(dips), np.max(dips), np.mean(dips)),
        "width      : %.4f - %.4f km (mean %.4f km, total %.4f km)"
        % (np.min(widths), np.max(widths), np.mean(widths),
           np.sum(widths)),
        "dip dir    : %.2f - %.2f deg (clockwise from north)"
        % (np.min(qc['dip_az']), np.max(qc['dip_az'])),
        "convention : %s" % mode_txt,
        "width_km depth_km dip_deg lat0_deg lon0_deg lat1_deg lon1_deg "
        "slp_strk_mm slp_ddip_mm slp_tens_mm",
    ]
    with open(args.output, 'w') as fid:
        if not args.no_header:
            for line in header:
                fid.write("# " + line + "\n")
        np.savetxt(fid, out, fmt=fmt)
    return header


def print_summary(args, cfg, qc):
    dips = cfg['dips']
    widths = cfg['widths']
    z_nodes = cfg['z_nodes']
    n_out = cfg['n_out']
    nsegs = qc['nsegs']
    dx = widths * np.cos(np.radians(dips))

    print("")
    print("trace      : %s" % args.trace)
    print("  %d segments along strike, total length %.3f km"
          % (nsegs, qc['total_len']))
    print("  along-strike direction A -> B : %.2f deg (mean)"
          % qc['az0'])
    if qc['reverse']:
        print("  trace written in reverse order (--dip-side left) : "
              "A -> B = %.2f deg" % ((qc['az0'] + 180.0) % 360.0))
    print("  down-dip azimuth              : %.2f - %.2f deg"
          % (np.min(qc['dip_az']), np.max(qc['dip_az'])))
    print("  dip direction relative to A->B: right-hand side")
    print("")
    print("geometry   : nsegs=%d ndeps=%d (%d patches)"
          % (nsegs, n_out, nsegs * n_out))
    print("  depth        : %.4f - %.4f km" % (z_nodes[0], z_nodes[-1]))
    print("  dip          : %.4f - %.4f deg (mean %.4f deg)"
          % (np.min(dips), np.max(dips), np.mean(dips)))
    print("  patch width  : %.4f - %.4f km (mean %.4f km, total %.4f km)"
          % (np.min(widths), np.max(widths), np.mean(widths),
             np.sum(widths)))
    print("  horizontal   : %.4f km over the full down-dip extent"
          % np.sum(dx))
    if cfg['width_mode']:
        print("  mode         : width-controlled (uniform patches)")
    else:
        print("  mode         : depth-controlled")
    print("  convention   : %s"
          % ('shallow edge given, dip = 180 + angle'
             if args.dip_mode == 'reverse' else
             'deep edge given, dip = true angle'))
    print("")
    print("  row   top(km)  bottom(km)   dip(deg)  width(km)   dx(km)")
    for i in range(n_out):
        print("  %3d %9.3f %10.3f %10.3f %10.3f %8.3f"
              % (i, z_nodes[i], z_nodes[i + 1], dips[i], widths[i], dx[i]))
    print("")
    print("output written to %s" % args.output)


# ---------------------------------------------------------------------------
# QC figure
# ---------------------------------------------------------------------------
SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK_2 = '#52514e'
GRID = '#d9d8d4'
C_BLUE = '#2a78d6'
C_ORANGE = '#eb6834'
BLUE_RAMP = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5',
             '#256abf', '#184f95', '#0d366b']


def make_plot(filename, args, cfg, qc):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    proj = qc['proj']
    nsegs = qc['nsegs']
    column = args.plot_column
    if column < 0 or column >= nsegs:
        raise ValueError("--plot-column must be in [0, %d)" % nsegs)

    dips = cfg['dips']
    widths = cfg['widths']
    z_nodes = cfg['z_nodes']
    n_out = cfg['n_out']
    dx = widths * np.cos(np.radians(dips))
    x_nodes = np.r_[0.0, np.cumsum(dx)]
    z_mid = 0.5 * (z_nodes[:-1] + z_nodes[1:])

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6),
                             facecolor=SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
        for spine in ax.spines.values():
            spine.set_color(GRID)
        ax.tick_params(colors=INK_2, labelsize=9)
        ax.grid(True, color=GRID, linestyle=':', linewidth=0.6)

    # ---- (1) map view ---------------------------------------------------
    ax = axes[0]
    lon0, lat0 = proj.to_ll(qc['A'][:, 0], qc['A'][:, 1])
    lon1, lat1 = proj.to_ll(qc['B'][:, 0], qc['B'][:, 1])
    tr_lon = np.r_[lon0, lon1[-1]]
    tr_lat = np.r_[lat0, lat1[-1]]
    ax.plot(tr_lon, tr_lat, color=C_BLUE, linewidth=2.0,
            label='surface trace')
    deep = qc['rows_xy'][-1]
    dlon, dlat = proj.to_ll(deep[:, 0], deep[:, 1])
    ax.plot(np.r_[dlon, dlon[0]], np.r_[dlat, dlat[0]], color=C_BLUE,
            linewidth=1.4, linestyle='--', label='deepest row')

    # dip-direction ticks at a few columns
    step = max(1, nsegs // 8)
    tick_len = 0.25 * max(1.0, np.sum(dx))
    for j in range(0, nsegs, step):
        e0 = 0.5 * (qc['A'][j, 0] + qc['B'][j, 0])
        n0 = 0.5 * (qc['A'][j, 1] + qc['B'][j, 1])
        e1 = e0 + qc['uvec'][j, 0] * tick_len
        n1 = n0 + qc['uvec'][j, 1] * tick_len
        tlon, tlat = proj.to_ll([e0, e1], [n0, n1])
        ax.plot(tlon, tlat, color=C_ORANGE, linewidth=1.2,
                solid_capstyle='round',
                label='dip direction' if j == 0 else None)

    ax.set_aspect(1.0 / max(0.2, np.cos(np.radians(np.mean(lat0)))))
    ax.set_xlabel('longitude (deg)', color=INK_2)
    ax.set_ylabel('latitude (deg)', color=INK_2)
    ax.set_title('Map view', color=INK, fontsize=12, loc='left')
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2,
              loc='best')

    # ---- (2) cross-section ---------------------------------------------
    ax = axes[1]
    ax.plot(x_nodes, z_nodes, color=C_BLUE, linewidth=2.0,
            marker='o', markersize=3.5, label='fault surface')
    ax.plot([0.0], [z_nodes[0]], marker='o', markersize=8,
            markerfacecolor=C_ORANGE, markeredgecolor=SURFACE,
            markeredgewidth=1.5, linestyle='none', label='surface trace')
    ax.invert_yaxis()
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlabel('horizontal distance along dip (km)', color=INK_2)
    ax.set_ylabel('depth (km)', color=INK_2)
    ax.set_title('Cross-section (column %d)' % column, color=INK,
                 fontsize=12, loc='left')
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc='best')
    if n_out <= 20:
        for i in range(n_out + 1):
            ax.annotate('%d' % i, (x_nodes[i], z_nodes[i]),
                        textcoords='offset points', xytext=(4, -8),
                        fontsize=7, color=INK_2)

    # ---- (3) dip profile -----------------------------------------------
    ax = axes[2]
    ax.plot(dips, z_mid, color=C_BLUE, linewidth=2.0, marker='o',
            markersize=4, label='dip profile')
    ax.invert_yaxis()
    ax.set_xlabel('dip (degrees)', color=INK_2)
    ax.set_ylabel('depth (km)', color=INK_2)
    ax.set_title('Dip versus depth', color=INK, fontsize=12, loc='left')
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc='best')

    fig.tight_layout()
    fig.savefig(filename, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3D patch view: rendering is delegated to pyginvc's Fault.plot_faultgeom
# ---------------------------------------------------------------------------
def make_plot3d(filename, args, cfg, qc, out):
    """3D view of the built patches with pyginvc's own plotter.

    Patch polygons, patch indices and the aftershock overlay are drawn by
    ``pyginvc.Geometry.Patch.Fault.plot_faultgeom``; this wrapper only picks
    the per-patch colour value and saves the resulting figure.
    """
    import matplotlib
    if not args.plot3d_show:
        matplotlib.use('Agg')
    import warnings

    import matplotlib.pyplot as plt
    from pyginvc.Geometry.Patch import Fault

    n = qc['nsegs']
    zmid = 0.5 * (cfg['z_nodes'][:-1] + cfg['z_nodes'][1:])
    if args.color_by == 'depth':
        value = np.repeat(zmid, n)
        clabel = 'depth (km)'
    elif args.color_by == 'dip':
        value = np.repeat(cfg['dips'], n)
        clabel = 'dip (deg)'
    elif args.color_by == 'width':
        value = np.repeat(cfg['widths'], n)
        clabel = 'patch width (km)'
    else:
        value = np.hypot(np.hypot(out[:, 7], out[:, 8]), out[:, 9])
        clabel = 'slip (mm)'
        if not np.any(value):
            print("WARNING: all slip values are zero; colouring by depth "
                  "instead.", file=sys.stderr)
            args.color_by = 'depth'
            value = np.repeat(zmid, n)
            clabel = 'depth (km)'

    fault = Fault(args.output, n, cfg['n_out'], False)
    fault.load_fault()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore', message='FigureCanvasAgg is non-interactive')
        fault.plot_faultgeom(value=value, afsfile=args.aftershock,
                             coordtype='llh', show_idx=args.show_index,
                             azimuth=args.view_azimuth,
                             elevation=args.view_elevation)

    # the built-in colourbar is hard-wired to the file's slip values, so
    # relabel/re-scale it when another quantity was used for the patches
    fig = plt.gcf()
    fig.set_size_inches(10.0, 8.0)
    if args.color_by != 'slip':
        cb = fig.axes[-1]
        ax3 = fig.axes[0]
        vmax = float(np.max(np.abs(value))) or 1.0
        cb.remove()
        sm = plt.cm.ScalarMappable(norm=plt.Normalize(0.0, vmax),
                                   cmap=plt.cm.jet)
        fig.colorbar(sm, ax=ax3, shrink=0.6, pad=0.1, label=clabel)

    fig.savefig(filename, dpi=160, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)


# ---------------------------------------------------------------------------
def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        cfg = resolve_arguments(args)
        out, cfg, qc = build_geometry(args, cfg)
        write_output(args, cfg, out, qc)
        print_summary(args, cfg, qc)
        if args.plot is not None:
            plotfile = args.plot if args.plot else 'trace2faultgeom.png'
            make_plot(plotfile, args, cfg, qc)
            print("QC figure written to %s" % plotfile)
        if args.plot3d is not None:
            plotfile3d = args.plot3d if args.plot3d else \
                os.path.splitext(args.output)[0] + '_3d.png'
            make_plot3d(plotfile3d, args, cfg, qc, out)
            print("3D figure written to %s" % plotfile3d)
    except (ValueError, IOError, ImportError) as err:
        print("ERROR: %s" % err, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
