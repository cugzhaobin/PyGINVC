#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
faultgeom2curvedfault_profile.py

Build a curved (listric) fault from a planar / gridded faultgeom file while
preserving either the top and bottom depth of the fault (default) or its
down-dip width (--keep-width).

Depth-parameterised construction (default)
------------------------------------------
    z     = depth nodes from z_top to z_bot            (nrows + 1 values)
    dip_i = dip of row i                               (profile or --dips)
    dz_i  = z[i+1] - z[i]
    width_i = dz_i / sin(dip_i)                        (down-dip width)
    dx_i    = dz_i / tan(dip_i)                        (horizontal step)

The fault surface is continuous: the bottom edge of row i is the top edge of
row i+1.  The map position of the top edge of the whole fault is kept fixed.

Width-preserving construction (--keep-width)
--------------------------------------------
    width_i = down-dip width of output row i (taken from the input and
              resampled linearly when the number of output rows changes)
    dx_i    = width_i * cos(dip_i)                     (horizontal step)
    z[0]    = top depth
    z[i+1]  = z[i] + width_i * sin(dip_i)

Here the top depth and the input down-dip widths are kept, so the bottom
depth follows from the widths and the dip values.  t in the dip profile is
then the normalised down-dip distance (cumulative width) instead of
normalised depth, and --t0_depth / --s_depth do not apply.

Dip profiles (t = normalised depth, 0 at the top, 1 at the bottom):

    linear   f(t) = 1 - t
    power    f(t) = 1 - t**p
             p > 1 : upper part stays steep, the main change is at depth
             p < 1 : upper part flattens quickly
    logistic h(t) = 1 / (1 + exp((t - t0) / s))
             f(t) = (h(t) - h(1)) / (h(0) - h(1))
             t0 = centre of the steep-to-gentle transition
             s  = width of the transition zone (smaller s = sharper bend)

    dip(t) = bottom_dip + (top_dip - bottom_dip) * f(t)

Alternatively, the true dip of every output row can be given directly:

    --dips 80,74,66,55,45,38,33,30

Input / output faultgeom format (same as the rest of PyGINVC):

    width_km depth_km dip_deg lat0 lon0 lat1 lon1 ss ds op

The input must be a regular nsegs x ndeps grid written row by row along
strike.  The input dip convention is preserved by default: if the input uses
180 + dip (for example 187 for a 7-degree dip), the output does too.  Use
--dip-mode angle to always write the true dip angle (0-90).

The projection and faultgeom reader are taken from the PyGINVC library
(pyginvc.libs.geotools and pyginvc.Geometry.Patch).  Dependencies: numpy,
pyproj (through geotools), plus the usual PyGINVC imports (matplotlib only
for --plot).
"""

import argparse
import sys

import numpy as np

from pyginvc.Geometry.Patch import Fault
from pyginvc.libs import geotools as gt


# ---------------------------------------------------------------------------
# parsing helpers
# ---------------------------------------------------------------------------
def load_faultgeom(filename):
    """Read a faultgeom file using the PyGINVC Fault geometry loader."""
    dat = Fault(filename, 1, 1, False).LoadFaultGeom()
    if dat.size == 0:
        raise IOError("cannot read fault geometry from %s" % filename)
    if dat.shape[1] < 10:
        raise ValueError("faultgeom file needs at least 10 columns, got %d"
                         % dat.shape[1])
    return dat[:, :10]


def cluster_depths(dep, tol=1.0e-3):
    """Group row indices by depth (tol in km), shallowest group first.

    Within one group the original file order (along strike) is preserved.
    """
    order = np.argsort(dep, kind='stable')
    d = dep[order]
    cuts = np.where(np.diff(d) > tol)[0] + 1
    bounds = np.r_[0, cuts, len(d)]
    return [order[bounds[k]:bounds[k + 1]] for k in range(len(bounds) - 1)]


def true_dip(raw_dip):
    """Fold a dip value into (0, 90] degrees."""
    d = np.mod(raw_dip, 360.0)
    d = np.where(d > 180.0, d - 180.0, d)
    d = np.where(d > 90.0, 180.0 - d, d)
    return d


def parse_floats(text):
    return np.array([float(v) for v in text.replace(' ', '').split(',')
                     if v != ''], dtype=float)


def dip_profile(t, top_dip, bottom_dip, mode, power, t0, s):
    """Return dip(t) for t in [0, 1]."""
    t = np.asarray(t, dtype=float)
    if mode == 'linear':
        f = 1.0 - t
    elif mode == 'power':
        if power <= 0.0:
            raise ValueError("--power must be > 0")
        f = 1.0 - np.power(t, power)
    elif mode == 'logistic':
        if s <= 0.0:
            raise ValueError("--s / --s-depth must be > 0")
        if not (0.0 < t0 < 1.0):
            raise ValueError("--t0 / --t0-depth must be inside the fault "
                             "(0 < t0 < 1)")
        h0 = 1.0 / (1.0 + np.exp(-t0 / s))
        h1 = 1.0 / (1.0 + np.exp((1.0 - t0) / s))
        h = 1.0 / (1.0 + np.exp((t - t0) / s))
        f = (h - h1) / (h0 - h1)
    else:
        raise ValueError("unknown profile: %s" % mode)
    return bottom_dip + (top_dip - bottom_dip) * f


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Build a curved fault from a planar faultgeom file.  By "
                    "default the top and bottom depth are preserved; with "
                    "--keep-width the down-dip width is preserved instead.")

    parser.add_argument('--faultfile', type=str, required=True,
                        help='input faultgeom file (width depth dip lat0 lon0 '
                             'lat1 lon1 ss ds op)')
    parser.add_argument('--output', type=str, default='curved.faultgeom',
                        help='output faultgeom file (default: curved.faultgeom)')

    # input grid size (optional, inferred from the depth column if omitted)
    parser.add_argument('--nsegs', type=int, default=0,
                        help='number of patches along strike (default: infer)')
    parser.add_argument('--ndeps', type=int, default=0,
                        help='number of input rows along dip (default: infer)')

    # output resolution
    parser.add_argument('--out_ndeps', type=int, default=0,
                        help='number of output rows along dip (default: same '
                             'as input; 10 for a single-row input)')

    # dip description
    parser.add_argument('--top_dip', type=float, default=None,
                        help='dip at the top of the fault (degrees)')
    parser.add_argument('--bottom_dip', type=float, default=None,
                        help='dip at the bottom of the fault (degrees)')
    parser.add_argument('--profile', type=str, default='power',
                        choices=['linear', 'power', 'logistic'],
                        help='dip profile (default: power)')
    parser.add_argument('--power', type=float, default=1.0,
                        help='exponent p of the power profile (default: 1.0)')
    parser.add_argument('--t0', type=float, default=0.5,
                        help='logistic transition centre, normalised depth '
                             '(normalised down-dip distance with --keep-width) '
                             '(default: 0.5)')
    parser.add_argument('--s', type=float, default=0.1,
                        help='logistic transition width, normalised depth '
                             '(normalised down-dip distance with --keep-width) '
                             '(default: 0.1)')
    parser.add_argument('--t0_depth', type=float, default=None,
                        help='logistic transition centre in km (overrides '
                             '--t0; not available with --keep-width)')
    parser.add_argument('--s_depth', type=float, default=None,
                        help='logistic transition width in km (overrides --s; '
                             'not available with --keep-width)')
    parser.add_argument('--dips', type=str, default=None,
                        help='explicit true dip for every output row, '
                             'comma-separated, top to bottom; overrides the '
                             'profile')

    # depth description
    parser.add_argument('--top_depth', type=float, default=None,
                        help='override the top depth read from the input (km)')
    parser.add_argument('--bottom_depth', type=float, default=None,
                        help='override the bottom depth read from the input '
                             '(km; not available with --keep-width)')
    parser.add_argument('--depth_nodes', type=str, default=None,
                        help='explicit depth nodes (km), comma-separated; '
                             'needs out_ndeps + 1 values and overrides '
                             'top_depth / bottom_depth (not available with '
                             '--keep-width)')
    parser.add_argument('--keep-width', action='store_true',
                        help='keep the input down-dip width of every row '
                             'instead of the top/bottom depth; the bottom '
                             'depth then follows from the widths and dips')

    # conventions / extras
    parser.add_argument('--dip-mode', type=str, default='preserve',
                        choices=['preserve', 'angle', 'reverse'],
                        help='output dip convention: preserve the input '
                             'convention (default), always write the true '
                             'angle (angle), or always write 180 + angle '
                             '(reverse)')
    parser.add_argument('--slip', type=str, default='copy',
                        choices=['copy', 'zero'],
                        help='copy the input slip columns when the number of '
                             'rows is unchanged, otherwise write zeros '
                             '(default: copy)')
    parser.add_argument('--flip-dip-direction', action='store_true',
                        help='flip the inferred down-dip direction (use this '
                             'if the curve bends to the wrong side for a '
                             'single-row input)')
    parser.add_argument('--plot', type=str, nargs='?', const='',
                        default=None,
                        help='save a QC cross-section / dip-profile figure '
                             '(default name: curvedfault_profile.png)')
    parser.add_argument('--plot-column', type=int, default=0,
                        help='along-strike column used for the QC figure '
                             '(default: 0)')

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. read input and find the regular nsegs x ndeps grid
    # ------------------------------------------------------------------
    dat = load_faultgeom(args.faultfile)
    nrow = dat.shape[0]

    w_in = dat[:, 0].copy()
    dep_in = dat[:, 1].copy()
    raw_dip = dat[:, 2].copy()
    lat0, lon0 = dat[:, 3].copy(), dat[:, 4].copy()
    lat1, lon1 = dat[:, 5].copy(), dat[:, 6].copy()
    slip_in = dat[:, 7:10].copy()

    groups = cluster_depths(dep_in)
    in_ndeps = len(groups)
    in_nsegs = len(groups[0])
    if any(len(g) != in_nsegs for g in groups):
        raise ValueError("input is not a regular grid: depth groups have "
                         "different lengths")
    if in_nsegs * in_ndeps != nrow:
        raise ValueError("internal error while grouping the input rows")
    if args.nsegs and args.nsegs != in_nsegs:
        raise ValueError("--nsegs=%d but the input has %d patches per row"
                         % (args.nsegs, in_nsegs))
    if args.ndeps and args.ndeps != in_ndeps:
        raise ValueError("--ndeps=%d but the input has %d depth rows"
                         % (args.ndeps, in_ndeps))
    nsegs = in_nsegs

    # ------------------------------------------------------------------
    # 2. local projection and patch geometry
    # ------------------------------------------------------------------
    origin_lat = float(np.mean(np.r_[lat0, lat1]))
    origin_lon = float(np.mean(np.r_[lon0, lon1]))
    origin = np.array([origin_lat, origin_lon])

    e0 = np.zeros(nrow)
    n0 = np.zeros(nrow)
    e1 = np.zeros(nrow)
    n1 = np.zeros(nrow)
    for k in range(nrow):
        e0[k], n0[k] = gt.llh2utm([lat0[k], lon0[k]], origin)
        e1[k], n1[k] = gt.llh2utm([lat1[k], lon1[k]], origin)

    theta = np.deg2rad(true_dip(raw_dip))

    # depth of the two edges of each input patch.
    # 'dep_in' is the depth of the edge whose coordinates are given; the other
    # edge follows the same rule used by Geometry/Fault.py.
    z_given = dep_in
    z_other = dep_in - w_in * np.sin(np.deg2rad(raw_dip))
    z_top_patch = np.minimum(z_given, z_other)
    z_bot_patch = np.maximum(z_given, z_other)

    z_top_in = float(np.min(z_top_patch))
    z_bot_in = float(np.max(z_bot_patch))

    # check that the shallow and deep edges are roughly flat along strike
    top_row = z_top_patch[groups[0]]
    bot_row = z_bot_patch[groups[-1]]
    if np.ptp(top_row) > 0.1:
        print("WARNING: the shallowest edge varies by %.3f km along strike; "
              "a single top depth is used." % np.ptp(top_row),
              file=sys.stderr)
    if np.ptp(bot_row) > 0.1:
        print("WARNING: the deepest edge varies by %.3f km along strike; "
              "a single bottom depth is used." % np.ptp(bot_row),
              file=sys.stderr)

    # centres of the input patches
    centres = np.zeros((in_ndeps, nsegs, 2))
    for i, g in enumerate(groups):
        for j, k in enumerate(g):
            centres[i, j, 0] = 0.5 * (e0[k] + e1[k])
            centres[i, j, 1] = 0.5 * (n0[k] + n1[k])

    # ------------------------------------------------------------------
    # 3. infer the horizontal down-dip unit vector for every column
    # ------------------------------------------------------------------
    u = np.zeros((nsegs, 2))
    good = np.zeros(nsegs, dtype=bool)
    if in_ndeps >= 2:
        shifts = np.diff(centres, axis=0)          # (ndeps-1, nsegs, 2)
        v = np.mean(shifts, axis=0)                # (nsegs, 2)
        norm = np.linalg.norm(v, axis=1)
        good = norm > 1.0e-9
        u[good] = v[good] / norm[good, None]

    # fallback for a single-row input (or a degenerate shift): use the
    # along-strike direction rotated clockwise by 90 degrees.
    for j in range(nsegs):
        if not good[j]:
            k = groups[0][j]
            sx = e1[k] - e0[k]
            sy = n1[k] - n0[k]
            sn = np.hypot(sx, sy)
            if sn < 1.0e-9:
                raise ValueError("degenerate along-strike length in column %d"
                                 % j)
            u[j, 0] = sy / sn
            u[j, 1] = -sx / sn
            good[j] = True

    if args.flip_dip_direction:
        u = -u

    # consistency check on the inferred direction
    if in_ndeps >= 2 and nsegs > 1:
        bearings = np.degrees(np.arctan2(u[:, 0], u[:, 1]))
        ang = np.angle(np.mean(np.exp(1j * np.radians(bearings))))
        spread = np.degrees(np.std(np.angle(
            np.exp(1j * (np.radians(bearings) - ang)))))
        if spread > 5.0:
            print("WARNING: the inferred down-dip direction varies by %.1f "
                  "degrees along strike; per-column directions are used."
                  % spread, file=sys.stderr)

    # ------------------------------------------------------------------
    # 4. top trace of the fault (map coordinates and depth)
    # ------------------------------------------------------------------
    top_l = np.zeros((nsegs, 2))
    top_r = np.zeros((nsegs, 2))
    for j, k in enumerate(groups[0]):
        if z_given[k] <= z_other[k]:
            # the given coordinates belong to the shallow edge
            top_l[j] = [e0[k], n0[k]]
            top_r[j] = [e1[k], n1[k]]
        else:
            # the given coordinates belong to the deep edge: move up-dip
            off = w_in[k] * np.cos(theta[k])
            top_l[j] = [e0[k] - u[j, 0] * off, n0[k] - u[j, 1] * off]
            top_r[j] = [e1[k] - u[j, 0] * off, n1[k] - u[j, 1] * off]

    # ------------------------------------------------------------------
    # 5. output depth nodes and dip values
    # ------------------------------------------------------------------
    if args.keep_width and args.depth_nodes:
        raise ValueError("--keep-width cannot be combined with --depth_nodes")
    if args.keep_width and args.bottom_depth is not None:
        raise ValueError("--keep-width cannot be combined with "
                         "--bottom_depth: the bottom depth follows from the "
                         "preserved width and the dip values")
    if args.keep_width and (args.t0_depth is not None or
                            args.s_depth is not None):
        raise ValueError("--keep-width parameterises the dip profile by the "
                         "normalised down-dip distance, so --t0_depth and "
                         "--s_depth do not apply")

    if args.depth_nodes:
        z_nodes = parse_floats(args.depth_nodes)
        if len(z_nodes) < 2:
            raise ValueError("--depth_nodes needs at least two values")
        if np.any(np.diff(z_nodes) <= 0.0):
            raise ValueError("--depth_nodes must be strictly increasing")
        out_from_nodes = len(z_nodes) - 1
    else:
        z_nodes = None
        out_from_nodes = None

    if args.dips:
        dips = parse_floats(args.dips)
        out_from_dips = len(dips)
    else:
        dips = None
        out_from_dips = None

    candidates = [v for v in (args.out_ndeps, out_from_dips,
                              out_from_nodes) if v]
    if len(set(candidates)) > 1:
        raise ValueError("--out_ndeps, --dips and --depth_nodes disagree "
                         "about the number of output rows: %s" % candidates)
    if candidates:
        out_ndeps = candidates[0]
    else:
        out_ndeps = in_ndeps if in_ndeps > 1 else 10

    if out_ndeps < 2:
        raise ValueError("need at least two output rows along dip")

    # representative input width and dip of every input row
    w_row = np.array([np.mean(w_in[g]) for g in groups])
    theta_row = np.array([np.mean(theta[g]) for g in groups])

    if args.keep_width:
        z_top = z_top_in if args.top_depth is None else float(args.top_depth)

        # keep the input down-dip widths.  If the number of output rows
        # changes, resample the cumulative width profile linearly: this
        # preserves the total width and reduces to the input widths when the
        # row count is unchanged.
        for i, g in enumerate(groups):
            if np.ptp(w_in[g]) > 0.1:
                print("WARNING: the down-dip width of row %d varies by "
                      "%.3f km along strike; its mean value is used."
                      % (i, np.ptp(w_in[g])), file=sys.stderr)
        if out_ndeps == in_ndeps:
            width = w_row.copy()
        else:
            print("NOTE: the input has %d dip row(s); the preserved width is "
                  "redistributed over %d output rows."
                  % (in_ndeps, out_ndeps), file=sys.stderr)
            cum = np.r_[0.0, np.cumsum(w_row)]
            u_nodes = np.linspace(0.0, in_ndeps, out_ndeps + 1)
            width = np.diff(np.interp(u_nodes,
                                      np.arange(in_ndeps + 1), cum))

        # t is the normalised down-dip distance; unlike the depth mode this
        # is known before the dips because the widths come from the input
        s_mid = np.cumsum(width) - 0.5 * width
        t_norm = s_mid / np.sum(width)
        t0, s = args.t0, args.s
    else:
        if z_nodes is None:
            z_top = z_top_in if args.top_depth is None else float(args.top_depth)
            z_bot = z_bot_in if args.bottom_depth is None else float(args.bottom_depth)
            if z_bot <= z_top:
                raise ValueError("bottom depth (%.3f) must be larger than top "
                                 "depth (%.3f)" % (z_bot, z_top))
            z_nodes = np.linspace(z_top, z_bot, out_ndeps + 1)
        else:
            z_top, z_bot = float(z_nodes[0]), float(z_nodes[-1])

        dz = np.diff(z_nodes)
        t_mid = 0.5 * (z_nodes[:-1] + z_nodes[1:])
        t_norm = (t_mid - z_top) / (z_bot - z_top)
        t0, s = args.t0, args.s
        if args.t0_depth is not None:
            t0 = (args.t0_depth - z_top) / (z_bot - z_top)
        if args.s_depth is not None:
            s = args.s_depth / (z_bot - z_top)

    if dips is None:
        if args.top_dip is None or args.bottom_dip is None:
            raise ValueError("--top_dip and --bottom_dip are required unless "
                             "--dips is used")
        dips = dip_profile(t_norm, args.top_dip, args.bottom_dip,
                           args.profile, args.power, t0, s)
    elif len(dips) != out_ndeps:
        raise ValueError("--dips has %d values but %d output rows are "
                         "requested" % (len(dips), out_ndeps))

    dips = np.asarray(dips, dtype=float)
    if np.any(dips <= 0.0) or np.any(dips > 90.0):
        raise ValueError("dip values must satisfy 0 < dip <= 90 degrees")

    theta_out = np.deg2rad(dips)

    if args.keep_width:
        hstep = width * np.cos(theta_out)
        z_nodes = z_top + np.r_[0.0, np.cumsum(width * np.sin(theta_out))]
        z_bot = float(z_nodes[-1])
    else:
        width = dz / np.sin(theta_out)
        # dz / tan(dip) = dz * cos(dip) / sin(dip); 0 for a vertical fault
        hstep = dz * np.cos(theta_out) / np.sin(theta_out)
        hstep = np.where(np.abs(np.cos(theta_out)) < 1.0e-12, 0.0, hstep)

    # ------------------------------------------------------------------
    # 6. output dip convention
    # ------------------------------------------------------------------
    input_reversed = np.mean(raw_dip > 180.0) > 0.5
    if np.any(raw_dip > 180.0) and not input_reversed:
        print("WARNING: the input mixes dip conventions; the majority "
              "convention is used.", file=sys.stderr)

    if args.dip_mode == 'angle':
        out_reversed = False
    elif args.dip_mode == 'reverse':
        out_reversed = True
    else:
        out_reversed = input_reversed

    # ------------------------------------------------------------------
    # 7. build the output grid
    # ------------------------------------------------------------------
    out = np.zeros((out_ndeps * nsegs, 10))

    # slip is copied only when the row structure is unchanged
    copy_slip = (args.slip == 'copy') and (out_ndeps == in_ndeps)
    if args.slip == 'copy' and not copy_slip:
        print("NOTE: the number of rows changes, so the input slip is not "
              "copied; zeros are written.", file=sys.stderr)
    if copy_slip and np.any(np.abs(slip_in) > 0.0):
        print("NOTE: non-zero input slip is copied row by row.",
              file=sys.stderr)

    top_cur_l = top_l.copy()
    top_cur_r = top_r.copy()

    for i in range(out_ndeps):
        for j in range(nsegs):
            tl = top_cur_l[j]
            tr = top_cur_r[j]
            bl = tl + u[j] * hstep[i]
            br = tr + u[j] * hstep[i]

            if out_reversed:
                # the given edge is the shallow edge
                p0, p1 = tl, tr
                dep_out = z_nodes[i]
                dip_out = 180.0 + dips[i]
            else:
                # the given edge is the deep edge
                p0, p1 = bl, br
                dep_out = z_nodes[i + 1]
                dip_out = dips[i]

            # gt.utm2llh expects [north, east] and returns [lat, lon]
            la0, lo0 = gt.utm2llh([[p0[1], p0[0]]], origin)[0]
            la1, lo1 = gt.utm2llh([[p1[1], p1[0]]], origin)[0]

            if copy_slip:
                k = groups[i][j]
                ss, ds, op = slip_in[k]
            else:
                ss, ds, op = 0.0, 0.0, 0.0

            out[i * nsegs + j, :] = np.array([
                width[i], dep_out, dip_out,
                la0, lo0, la1, lo1, ss, ds, op
            ])

            top_cur_l[j] = bl
            top_cur_r[j] = br

    # ------------------------------------------------------------------
    # 8. write the output file
    # ------------------------------------------------------------------
    fmt = ' '.join(['%12.4f'] * 3 + ['%15.6f'] * 4 + ['%15.6f'] * 3)

    header_lines = [
        "curved fault generated by faultgeom2curvedfault_profile.py",
        "input      : %s" % args.faultfile,
        "grid       : nsegs=%d ndeps_out=%d" % (nsegs, out_ndeps),
        "mode       : %s" % ("keep down-dip width"
                             if args.keep_width else "keep top/bottom depth"),
        "depth      : top=%.4f km bottom=%.4f km" % (z_top, z_bot),
        "dip source : %s" % (("dips=" + args.dips) if args.dips else
                             ("profile=%s" % args.profile)),
        "dip mode   : %s" % ("180+angle" if out_reversed else "angle"),
        "width_km depth_km dip_deg lat0 lon0 lat1 lon1 ss ds op",
    ]
    with open(args.output, 'w') as fid:
        for line in header_lines:
            fid.write("# " + line + "\n")
        np.savetxt(fid, out, fmt=fmt)

    # ------------------------------------------------------------------
    # 9. summary
    # ------------------------------------------------------------------
    in_width = float(np.sum(w_row))
    in_dx = float(np.sum(w_row * np.cos(theta_row)))

    print("")
    print("input : nsegs=%d ndeps=%d" % (nsegs, in_ndeps))
    print("        depth column       = %.3f - %.3f km"
          % (np.min(dep_in), np.max(dep_in)))
    print("        true edges         = %.3f (top) - %.3f (bottom) km"
          % (z_top_in, z_bot_in))
    print("        total down-dip width = %.2f km, horizontal extent = %.2f km"
          % (in_width, in_dx))
    dep_out_col = z_nodes[:-1] if out_reversed else z_nodes[1:]
    print("output: nsegs=%d ndeps=%d" % (nsegs, out_ndeps))
    print("        depth column       = %.3f - %.3f km"
          % (np.min(dep_out_col), np.max(dep_out_col)))
    print("        true edges         = %.3f (top) - %.3f (bottom) km"
          % (z_nodes[0], z_nodes[-1]))
    print("        total down-dip width = %.2f km, horizontal extent = %.2f km"
          % (np.sum(width), np.sum(hstep)))
    print("        dip convention      = %s"
          % ("180 + angle" if out_reversed else "true angle"))
    print("        construction        = %s"
          % ("keep width (bottom depth computed)"
             if args.keep_width else "keep depth (width computed)"))
    out_width = float(np.sum(width))
    if args.keep_width:
        if abs(z_bot - z_bot_in) > 0.2 * max(abs(z_bot_in), 1.0e-9):
            print("NOTE: --keep-width preserves the down-dip width, so the "
                  "bottom depth changes from %.2f to %.2f km."
                  % (z_bot_in, z_bot))
    elif in_width > 0.0 and abs(out_width - in_width) / in_width > 0.2:
        print("NOTE: depth mode keeps the top/bottom depth fixed, so the "
              "down-dip width changes from %.2f to %.2f km (%.0f%%)."
              % (in_width, out_width,
                 100.0 * (out_width - in_width) / in_width))
    print("")
    print(" row    top(km)   bottom(km)   dip(deg)   width(km)    dx(km)")
    for i in range(out_ndeps):
        print("%4d %10.3f %11.3f %10.3f %11.3f %9.3f"
              % (i, z_nodes[i], z_nodes[i + 1], dips[i], width[i], hstep[i]))
    print("")
    print("output written to %s" % args.output)

    # ------------------------------------------------------------------
    # 10. optional QC figure
    # ------------------------------------------------------------------
    if args.plot is not None:
        if args.plot_column < 0 or args.plot_column >= nsegs:
            raise ValueError("--plot-column must be in [0, %d)" % nsegs)
        plotfile = args.plot if args.plot else 'curvedfault_profile.png'
        make_plot(plotfile, args.plot_column, groups, w_in, theta,
                  z_top_in, in_ndeps, dips, z_nodes, hstep, width,
                  out_ndeps)
        print("QC figure written to %s" % plotfile)


def make_plot(filename, column, groups, w_in, theta_in, z_top_in, in_ndeps,
              dips, z_nodes, hstep, width, out_ndeps):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 6.0))

    # --- cross-section -----------------------------------------------------
    ax = axes[0]
    x_in = [0.0]
    z_in = [z_top_in]
    for i in range(in_ndeps):
        k = groups[i][column]
        x_in.append(x_in[-1] + w_in[k] * np.cos(theta_in[k]))
        z_in.append(z_in[-1] + w_in[k] * np.sin(theta_in[k]))
    ax.plot(x_in, z_in, 'k--o', ms=3, lw=1.0, label='input')

    x_out = [0.0]
    for h in hstep:
        x_out.append(x_out[-1] + h)
    ax.plot(x_out, z_nodes, 'r-o', ms=4, lw=1.6, label='curved (output)')
    ax.invert_yaxis()
    ax.set_xlabel('horizontal distance along dip (km)')
    ax.set_ylabel('depth (km)')
    ax.set_title('fault cross-section (column %d)' % column)
    ax.grid(True, ls=':', lw=0.5)
    ax.legend(loc='best')

    # --- dip profile -------------------------------------------------------
    ax2 = axes[1]
    dip_in = [theta_in[groups[i][column]] * 180.0 / np.pi
              for i in range(in_ndeps)]
    # cumulative mid-depth of the input rows
    z_in_mid = []
    zc = z_top_in
    for i in range(in_ndeps):
        k = groups[i][column]
        z_in_mid.append(zc + 0.5 * w_in[k] * np.sin(theta_in[k]))
        zc += w_in[k] * np.sin(theta_in[k])

    ax2.plot(dip_in, z_in_mid, 'k--o', ms=3, lw=1.0, label='input')
    ax2.plot(dips, 0.5 * (z_nodes[:-1] + z_nodes[1:]), 'r-o', ms=4, lw=1.6,
             label='curved (output)')
    ax2.invert_yaxis()
    ax2.set_xlabel('dip (degrees)')
    ax2.set_ylabel('depth (km)')
    ax2.set_title('dip profile')
    ax2.grid(True, ls=':', lw=0.5)
    ax2.legend(loc='best')

    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    main()
