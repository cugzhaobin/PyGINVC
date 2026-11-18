#!/usr/bin/env python3
"""
filter_near_fault_points.py -- remove near-fault observations whose sign is
opposite to the model prediction.

Idea
----
Because the fault trace is simplified, a pixel that really belongs to one
wall of the fault can fall on the other wall in the model.  The predicted LOS
then has the wrong sign.  Such points are easiest to recognise directly from
the data: look only at points within a narrow band on both sides of the
surface trace and keep (i.e. delete) the points where

        sign(observation) != sign(model prediction)

and both amplitudes are above the noise.

The surface trace is read from the .faultgeom file itself -- no separate
detailed trace (e.g. GaoFen_Fault7.txt) is used.  For a PyGINVC fault file the
first ``ndeps`` rows are the surface patches, and their shallow edges
(SLon,SLat)->(ELon,ELat) form the trace.  A combined file with several faults
is split into its depth blocks automatically.  For a triangular fault give
``--solution solutions.h5`` instead: the .tri slip files store only node
indices, so the trace is built from ``flt/vertex_llh`` and ``flt/element``
(the boundary edge chain whose nodes sit at the surface).

Inputs (default names, all relative to --datadir)
------------------------------------------------
sar_obs_<track>.xyz          observed LOS (mm)
sm000/sar_mod_<track>.xyz    predicted LOS (mm), ramp already included
sm000/sar_res_<track>.xyz    residual (mm), report only
gps1_obs.gmtvec              GPS east/north (mm)
gps1_obs_up.gmtvec           GPS up (mm)
sm000/gps1_mod.gmtvec        modelled east/north (mm)
sm000/gps1_mod_up.gmtvec     modelled up (mm)

For a GPS station the same sign test is applied to E, N and Up; a station is
removed when any component flips.

Outputs (in --outdir, default clean_obs/)
-----------------------------------------
sar_obs_<track>.xyz          cleaned observation files
gps1_obs.gmtvec
gps1_obs_up.gmtvec
<head of --los>              cleaned source .los, if --los is given
removed_points.csv           every removed point with its geometry
near_fault_qc.png            map + obs/model scatter
filter_summary.txt           summary report

Example
-------
python filter_near_fault_points.py \
    --fault ../../../fault/trace_tdip_58_bdip30_decay10.faultgeom \
    --fault ../../../fault/seg2.faultgeom \
    --los ../../../data/insar/WX_downsample/T112A_los_lltnde.los
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KM_PER_DEG = 111.32


# --------------------------------------------------------------------------
# fault geometry
# --------------------------------------------------------------------------
def _read_header_int(path, key):
    with open(path, "r", errors="ignore") as f:
        head = f.read(4096)
    m = re.search(key + r"\s*=\s*(\d+)", head)
    return int(m.group(1)) if m else None


def read_fault_traces(path):
    """Surface traces of a PyGINVC .faultgeom file.

    The rows may hold one fault (ndeps blocks of nsegs patches) or several
    faults concatenated (PyGINVC's combined ``slipmodel_001.faultgeom``);
    a new fault block is detected where the row depth drops back.  Each
    fault contributes its surface row: the shallow edges
    (lon0,lat0)->(lon1,lat1) of the first nsegs rows.  Returns a list of
    (lon, lat, info) tuples.
    """
    nsegs_h = _read_header_int(path, "nsegs")
    ndeps_h = _read_header_int(path, "ndeps")

    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 7:
        raise ValueError("%s: not a faultgeom table" % path)

    dep = np.round(data[:, 1], 4)
    cuts = [0] + list(np.where(np.diff(dep) < -1e-6)[0] + 1) + [len(data)]

    traces = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        blk = data[a:b]
        nsegs = nsegs_h
        if not nsegs or len(blk) % nsegs:
            nsegs = int(np.sum(np.isclose(np.round(blk[:, 1], 4),
                                          np.round(blk[0, 1], 4))))
        if nsegs <= 0 or len(blk) % nsegs:
            raise ValueError("%s: could not determine the patch grid"
                             % path)
        ndeps = len(blk) // nsegs
        if ndeps_h and ndeps_h * nsegs == len(blk):
            ndeps = ndeps_h

        surf = blk[:nsegs]
        lon = np.r_[surf[:, 4], surf[-1, 6]]
        lat = np.r_[surf[:, 3], surf[-1, 5]]

        # continuity: end of patch i must equal start of patch i+1
        coslat = np.cos(np.radians(np.mean(lat)))
        jump = np.hypot((surf[:-1, 6] - surf[1:, 4]) * coslat,
                        surf[:-1, 5] - surf[1:, 3]) * KM_PER_DEG
        gap = float(jump.max()) if len(jump) else 0.0
        if gap > 0.05:
            print("WARNING: %s: surface row is discontinuous by %.1f km; "
                  "check row order" % (path, gap))

        info = dict(path=path, nsegs=nsegs, ndeps=ndeps, nnodes=len(lon),
                    jump_max_km=gap)
        traces.append((lon, lat, info))
    return traces


def read_fault_trace(path):
    """First surface trace of a .faultgeom file (single-fault files)."""
    return read_fault_traces(path)[0]


def read_triangle_traces(solutionfile):
    """Surface traces of a triangular fault from a solutions.h5 file.

    The ``*.tri`` slip-model files only store node indices and slip, so the
    vertex coordinates are taken from the solution file (``flt/vertex_llh``,
    columns lon/lat/depth, and ``flt/element``, 1-based triangles).  The
    surface trace is the boundary edge chain whose two end nodes sit at the
    minimum depth; several disconnected chains yield several traces.
    Returns a list of (lon, lat, info) tuples.
    """
    import h5py
    from collections import Counter, defaultdict

    with h5py.File(solutionfile, "r") as h5:
        flttype = h5["flt/flttype"][()]
        if isinstance(flttype, bytes):
            flttype = flttype.decode()
        if flttype != "triangle":
            raise ValueError("%s: fault type is %r, not triangle"
                             % (solutionfile, flttype))
        llh = np.asarray(h5["flt/vertex_llh"][()], float)   # lat lon up
        elem = np.asarray(h5["flt/element"][()], dtype=int) - 1

    edges = Counter()
    for tri in elem:
        for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            edges[(min(u, v), max(u, v))] += 1

    # vertex_llh[:, 2] is the up coordinate (positive up, 0 at the surface),
    # so the surface nodes sit at the MAXIMUM value.
    zsurf = llh[:, 2].max()
    surf = [e for e, n in edges.items() if n == 1
            and abs(llh[e[0], 2] - zsurf) < 1e-3
            and abs(llh[e[1], 2] - zsurf) < 1e-3]
    if not surf:
        raise ValueError("%s: no surface boundary edges found"
                         % solutionfile)

    # order the boundary edges into open/closed chains
    adj = defaultdict(list)
    for u, v in surf:
        adj[u].append(v)
        adj[v].append(u)

    seen = set()
    chains = []
    ends = [n for n, nb in adj.items() if len(nb) == 1]
    starts = ends + [n for n in adj if n not in ends]
    for start in starts:
        if start in seen:
            continue
        chain = [start]
        seen.add(start)
        cur, prev = start, None
        while True:
            nxt = [n for n in adj[cur] if n != prev and n not in seen]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            chain.append(cur)
            seen.add(cur)
        if len(chain) > 1:
            chains.append(chain)

    traces = []
    for chain in chains:
        lat = llh[chain, 0]
        lon = llh[chain, 1]
        coslat = np.cos(np.radians(np.mean(lat)))
        jump = np.hypot(np.diff(lon) * coslat, np.diff(lat)) * KM_PER_DEG
        info = dict(path=solutionfile, nsegs=len(chain) - 1, ndeps="tri",
                    nnodes=len(chain),
                    jump_max_km=float(jump.max()) if len(jump) else 0.0)
        traces.append((lon, lat, info))
    return traces


class FaultTrace:
    """Nearest-segment distance and hanging/foot-wall side of a polyline."""

    def __init__(self, path, name=None, block=0):
        self.path = path
        traces = read_fault_traces(path)
        self.lon, self.lat, self.info = traces[block]
        if name is None:
            base = os.path.splitext(os.path.basename(path))[0]
            name = base + ('_f%d' % (block + 1) if len(traces) > 1 else '')
        self.name = name
        self._setup()

    @classmethod
    def from_polyline(cls, lon, lat, info, name, path=""):
        """Build a trace from explicit polyline arrays (triangular faults)."""
        self = cls.__new__(cls)
        self.path = path
        self.name = name
        self.lon = np.asarray(lon, float)
        self.lat = np.asarray(lat, float)
        self.info = info
        self._setup()
        return self

    def _setup(self):
        self.lat0 = np.radians(np.mean(self.lat))

        self.lat0 = np.radians(np.mean(self.lat))
        self.kx = KM_PER_DEG * np.cos(self.lat0)
        self.ky = KM_PER_DEG

        self.x = self.lon * self.kx
        self.y = self.lat * self.ky

        # dip azimuth = strike + 90 deg (pyginvc convention; the dip column is
        # 180+angle, so it cannot be used directly as an azimuth)
        dx = np.diff(self.x)
        dy = np.diff(self.y)
        strike = (90.0 - np.degrees(np.arctan2(dy, dx))) % 360.0
        self.dipaz = (strike + 90.0) % 360.0
        az = np.radians(self.dipaz)
        self.ux = np.sin(az)
        self.uy = np.cos(az)

    def nearest(self, lon, lat):
        """Per point: distance (km), nearest segment, side (+1 HW, -1 FW)."""
        qx = np.asarray(lon) * self.kx
        qy = np.asarray(lat) * self.ky
        n = len(qx)
        dist = np.full(n, np.inf)
        seg = np.full(n, -1, dtype=int)
        side = np.zeros(n, dtype=int)

        for i in range(len(self.x) - 1):
            vx = self.x[i + 1] - self.x[i]
            vy = self.y[i + 1] - self.y[i]
            L2 = vx * vx + vy * vy
            if L2 == 0.0:
                continue
            t = np.clip(((qx - self.x[i]) * vx + (qy - self.y[i]) * vy) / L2,
                        0.0, 1.0)
            cx = self.x[i] + t * vx
            cy = self.y[i] + t * vy
            d = np.hypot(qx - cx, qy - cy)
            m = d < dist
            if m.any():
                dist[m] = d[m]
                seg[m] = i
                side[m] = np.sign((qx[m] - cx[m]) * self.ux[i]
                                  + (qy[m] - cy[m]) * self.uy[i])
        return dist, seg, side


def nearest_fault(traces, lon, lat):
    """Nearest trace among several; returns dist, trace index, seg, side."""
    best_d = np.full(len(lon), np.inf)
    best_f = np.full(len(lon), -1, dtype=int)
    best_s = np.zeros(len(lon), dtype=int)
    best_side = np.zeros(len(lon), dtype=int)
    for k, tr in enumerate(traces):
        d, s, side = tr.nearest(lon, lat)
        m = d < best_d
        best_d[m] = d[m]
        best_f[m] = k
        best_s[m] = s[m]
        best_side[m] = side[m]
    return best_d, best_f, best_s, best_side


# --------------------------------------------------------------------------
# text tables: keep original lines so cleaned files preserve the formatting
# --------------------------------------------------------------------------
def read_table(path, ncols=None):
    """Return numeric rows, original lines, data-line indices, raw tokens.

    ncols=None  parse every token as a number (SAR / .los files)
    ncols=N     parse only the first N tokens; trailing text such as the GPS
                site name is kept in the returned raw token list
    """
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()
    rows, idx, raw = [], [], []
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if ncols is None:
            vals = [float(v) for v in parts]
        else:
            vals = [float(v) for v in parts[:ncols]]
        rows.append(vals)
        raw.append(parts)
        idx.append(i)
    arr = np.asarray(rows, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr, lines, np.asarray(idx, dtype=int), raw


def write_filtered(src, dst, lines, data_idx, keep):
    """Copy src to dst, dropping the data rows where keep is False."""
    keep = np.asarray(keep, dtype=bool)
    if len(keep) != len(data_idx):
        raise ValueError("keep mask does not match the number of data rows")
    drop = set(data_idx[~keep].tolist())
    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    with open(dst, "w") as f:
        for i, line in enumerate(lines):
            if i not in drop:
                f.write(line)


def same_grid(a, b, tol=2e-4):
    if a.shape[0] != b.shape[0]:
        return False
    lat0 = np.radians(np.mean(a[:, 1]))
    d = np.hypot((a[:, 0] - b[:, 0]) * np.cos(lat0), a[:, 1] - b[:, 1])
    return float(d.max()) < tol


# --------------------------------------------------------------------------
# criteria
# --------------------------------------------------------------------------
def sign_flip(obs, mod, min_obs, min_model):
    return ((np.sign(obs) != np.sign(mod))
            & (np.abs(obs) >= min_obs)
            & (np.abs(mod) >= min_model))


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Remove near-fault observations whose sign is opposite "
                    "to the model prediction.")
    ap.add_argument("--fault", action="append", default=None,
                    help="PyGINVC .faultgeom file (repeat for several faults)")
    ap.add_argument("--solution", default=None,
                    help="solutions.h5 of a triangular-fault run; the surface "
                         "trace is taken from flt/vertex_llh + flt/element "
                         "(the .tri slip files store only node indices)")
    ap.add_argument("--datadir", default=".",
                    help="directory holding the observation files (default: .)")
    ap.add_argument("--modeldir", default="sm000",
                    help="model output directory relative to datadir "
                         "(default: sm000)")
    ap.add_argument("--outdir", default="clean_obs",
                    help="output directory (default: clean_obs)")
    ap.add_argument("--dist", type=float, default=1.0,
                    help="half-width of the trace band in km (default: 1.0)")
    ap.add_argument("--min-obs", type=float, default=10.0,
                    help="minimum |observation| for a sign flip, mm "
                         "(default: 10)")
    ap.add_argument("--min-model", type=float, default=30.0,
                    help="minimum |model| for a sign flip, mm "
                         "(default: 30)")
    ap.add_argument("--los", default=None,
                    help="optional source .los file to filter as well")
    ap.add_argument("--los-track", type=int, default=1,
                    help="SAR track number the .los file corresponds to "
                         "(default: 1)")
    args = ap.parse_args(argv)
    if not args.fault and not args.solution:
        ap.error("give --fault and/or --solution")

    datadir = os.path.abspath(args.datadir)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    traces = []
    if args.solution:
        import h5py
        with h5py.File(args.solution, "r") as h5:
            ftype = h5["flt/flttype"][()]
            if isinstance(ftype, bytes):
                ftype = ftype.decode()
        if ftype == "triangle":
            tri = read_triangle_traces(args.solution)
            for k, (lon, lat, info) in enumerate(tri):
                name = ("triangle" if len(tri) == 1
                        else "triangle_%d" % (k + 1))
                traces.append(FaultTrace.from_polyline(
                    lon, lat, info, name, args.solution))
        else:
            print("--solution: fault type is %r (not triangle); the "
                  "rectangular traces come from --fault" % ftype)

    for p in args.fault or []:
        nblock = len(read_fault_traces(p))
        for ib in range(nblock):
            traces.append(FaultTrace(p, block=ib))
    for tr in traces:
        print("fault %-12s %3d nodes, nsegs=%s ndeps=%s, max node gap %.2f km"
              % (tr.name, tr.info["nnodes"], tr.info["nsegs"],
                 tr.info["ndeps"], tr.info["jump_max_km"]))

    summary = []
    summary.append("near-fault sign-flip filter")
    summary.append("band half-width : %.2f km" % args.dist)
    summary.append("|obs| threshold : %.1f mm" % args.min_obs)
    summary.append("|model| threshold: %.1f mm" % args.min_model)
    summary.append("")
    summary.append("faults:")
    for tr in traces:
        summary.append("  %-12s %s  (%d nodes, nsegs=%s, ndeps=%s)"
                       % (tr.name, tr.path, tr.info["nnodes"],
                          tr.info["nsegs"], tr.info["ndeps"]))
    summary.append("")

    removed_rows = []
    qc = dict(sar=[])
    sar_clean = []          # (track, lon, lat, dist, flip, keep)
    outside_flips = []

    # ------------------------------------------------------------------
    # InSAR tracks
    # ------------------------------------------------------------------
    obs_files = sorted(
        glob.glob(os.path.join(datadir, "sar_obs_*.xyz")),
        key=lambda p: int(re.search(r"sar_obs_(\d+)", p).group(1))
        if re.search(r"sar_obs_(\d+)", p) else 0)

    for obs_path in obs_files:
        tag = re.search(r"sar_obs_(\d+)\.xyz$", obs_path).group(1)
        mod_path = os.path.join(datadir, args.modeldir,
                                "sar_mod_%s.xyz" % tag)
        res_path = os.path.join(datadir, args.modeldir,
                                "sar_res_%s.xyz" % tag)
        if not os.path.exists(mod_path):
            print("skip track %s: no %s" % (tag, mod_path))
            continue

        obs, lines, idx, _ = read_table(obs_path, ncols=3)
        mod = read_table(mod_path, ncols=3)[0]
        if not same_grid(obs, mod):
            raise ValueError("track %s: obs and model grids differ" % tag)

        obs_v = obs[:, 2]
        mod_v = mod[:, 2]
        if os.path.exists(res_path):
            res_v = read_table(res_path, ncols=3)[0][:, 2]
        else:
            res_v = obs_v - mod_v

        dist, fi, si, side = nearest_fault(traces, obs[:, 0], obs[:, 1])
        flip = sign_flip(obs_v, mod_v, args.min_obs, args.min_model)
        near = dist <= args.dist
        remove = flip & near

        for i in np.flatnonzero(flip & ~near):
            outside_flips.append((tag, i, dist[i], obs_v[i], mod_v[i]))
        for i in np.flatnonzero(remove):
            tr = traces[fi[i]]
            removed_rows.append(dict(
                kind="SAR", track=tag, index=int(i),
                lon=obs[i, 0], lat=obs[i, 1], dist_km=dist[i],
                fault=tr.name, segment=int(si[i]) + 1,
                side="HW" if side[i] > 0 else ("FW" if side[i] < 0 else "--"),
                obs=obs_v[i], mod=mod_v[i], res=res_v[i],
                obs_sign=int(np.sign(obs_v[i])),
                mod_sign=int(np.sign(mod_v[i]))))

        dst = os.path.join(outdir, "sar_obs_%s.xyz" % tag)
        write_filtered(obs_path, dst, lines, idx, ~remove)

        qc["sar"].append(dict(tag=tag, lon=obs[:, 0], lat=obs[:, 1],
                              dist=dist, flip=flip, near=near,
                              remove=remove, obs=obs_v, mod=mod_v,
                              res=res_v))
        sar_clean.append((int(tag), obs_path, dst, obs, idx, lines, remove))

        n_near = int(near.sum())
        summary.append(
            "InSAR track %s: %d points, %d within %.1f km, %d sign flips, "
            "%d removed"
            % (tag, len(obs), n_near, args.dist, int(flip.sum()),
               int(remove.sum())))
        if n_near:
            dd = dist[near]
            summary.append("  band distance  : %.3f - %.3f km (median %.3f)"
                           % (dd.min(), dd.max(), np.median(dd)))
        if remove.any():
            summary.append("  removed |obs|  : %.1f - %.1f mm"
                           % (np.abs(obs_v[remove]).min(),
                              np.abs(obs_v[remove]).max()))
            summary.append("  removed |model|: %.1f - %.1f mm"
                           % (np.abs(mod_v[remove]).min(),
                              np.abs(mod_v[remove]).max()))
            summary.append("  removed |res|  : %.1f - %.1f mm"
                           % (np.abs(res_v[remove]).min(),
                              np.abs(res_v[remove]).max()))
            hw = int((side[remove] > 0).sum())
            summary.append("  side           : %d HW, %d FW"
                           % (hw, int(remove.sum()) - hw))
        summary.append("")

    # ------------------------------------------------------------------
    # GPS
    # ------------------------------------------------------------------
    gps_specs = [
        ("gps1_obs.gmtvec", "gps1_mod.gmtvec", "EN", (2, 3)),
        ("gps1_obs_up.gmtvec", "gps1_mod_up.gmtvec", "U", (3,)),
    ]
    gps_data = []
    for obs_name, mod_name, label, cols in gps_specs:
        obs_path = os.path.join(datadir, obs_name)
        mod_path = os.path.join(datadir, args.modeldir, mod_name)
        if not (os.path.exists(obs_path) and os.path.exists(mod_path)):
            continue

        obs, lines, idx, raw = read_table(obs_path, ncols=7)
        mod = read_table(mod_path, ncols=7)[0]
        if not same_grid(obs, mod):
            raise ValueError("%s: obs and model grids differ" % obs_name)

        dist, fi, si, side = nearest_fault(traces, obs[:, 0], obs[:, 1])
        near = dist <= args.dist
        flip = np.zeros(len(obs), dtype=bool)
        comp_flip = {}
        for c in cols:
            f = sign_flip(obs[:, c], mod[:, c], args.min_obs, args.min_model)
            comp_flip[c] = f
            flip |= f

        gps_data.append(dict(name=obs_name, path=obs_path, label=label,
                             cols=cols, obs=obs, mod=mod, lines=lines,
                             idx=idx, raw=raw, dist=dist, near=near,
                             flip=flip, comp_flip=comp_flip, fi=fi, si=si,
                             side=side))

    if gps_data:
        nsta = len(gps_data[0]["obs"])
        for g in gps_data[1:]:
            if len(g["obs"]) != nsta:
                raise ValueError("GPS files do not have the same stations")
        remove = np.zeros(nsta, dtype=bool)
        for g in gps_data:
            remove |= (g["flip"] & g["near"])

        for g in gps_data:
            obs = g["obs"]
            for i in np.flatnonzero(g["flip"] & g["near"]):
                tr = traces[g["fi"][i]]
                parts = g["raw"][i]
                site = parts[7] if len(parts) > 7 else "-"
                c = [c for c in g["cols"] if g["comp_flip"][c][i]][0]
                comp = "E" if c == 2 else "N" if c == 3 else "U"
                removed_rows.append(dict(
                    kind="GPS-" + comp, track="-", index=int(i), site=site,
                    lon=obs[i, 0], lat=obs[i, 1], dist_km=g["dist"][i],
                    fault=tr.name, segment=int(g["si"][i]) + 1,
                    side=("HW" if g["side"][i] > 0
                          else "FW" if g["side"][i] < 0 else "--"),
                    obs=obs[i, c], mod=g["mod"][i, c],
                    res=obs[i, c] - g["mod"][i, c],
                    obs_sign=int(np.sign(obs[i, c])),
                    mod_sign=int(np.sign(g["mod"][i, c]))))

            dst = os.path.join(outdir, g["name"])
            write_filtered(g["path"], dst, g["lines"], g["idx"], ~remove)

            summary.append(
                "GPS %-3s: %d stations, min distance %.2f km, %d within "
                "%.1f km, %d sign flips, %d stations removed (union of "
                "components)"
                % (g["label"], nsta, float(g["dist"].min()),
                   int(g["near"].sum()), args.dist,
                   int((g["flip"] & g["near"]).sum()), int(remove.sum())))
        summary.append("")

    # ------------------------------------------------------------------
    # source .los (same points as one SAR track; match by coordinates)
    # ------------------------------------------------------------------
    los_removed = 0
    if args.los:
        track = args.los_track
        rec = [r for r in sar_clean if r[0] == track]
        if not rec:
            print("WARNING: no SAR track %d to match %s" % (track, args.los))
        else:
            _, _, _, s_obs, _, _, s_remove = rec[0]
            los, lines, idx, _ = read_table(args.los)
            if len(los) == len(s_obs):
                lat0 = np.radians(np.mean(s_obs[:, 1]))
                d = np.hypot((los[:, 0] - s_obs[:, 0]) * np.cos(lat0),
                             los[:, 1] - s_obs[:, 1])
                if d.max() < 2e-4:
                    mapping = np.arange(len(los))
                else:
                    mapping = None
            else:
                mapping = None

            if mapping is None:
                # nearest neighbour fallback
                from scipy.spatial import cKDTree
                kx = KM_PER_DEG * np.cos(np.radians(np.mean(s_obs[:, 1])))
                tree = cKDTree(np.c_[s_obs[:, 0] * kx, s_obs[:, 1] * KM_PER_DEG])
                dd, mapping = tree.query(
                    np.c_[los[:, 0] * kx, los[:, 1] * KM_PER_DEG], k=1)
                if dd.max() > 0.05:
                    raise ValueError(
                        "%s: no matching SAR points within 50 m" % args.los)

            los_remove = np.zeros(len(los), dtype=bool)
            los_remove[:] = s_remove[mapping]
            los_removed = int(los_remove.sum())

            dst = os.path.join(outdir, os.path.basename(args.los))
            write_filtered(args.los, dst, lines, idx, ~los_remove)
            summary.append(
                ".los %s: %d rows, %d removed (matched to track %d)"
                % (os.path.basename(args.los), len(los), los_removed, track))
            summary.append("")

    # ------------------------------------------------------------------
    # CSV report
    # ------------------------------------------------------------------
    csv_path = os.path.join(outdir, "removed_points.csv")
    with open(csv_path, "w") as f:
        f.write("kind,track,index,site,lon,lat,dist_km,fault,segment,side,"
                "obs_mm,mod_mm,res_mm,sign_obs,sign_mod\n")
        for r in removed_rows:
            f.write("%s,%s,%d,%s,%.6f,%.6f,%.4f,%s,%d,%s,%.2f,%.2f,%.2f,"
                    "%+d,%+d\n"
                    % (r["kind"], r["track"], r["index"], r.get("site", "-"),
                       r["lon"], r["lat"], r["dist_km"], r["fault"],
                       r["segment"], r["side"], r["obs"], r["mod"],
                       r["res"], r["obs_sign"], r["mod_sign"]))

    # ------------------------------------------------------------------
    # QC figure
    # ------------------------------------------------------------------
    if qc["sar"]:
        fig = plt.figure(figsize=(15, 5.2))
        fig.patch.set_facecolor("white")

        ax = fig.add_subplot(1, 3, 1)
        for k, tr in enumerate(traces):
            ax.plot(tr.lon, tr.lat, "k-", lw=1.4,
                    label=("trace " + tr.name) if k < 4 else None)
        first = True
        for s in qc["sar"]:
            ax.scatter(s["lon"], s["lat"], s=2, c="lightgray", rasterized=True,
                       label="all InSAR" if first else None)
            first = False
        for s in qc["sar"]:
            ax.scatter(s["lon"][s["near"]], s["lat"][s["near"]], s=10,
                       c="tab:blue", label="within band")
            ax.scatter(s["lon"][s["remove"]], s["lat"][s["remove"]], s=28,
                       facecolors="none", edgecolors="red", lw=0.9,
                       label="sign flip (removed)")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("near-fault points and removed sign flips")
        ax.legend(loc="lower right", fontsize=7, markerscale=0.8)
        ax.set_aspect("equal", adjustable="datalim")

        # zoom on the trace
        all_lon = np.concatenate([tr.lon for tr in traces])
        all_lat = np.concatenate([tr.lat for tr in traces])
        pad = 0.06
        ax2 = fig.add_subplot(1, 3, 2)
        ax2.set_xlim(all_lon.min() - pad, all_lon.max() + pad)
        ax2.set_ylim(all_lat.min() - pad, all_lat.max() + pad)
        for tr in traces:
            ax2.plot(tr.lon, tr.lat, "k-", lw=1.4)
            # tick marks showing the local dip direction
            for i in range(0, len(tr.lon) - 1, 2):
                ax2.annotate("", xy=(tr.lon[i] + 0.006 * tr.ux[i],
                                     tr.lat[i] + 0.006 * tr.uy[i]),
                             xytext=(tr.lon[i], tr.lat[i]),
                             arrowprops=dict(arrowstyle="->", lw=0.5,
                                             color="0.4"))
        for s in qc["sar"]:
            m = s["dist"] <= args.dist
            ax2.scatter(s["lon"][m], s["lat"][m], s=12,
                        c=np.where(s["flip"][m], "red", "tab:blue"),
                        zorder=3)
            for i in np.flatnonzero(s["remove"]):
                ax2.annotate(str(i), (s["lon"][i], s["lat"][i]),
                             xytext=(3, 3), textcoords="offset points",
                             fontsize=6, color="darkred")
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")
        ax2.set_title("zoom: band points (red = sign flip)")
        ax2.set_aspect("equal", adjustable="datalim")

        # obs vs model for the band points
        ax3 = fig.add_subplot(1, 3, 3)
        for s in qc["sar"]:
            m = s["near"]
            ax3.scatter(s["mod"][m], s["obs"][m], s=12,
                        c=np.where(s["flip"][m], "red", "tab:blue"),
                        label=None)
            for i in np.flatnonzero(s["remove"]):
                ax3.annotate(str(i), (s["mod"][i], s["obs"][i]),
                             xytext=(3, 3), textcoords="offset points",
                             fontsize=6, color="darkred")
        lim = 1000.0
        ax3.axhline(0, color="k", lw=0.6)
        ax3.axvline(0, color="k", lw=0.6)
        ax3.axhline(args.min_obs, color="0.6", lw=0.5, ls="--")
        ax3.axhline(-args.min_obs, color="0.6", lw=0.5, ls="--")
        ax3.axvline(args.min_model, color="0.6", lw=0.5, ls="--")
        ax3.axvline(-args.min_model, color="0.6", lw=0.5, ls="--")
        ax3.set_xlim(-lim, lim)
        ax3.set_ylim(-lim, lim)
        ax3.set_xlabel("model LOS (mm)")
        ax3.set_ylabel("observed LOS (mm)")
        ax3.set_title("band points: sign flips are in Q2/Q4")
        ax3.grid(alpha=0.25)

        fig.tight_layout()
        fig_path = os.path.join(outdir, "near_fault_qc.png")
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)
        summary.append("figure: %s" % fig_path)

    # ------------------------------------------------------------------
    txt_path = os.path.join(outdir, "filter_summary.txt")
    with open(txt_path, "w") as f:
        f.write("\n".join(summary) + "\n")

    print("\n".join(summary))
    print("\nwrote %s" % csv_path)
    print("wrote %s" % txt_path)

    if outside_flips:
        print("\nsign flips outside the %.1f km band (kept):" % args.dist)
        for tag, i, d, o, m in outside_flips:
            print("  track %s idx %d d=%.2f km obs=%.1f mod=%.1f"
                  % (tag, i, d, o, m))
    if args.los:
        print("source .los rows removed: %d" % los_removed)
    print("cleaned observations: %s" % outdir)


if __name__ == "__main__":
    sys.exit(main())
