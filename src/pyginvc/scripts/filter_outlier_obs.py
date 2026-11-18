#!/usr/bin/env python3
"""
filter_outlier_obs.py -- remove near-fault outliers whose observed sign is
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
sar_obs_<track>.xyz          observed LOS (mm), either ``lon lat los`` (as
                             written by PyGINVC) or ``lon lat los E N U``
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
sar_obs_<track>.xyz          cleaned observation files in the PyGINVC input
                             format ``lon lat los E N U``.  PyGINVC's
                             SARData reader requires all six columns, so the
                             E N U unit vectors are taken from, in order:
                             the input file itself, --los / --cfg, the run's
                             solutions.h5 (obs/llh_sar + obs/unit), a
                             6-column .los file found next to the run, or a
                             nearby solutions.h5; the points are always
                             matched by coordinates.
gps1_obs.gmtvec
gps1_obs_up.gmtvec
<name>.los                   cleaned copy of the matched source .los (or a
                             six-column file rebuilt from solutions.h5 when
                             no .los is available), so the inversion config
                             can point at it directly
removed_points.csv           every removed point with its geometry
near_fault_qc.png            map + obs/model scatter
filter_summary.txt           summary report

PyGINVC only writes three columns (``lon lat los``) into sar_obs_*.xyz, so a
plain run on a run directory has no unit vectors.  The script therefore
searches for the original six-column .los automatically.  Pass --los (or
--cfg) to pin a specific file; if nothing matches, the run aborts instead of
writing a three-column file that SARData cannot load.

Example
-------
python filter_outlier_obs.py \
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


def match_points(src_lon, src_lat, dst_lon, dst_lat, tol_km=0.05):
    """For each source point, the index of the matching destination point.

    The identity mapping is used when both sets have the same length and the
    points coincide row by row; otherwise a KD-tree nearest-neighbour lookup
    is done.  ``src`` and ``dst`` may therefore be either the .los file or a
    SAR track.
    """
    src_lon = np.asarray(src_lon, float)
    src_lat = np.asarray(src_lat, float)
    dst_lon = np.asarray(dst_lon, float)
    dst_lat = np.asarray(dst_lat, float)

    if len(src_lon) == len(dst_lon):
        lat0 = np.radians(np.mean(dst_lat))
        d = np.hypot((src_lon - dst_lon) * np.cos(lat0), src_lat - dst_lat)
        if d.max() < 2e-4:
            return np.arange(len(src_lon))

    from scipy.spatial import cKDTree
    kx = KM_PER_DEG * np.cos(np.radians(np.mean(dst_lat)))
    tree = cKDTree(np.c_[dst_lon * kx, dst_lat * KM_PER_DEG])
    dd, idx = tree.query(np.c_[src_lon * kx, src_lat * KM_PER_DEG], k=1)
    idx = np.asarray(idx)
    if dd.max() > tol_km:
        raise ValueError("points do not match within %.0f m"
                         % (tol_km * 1000.0))
    return idx


def write_sar_clean(dst, lon, lat, los, unit, keep,
                    header='Lon  Lat  LOS(mm)  E  N  U'):
    """Write the cleaned SAR track as ``lon lat los E N U``.

    This is the six-column form PyGINVC's ``SARData`` reads, so the cleaned
    file can replace the source .los file in the inversion config directly.
    """
    keep = np.asarray(keep, dtype=bool)
    data = np.column_stack((np.asarray(lon)[keep], np.asarray(lat)[keep],
                            np.asarray(los)[keep],
                            np.asarray(unit)[keep, 0],
                            np.asarray(unit)[keep, 1],
                            np.asarray(unit)[keep, 2]))
    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    np.savetxt(dst, data,
               fmt='%15.9f\t%15.9f\t%12.4f\t%15.10f\t%15.10f\t%15.10f',
               header=header)


# --------------------------------------------------------------------------
# recovering the six-column source .los
#
# PyGINVC writes sar_obs_*.xyz with only lon lat los, but its SARData reader
# hard-codes six columns (lon lat los E N U).  The unit vectors therefore
# have to be recovered from the original .los file, either the one named by
# --los/--cfg or one found near the run whose points match the observations.
# --------------------------------------------------------------------------
def read_los6(path):
    """Return (arr, lines, idx) for a six-column .los file, else None."""
    try:
        arr, lines, idx, _ = read_table(path)
    except (OSError, ValueError):
        return None
    if arr.ndim != 2 or arr.shape[1] < 6 or len(arr) == 0:
        return None
    return arr, lines, idx


def match_los_file(path, obs_lon, obs_lat, cache=None):
    """Return dict(path, arr, lines, idx, map) if path matches the obs grid.

    ``map`` maps every observation to its row in the .los file, so
    ``arr[map, 3:6]`` are the E N U unit vectors.
    """
    if cache is not None and path in cache:
        rec = cache[path]
    else:
        rec = read_los6(path)
        if cache is not None:
            cache[path] = rec
    if rec is None:
        return None
    arr, lines, idx = rec
    try:
        mapping = match_points(obs_lon, obs_lat, arr[:, 0], arr[:, 1])
    except ValueError:
        return None
    return dict(path=path, arr=arr, lines=lines, idx=idx, map=mapping)


def resolve_los(obs_lon, obs_lat, candidates, cache=None):
    """First candidate that covers the observation points.

    A file with exactly the same number of rows wins, since that is the
    original .los the run was built from rather than another downsampling.
    """
    best = None
    for path in candidates:
        rec = match_los_file(path, obs_lon, obs_lat, cache)
        if rec is None:
            continue
        exact = len(rec["arr"]) == len(obs_lon)
        if best is None or (exact and len(best["arr"]) != len(obs_lon)):
            best = rec
        if exact:
            break
    return best


def los_paths_from_cfg(cfg_path):
    """InSAR file paths listed in a PyGINVC YAML config."""
    try:
        import yaml
        with open(cfg_path, "r") as f:
            cfg = yaml.safe_load(f) or {}
    except (OSError, ImportError) as exc:
        print("WARNING: cannot read %s: %s" % (cfg_path, exc))
        return []
    bases = []
    d = os.path.dirname(os.path.abspath(cfg_path))
    while True:
        bases.append(d)
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    paths = []
    for entry in cfg.get("insar") or []:
        name = entry.get("file") if isinstance(entry, dict) else None
        if not name:
            continue
        for cand in [name] + [os.path.join(b, name) for b in bases]:
            if os.path.isfile(cand):
                paths.append(os.path.abspath(cand))
                break
    return paths


def search_root(datadir, levels=4):
    """Ancestor directory to scan for the source .los / solutions.h5 files."""
    root = os.path.abspath(datadir)
    home = os.path.abspath(os.path.expanduser("~"))
    stop = os.path.abspath(os.sep)
    for _ in range(levels):
        parent = os.path.dirname(root)
        if parent == root or parent in (home, stop):
            break
        root = parent
    return root


def find_data_candidates(root, max_depth=8):
    """(.los files, solutions.h5 files) below root, pruned of hidden dirs."""
    los, h5 = [], []
    for dirpath, dirnames, filenames in os.walk(root):
        rel = os.path.relpath(dirpath, root)
        depth = 0 if rel == "." else rel.count(os.sep) + 1
        if depth >= max_depth:
            dirnames[:] = []
        dirnames[:] = [d for d in dirnames
                       if not d.startswith(".") and d != "__pycache__"]
        if "solutions.h5" in filenames:
            h5.append(os.path.join(dirpath, "solutions.h5"))
        for name in filenames:
            if name.endswith(".los"):
                los.append(os.path.join(dirpath, name))
    return sorted(los), sorted(h5)


def read_solution_units(path):
    """E N U unit vectors and coordinates stored in solutions.h5."""
    try:
        import h5py
        with h5py.File(path, "r") as h5:
            if "obs/llh_sar" not in h5 or "obs/unit" not in h5:
                return None
            llh = np.asarray(h5["obs/llh_sar"][()], dtype=float)
            unit = np.asarray(h5["obs/unit"][()], dtype=float)
            sar_files = []
            if "obs/sar_files" in h5:
                raw = h5["obs/sar_files"][()]
                sar_files = [x.decode() if isinstance(x, bytes) else str(x)
                             for x in np.atleast_1d(raw)]
    except (OSError, KeyError, ImportError):
        return None
    if (llh.ndim != 2 or llh.shape[0] != unit.shape[0]
            or len(llh) == 0 or unit.shape[1] < 3):
        return None
    # llh_sar is stored as (lat, lon), same as SARData.llh_sar
    return dict(path=path, lon=llh[:, 1], lat=llh[:, 0], unit=unit,
                sar_files=sar_files)


def match_solution_file(path, obs_lon, obs_lat):
    """Return dict(path, unit, map, sar_files) if the h5 grid matches obs."""
    rec = read_solution_units(path)
    if rec is None:
        return None
    try:
        mapping = match_points(obs_lon, obs_lat, rec["lon"], rec["lat"])
    except ValueError:
        return None
    return dict(path=path, unit=rec["unit"], map=mapping,
                sar_files=rec["sar_files"])


def resolve_solution(obs_lon, obs_lat, candidates, cache=None):
    """First solutions.h5 whose SAR points cover the observations."""
    for path in candidates:
        if cache is not None and path in cache:
            rec = cache[path]
        else:
            rec = match_solution_file(path, obs_lon, obs_lat)
            if cache is not None:
                cache[path] = rec
        if rec is not None:
            return rec
    return None


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
                    help="optional source .los file to filter as well; its "
                         "E N U unit vectors are attached to the cleaned SAR "
                         "track (lon lat los E N U)")
    ap.add_argument("--los-track", type=int, default=1,
                    help="SAR track number the .los file corresponds to "
                         "(default: 1)")
    ap.add_argument("--cfg", default=None,
                    help="PyGINVC YAML config; its insar .los files are used "
                         "to attach E N U unit vectors and are cleaned too")
    ap.add_argument("--no-search-los", action="store_true",
                    help="do not scan the directories around --datadir for "
                         "the source six-column .los")
    ap.add_argument("--allow-missing-unit", action="store_true",
                    help="write lon lat los even when no E N U unit vectors "
                         "can be recovered (such a file cannot be read by "
                         "PyGINVC's SARData)")
    args = ap.parse_args(argv)
    if not args.fault and not args.solution:
        ap.error("give --fault and/or --solution")

    datadir = os.path.abspath(args.datadir)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    cfg_los = los_paths_from_cfg(args.cfg) if args.cfg else []
    los_cache = {}
    sol_cache = {}
    auto_data = {"root": None, "los": [], "h5": []}

    def ensure_search():
        """Lazily scan the tree around --datadir for .los / solutions.h5."""
        if auto_data["root"] is None:
            root = search_root(datadir)
            auto_data["root"] = root
            los, h5 = find_data_candidates(root)
            # --no-search-los only disables the .los search; solutions.h5 is
            # still used because it is the authoritative copy of the run data
            auto_data["los"] = [] if args.no_search_los else los
            auto_data["h5"] = h5
            print("searching %d .los and %d solutions.h5 file(s) under %s"
                  % (len(auto_data["los"]), len(h5), root))
        return auto_data["root"]

    def auto_los_paths():
        ensure_search()
        return auto_data["los"]

    def auto_h5_paths():
        """Nearby solutions.h5 files, preferring the same run-directory name."""
        ensure_search()
        own = os.path.basename(datadir.rstrip(os.sep))
        paths = list(auto_data["h5"])
        paths.sort(key=lambda p: os.path.basename(os.path.dirname(p)) != own)
        return paths

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
    los_jobs = {}           # abs .los path -> (track, match record, remove)
    sol_jobs = {}           # out basename -> rebuilt .los from solutions.h5
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

        obs, lines, idx, _ = read_table(obs_path)
        if obs.shape[1] < 3:
            raise ValueError("%s: need at least lon lat los columns"
                             % obs_path)
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

        # E N U (LOS unit vector): from the obs file itself when present,
        # otherwise from a matching .los or from solutions.h5 (obs/llh_sar
        # + obs/unit store exactly the points and look vectors of this run).
        los_rec = None
        sol_rec = None
        if obs.shape[1] >= 6:
            unit = obs[:, 3:6]
        else:
            pinned = [args.los] if (args.los is not None
                                    and int(tag) == int(args.los_track)) else []
            candidates = pinned + [p for p in cfg_los if p not in pinned]
            los_rec = resolve_los(obs[:, 0], obs[:, 1], candidates,
                                  los_cache)

            # the run's own solutions.h5 first, then nearby ones
            if los_rec is None:
                own_h5 = os.path.join(datadir, "solutions.h5")
                if os.path.exists(own_h5):
                    sol_rec = resolve_solution(obs[:, 0], obs[:, 1],
                                               [own_h5], sol_cache)
            if los_rec is None and sol_rec is None:
                los_rec = resolve_los(obs[:, 0], obs[:, 1],
                                      auto_los_paths(), los_cache)
            if los_rec is None and sol_rec is None:
                sol_rec = resolve_solution(obs[:, 0], obs[:, 1],
                                           auto_h5_paths(), sol_cache)

            if los_rec is not None:
                unit = los_rec["arr"][los_rec["map"], 3:6]
            elif sol_rec is not None:
                unit = sol_rec["unit"][sol_rec["map"]]
            elif args.allow_missing_unit:
                unit = None
                print("WARNING: %s: no E N U unit vectors found; writing "
                      "lon lat los only" % obs_path)
            else:
                raise RuntimeError(
                    "%s: cannot recover the E N U unit vectors, and "
                    "PyGINVC's SARData requires six columns (lon lat los "
                    "E N U).  Pass --los FILE or --cfg CFG, or "
                    "--allow-missing-unit to write three columns."
                    % obs_path)

        dst = os.path.join(outdir, "sar_obs_%s.xyz" % tag)
        if unit is not None:
            write_sar_clean(dst, obs[:, 0], obs[:, 1], obs_v, unit, ~remove)
        else:
            write_filtered(obs_path, dst, lines, idx, ~remove)

        if los_rec is not None:
            los_jobs.setdefault(os.path.abspath(los_rec["path"]),
                                (int(tag), los_rec, remove))
        if sol_rec is not None:
            name = (os.path.basename(sol_rec["sar_files"][0])
                    if sol_rec["sar_files"] and sol_rec["sar_files"][0]
                    else "")
            name = name or "track%s_obs.los" % tag
            taken = (name in sol_jobs
                     or os.path.abspath(os.path.join(outdir, name)) in los_jobs)
            if taken:
                name = "track%s_%s" % (tag, name)
            sol_jobs[name] = (int(tag), sol_rec, remove, obs, obs_v, unit)

        qc["sar"].append(dict(tag=tag, lon=obs[:, 0], lat=obs[:, 1],
                              dist=dist, flip=flip, near=near,
                              remove=remove, obs=obs_v, mod=mod_v,
                              res=res_v))

        n_near = int(near.sum())
        summary.append(
            "InSAR track %s: %d points, %d within %.1f km, %d sign flips, "
            "%d removed"
            % (tag, len(obs), n_near, args.dist, int(flip.sum()),
               int(remove.sum())))
        summary.append("  output columns : %s"
                       % ("lon lat los E N U" if unit is not None
                          else "lon lat los (no unit vectors)"))
        if los_rec is not None:
            summary.append("  unit vector    : %s" % los_rec["path"])
        elif sol_rec is not None:
            summary.append("  unit vector    : %s (obs/unit)"
                           % sol_rec["path"])
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
    # source .los files used for the unit vectors: write cleaned copies in
    # the original six-column format with the removed rows dropped
    # ------------------------------------------------------------------
    los_removed = 0
    for los_path in sorted(los_jobs):
        track, rec, s_remove = los_jobs[los_path]
        los = rec["arr"]
        # rec["map"] maps each SAR point to its .los row
        los_remove = np.zeros(len(los), dtype=bool)
        los_remove[rec["map"]] = s_remove
        los_removed += int(los_remove.sum())

        dst = os.path.join(outdir, os.path.basename(los_path))
        write_filtered(los_path, dst, rec["lines"], rec["idx"], ~los_remove)
        summary.append(
            ".los %s: %d rows, %d removed (unit vectors for track %d)"
            % (os.path.basename(los_path), len(los), int(los_remove.sum()),
               track))

    # when the look vectors came from solutions.h5 there is no source .los
    # to copy, so rebuild the six-column file from obs + obs/unit instead
    for name in sorted(sol_jobs):
        track, _, s_remove, obs, obs_v, unit = sol_jobs[name]
        dst = os.path.join(outdir, name)
        write_sar_clean(dst, obs[:, 0], obs[:, 1], obs_v, unit, ~s_remove,
                        header='lon lat los look_E look_N look_U')
        los_removed += int(s_remove.sum())
        summary.append(
            ".los %s: %d rows, %d removed (rebuilt from solutions.h5, "
            "track %d)" % (name, len(obs), int(s_remove.sum()), track))
    if los_jobs or sol_jobs:
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
    if los_jobs:
        print("source .los rows removed: %d" % los_removed)
    print("cleaned observations: %s" % outdir)


if __name__ == "__main__":
    sys.exit(main())
