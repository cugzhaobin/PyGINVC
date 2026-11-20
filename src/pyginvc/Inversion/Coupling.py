import logging
import numpy as np
from pyginvc.Blocks.Blocks import Block

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s [%(filename)s:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")


class Coupling:
    """
    Locking (coupling) coefficient  chi = slip deficit rate / long-term slip rate.

    The long-term rate is *not* a single scalar taken from the along-strike
    projection.  It is the relative rigid motion of the two blocks flanking the
    fault, projected onto the fault's slip vector:

        v_rel(x) = e_ss(x) * (v_hang - v_foot) . e_ss
                 + e_up(x) * (v_hang - v_foot) . e_up

    where ``e_ss`` is the strike direction and ``e_up`` is the horizontal
    up-dip direction.  Both projections are sign-calibrated against the same
    Okada Green's function the inversion uses, so positive ``s_lt`` along a
    component is the long-term slip in the direction of a positive unit slip
    of that component.

    For a thrust this is the physically meaningful projection: the Nepal MHT
    converges at ~13-17 mm/yr down-dip, while the along-strike component is
    only ~1 mm/yr.  The old ``s_lt = v_rel . e_ss`` gave a denominator that was
    an order of magnitude too small.

    The surface trace
    -----------------
    Block polygons and Euler vectors are surface quantities.  ``s_lt`` is
    therefore evaluated at the *surface trace* of each patch (level 0 of the
    mesh), not at the depth-shifted patch midpoint.  The mesh is a translation
    copy of the surface trace, so along-strike index j of level 0 is the
    surface position of every patch in column j.

    Modes
    -----
    Linked model (default)
        The inversion estimates one chi per patch and forms the slip deficit as
        s(x) = chi(x) * s_lt(x).  This is the Meade block-backslip
        parameterisation: the fault slip direction and the long-term rate are
        fixed by the block model, and the data only decide how much of it is
        currently locked.

    Free-slip model (fallback)
        The inversion estimates a free slip vector ``s``; chi is derived here
        as the least-squares projection of ``s`` onto ``s_lt``:

            chi = (s . s_lt) / (s_lt . s_lt).

    Usage
    -----
        cpl = Coupling(blocks, fault, green=green, dict_coupling=cfg)
        cpl.local_basis()                    # once
        cpl.long_term_slip(sol.euler[-1])    # needs Euler vectors
        cpl.chi_from_slip(sol.slip[-1])      # free-slip mode only
        cpl.report()
    """

    # normal offsets used to find the two blocks on either side of the trace
    DEFAULT_OFFSETS_KM = (5., 10., 20., 40., 60., 80., 100., 120.)

    def __init__(self, blocks, fault, green=None, dict_coupling=None):
        self.blocks = blocks
        self.fault  = fault
        self.green  = green
        self.nf     = fault.nf
        self.nsegs  = fault.nsegs
        self.ndeps  = fault.ndeps

        cfg = dict(dict_coupling or {})
        self.offsets_km    = tuple(float(x) for x in cfg.get('offsets_km',
                                                            self.DEFAULT_OFFSETS_KM))
        self.footwall_block = int(cfg.get('footwall_block', -1))
        bp                  = cfg.get('block_pair') or []
        self.block_pair     = tuple(int(x) for x in bp) if len(bp) == 2 else None

        # surface trace, filled by surface_coords()
        self.lat_surf = None   # (nf,)
        self.lon_surf = None   # (nf,)

        # fault-local horizontal basis, filled by local_basis()
        self.e_ss   = None     # (nf, 2) unit along-strike vector (east, north)
        self.e_up   = None     # (nf, 2) unit horizontal UP-DIP vector
        self.e_dn   = None     # (nf, 2) unit horizontal DOWN-DIP vector
        self.strike = None     # (nf,)   azimuth, clockwise from north
        self.dip    = None     # (nf,)   physical dip angle, 0..180
        self.dip_okada = None  # (nf,)   the dip column as Okada reads it (180+angle)
        self.sign_ss = 1.0     # +1 if a positive Okada ss moves +e_ss on +y side
        self.sign_ds = 1.0     # +1 if a positive Okada ds moves the hanging wall up-dip

        # results of long_term_slip()
        self.v_rel    = None   # (nf, 2) rigid relative velocity (hanging - foot)
        self.s_lt_vec = None   # (nf, 3) long-term slip in the Okada [ss, ds, op] basis
        self.s_lt     = None   # (nf,)   magnitude sqrt(ss^2 + ds^2)
        self.mask     = None   # (nf,)   bool, True where chi is defined
        self.blk_left  = None  # (nf,)   block id on the footwall / -y side
        self.blk_right = None  # (nf,)   block id on the hanging / +y side
        self.exact  = None     # (nf,)   sampled inside the polygons
        self.offset = None     # (nf,)   normal offset used, km

        # fitted quantities
        self.chi    = None     # (nf,)   filled by chi_from_slip() or set_chi()
        self.slip   = None     # (nf, 3) the slip deficit chi was built from

    # ----------------------------------------------------------- surface trace
    def surface_coords(self):
        """
        Surface (lat, lon) of every patch.

        The mesh file is depth-major and uses the shallow-edge convention, so
        level 0 is the real surface trace.  Every deeper row is a translation
        copy of that trace; along-strike column j therefore has the same
        surface position at every depth.
        """
        g = self.fault.geom_grid
        if g.ndim != 2 or g.shape[0] != self.nf:
            raise ValueError('unexpected fault geometry shape: %r' % (g.shape,))
        if self.nf != self.nsegs * self.ndeps:
            raise ValueError('nf (%d) != nsegs*ndeps (%d*%d); the coupling '
                             'surface trace assumes a depth-major rectangular mesh'
                             % (self.nf, self.nsegs, self.ndeps))

        g0 = g.reshape(self.ndeps, self.nsegs, g.shape[1])[0]
        lat_s = 0.5*(g0[:, 3] + g0[:, 5])
        lon_s = 0.5*(g0[:, 4] + g0[:, 6])
        self.lat_surf = np.tile(lat_s, self.ndeps)
        self.lon_surf = np.tile(lon_s, self.ndeps)
        return self.lat_surf, self.lon_surf

    # ------------------------------------------------------------------ basis
    def local_basis(self):
        """
        Along-strike and up-dip unit vectors per patch, plus sign calibration.

        The mesh generator writes the trace so that the fault dips to the
        RIGHT of the strike direction (``dip_azimuth`` selection in
        trace2faultgeom.py), and writes the dip column as ``180 + angle``.
        The right normal is therefore the horizontal down-dip direction, and
        the physical dip is ``dip % 180``.

        Okada's strike-slip sign flips when ``cos(dip) < 0`` (the ``dip == 90
        or 270`` branch in Okada.py is not hit for 190).  Rather than hard-code
        that rule, the signs are measured with the same ``Okada.disloc`` the
        inversion uses; see ``_calibrate_signs``.
        """
        if self.lat_surf is None:
            self.surface_coords()

        g  = self.fault.geom_grid
        n  = self.nsegs
        g0 = g.reshape(self.ndeps, n, g.shape[1])[0]

        mlat = 0.5*(g0[:, 3] + g0[:, 5])
        de   = (g0[:, 6] - g0[:, 4]) * 111.0 * np.cos(np.deg2rad(mlat))
        dn   = (g0[:, 5] - g0[:, 3]) * 111.0
        L    = np.hypot(de, dn)
        if np.any(L <= 0):
            raise ValueError('degenerate patch geometry: zero-length top edge')

        strike_s = np.rad2deg(np.arctan2(de, dn)) % 360.0
        e_ss_s   = np.column_stack([de/L, dn/L])
        # right normal of a strike measured clockwise from north
        e_dn_s   = np.column_stack([np.cos(np.deg2rad(strike_s)),
                                    -np.sin(np.deg2rad(strike_s))])

        self.strike = np.tile(strike_s, self.ndeps)
        self.e_ss   = np.tile(e_ss_s,   (self.ndeps, 1))
        self.e_dn   = np.tile(e_dn_s,   (self.ndeps, 1))
        self.e_up   = -self.e_dn
        # dip is per patch and may vary with depth; strike and the normal are
        # tiled from the surface because the mesh is a translation copy.
        self.dip_okada = np.asarray(g[:, 2], dtype=float).copy()
        self.dip       = np.mod(self.dip_okada, 180.0)

        self._calibrate_signs()
        return self.e_ss, self.e_up

    def _calibrate_signs(self):
        """
        Measure the sign of the Okada ss/ds columns.

        A probe station is placed 20 km on the hanging-wall (+e_dn) side of a
        representative surface patch.  For a unit +ss, the displacement along
        +e_ss tells us the orientation of the strike-slip column; for a unit
        +ds, the displacement along the up-dip direction tells us whether
        positive dip-slip moves the hanging wall up-dip.

        Falls back to the analytic cos(dip) rule if Okada is unavailable or
        the probe displacement is degenerate.
        """
        # Fallback only.  The mesh writes dip = 180 + angle, which is what
        # Okada sees; cos(dip) < 0 flips the strike-slip column.  For the
        # reversed convention the +y side is the hanging wall and positive
        # dip-slip moves it up-dip; for an ordinary 0-90 dip it moves it
        # down-dip, hence the opposite sign when projected on e_up.
        flip = float(np.cos(np.deg2rad(self.dip_okada[0]))) < 0
        self.sign_ss = -1.0 if flip else 1.0
        self.sign_ds = +1.0 if flip else -1.0

        if self.green is None:
            return

        try:
            from pyginvc.Greens.Okada import Okada
            from pyginvc.libs import geotools as gt

            # Use a LEVEL-0 patch: its top edge is at the surface, so the
            # surface trace point and the Okada source are colocated.  Deeper
            # rows are translated NE along the dip and a probe measured from
            # their surface trace would land on the footwall.
            if self.nsegs < 1:
                return
            k = int(self.nsegs // 2)
            if not (np.isfinite(self.lat_surf[k]) and np.isfinite(self.lon_surf[k])):
                return

            nu   = float(self.green.dict_green.get('nu', 0.25))
            dg   = np.asarray(self.green.flt.dis_geom_grid[k], dtype=float)
            if not np.isfinite(dg).all() or dg[1] <= 0:
                return

            lat0, lon0 = self.lat_surf[k], self.lon_surf[k]
            e_dn = self.e_dn[k]
            off  = 20.0 / 111.0
            latp = lat0 + off*e_dn[1]
            lonp = lon0 + off*e_dn[0]/max(np.cos(np.deg2rad(lat0)), 1e-6)
            xy   = np.asarray(gt.llh2localxy([latp, lonp],
                                             self.green.flt.origin)).reshape(1, 2)

            src = dg.copy()
            src[7:10] = 0.0
            src[7] = 1.0
            u_ss = Okada.disloc(nu, src, xy)[0, :2]
            src[7:10] = 0.0
            src[8] = 1.0
            u_ds = Okada.disloc(nu, src, xy)[0, :2]

            p_ss = float(u_ss @ self.e_ss[k])
            p_ds = float(u_ds @ self.e_up[k])
            if abs(p_ss) > 1e-12:
                self.sign_ss = 1.0 if p_ss > 0 else -1.0
            if abs(p_ds) > 1e-12:
                self.sign_ds = 1.0 if p_ds > 0 else -1.0
            logging.info('coupling sign calibration: sign_ss=%+d sign_ds=%+d '
                         '(probe patch %d, dip=%.2f deg)',
                         int(self.sign_ss), int(self.sign_ds), k, self.dip[k])
        except Exception as err:       # never let sign calibration kill a run
            logging.warning('coupling sign calibration failed (%s); using '
                            'sign_ss=%+d sign_ds=%+d',
                            err, int(self.sign_ss), int(self.sign_ds))

    # ------------------------------------------------------- long-term rate
    def long_term_slip(self, euler, require_exact=False, min_slt=1.0,
                       verbose=True):
        """
        Long-term slip vector for every patch, from the fitted Euler vectors.

        euler         : (nblock, 3) or flat.  Same array the inversion returned.
        require_exact : drop patches whose flanking blocks had to be found with
                        the nearest-polygon fallback.
        min_slt       : |s_lt| below this (mm/yr) is treated as unresolved; the
                        ratio s/s_lt is meaningless there.
        verbose       : log the summary.  The linked inversion calls this every
                        iteration and passes verbose=False except for the last.

        Fills v_rel, s_lt_vec, s_lt, mask, blk_left, blk_right, exact, offset.
        """
        if self.e_ss is None:
            self.local_basis()
        if self.lat_surf is None:
            self.surface_coords()

        euler = np.asarray(euler, dtype=float).reshape(-1, 3)
        n     = self.nsegs

        # per along-strike column results, then tile to every depth row
        ss_col   = np.zeros(n)
        ds_col   = np.zeros(n)
        vrel_col = np.zeros((n, 2))
        left     = np.full(n, -1, dtype=int)
        right    = np.full(n, -1, dtype=int)
        exact    = np.zeros(n, dtype=bool)
        offset   = np.full(n, np.nan)

        for j in range(n):
            lon, lat = self.lon_surf[j], self.lat_surf[j]
            strike   = self.strike[j]

            if self.block_pair is not None:
                # config gives (footwall, hanging wall); +e_dn side is right
                foot, hang = self.block_pair
                left[j], right[j], exact[j], offset[j] = (
                    foot, hang, True, 0.0)
            else:
                l, r, ex, off = self.blocks.flanking_blocks(
                    lon, lat, strike, offsets_km=self.offsets_km)
                # flanking_blocks returns (left / -y = footwall,
                #                          right / +y = hanging wall)
                left[j], right[j], exact[j], offset[j] = l, r, ex, off

                if (self.footwall_block >= 0 and l >= 0 and r >= 0
                        and l != self.footwall_block and r != self.footwall_block):
                    # the two samples straddle a boundary that is not the
                    # subduction interface; do not use it
                    left[j] = right[j] = -1
                elif self.footwall_block >= 0:
                    if l == self.footwall_block:
                        pass
                    elif r == self.footwall_block:
                        left[j], right[j] = r, l

            if left[j] < 0 or right[j] < 0 or left[j] == right[j]:
                continue

            B      = Block.buildG(lon, lat)
            v_rel  = B @ (euler[right[j]] - euler[left[j]])   # hanging - footwall
            vrel_col[j] = v_rel
            ss_col[j]   = self.sign_ss * float(v_rel @ self.e_ss[j])
            ds_col[j]   = self.sign_ds * float(v_rel @ self.e_up[j])

        # tile along-strike columns to all depth rows
        self.v_rel    = np.tile(vrel_col, (self.ndeps, 1))
        s_lt_vec      = np.zeros((self.nf, 3))
        s_lt_vec[:,0] = np.tile(ss_col, self.ndeps)
        s_lt_vec[:,1] = np.tile(ds_col, self.ndeps)
        self.s_lt_vec = s_lt_vec
        self.s_lt     = np.hypot(s_lt_vec[:,0], s_lt_vec[:,1])

        self.blk_left  = np.tile(left, self.ndeps)
        self.blk_right = np.tile(right, self.ndeps)
        self.exact     = np.tile(exact, self.ndeps)
        self.offset    = np.tile(offset, self.ndeps)

        bad  = (self.blk_left < 0) | (self.blk_right < 0)
        ok   = (~bad) & (self.s_lt >= min_slt)
        if require_exact:
            ok = ok & self.exact
        self.mask = ok

        n_fail = int(bad.sum())
        n_weak = int(((~bad) & (self.s_lt < min_slt)).sum())
        fin    = self.offset[np.isfinite(self.offset) & (self.offset > 0)]
        if verbose:
            logging.info('s_lt: %d/%d patches defined (%.0f%%)', ok.sum(), self.nf,
                         100.0*ok.sum()/self.nf)
            logging.info('   dropped: %d where no block pair was found, %d where '
                         '|s_lt| < %g mm/yr', n_fail, n_weak, min_slt)
            if ok.any():
                logging.info('   median s_lt_ss = %+6.2f mm/yr   median s_lt_ds = '
                             '%+6.2f mm/yr   median |s_lt| = %5.2f mm/yr',
                             float(np.median(self.s_lt_vec[ok,0])),
                             float(np.median(self.s_lt_vec[ok,1])),
                             float(np.median(self.s_lt[ok])))
            if fin.size:
                logging.info('   normal offset used to find the pair: median %.1f km, '
                             'max %.1f km', np.median(fin), fin.max())
        return self.s_lt

    # ---------------------------------------------------------------- ratio
    def chi_from_slip(self, slip):
        """
        Free-slip mode: least-squares projection of the fitted slip deficit
        onto the long-term slip vector.

            chi = (s . s_lt) / (s_lt . s_lt)

        This is the scalar that best matches the *vector* s, not one arbitrarily
        chosen component.  For a pure thrust it reduces to s_ds / s_lt_ds.
        """
        s = np.asarray(slip, dtype=float).reshape(-1, 3)
        if len(s) < self.nf:
            raise ValueError('slip has %d rows, need at least nf = %d'
                             % (len(s), self.nf))
        s = s[:self.nf]

        self.slip = s
        self.chi  = np.full(self.nf, np.nan)
        if self.s_lt_vec is None:
            raise ValueError('call long_term_slip() before chi_from_slip()')
        num = np.einsum('ij,ij->i', s, self.s_lt_vec)
        den = np.einsum('ij,ij->i', self.s_lt_vec, self.s_lt_vec)
        ok  = self.mask & (den > 0)
        self.chi[ok] = num[ok] / den[ok]
        return self.chi

    def set_chi(self, chi):
        """Linked mode: the inversion fitted chi directly."""
        chi = np.asarray(chi, dtype=float).reshape(-1)
        if chi.size < self.nf:
            raise ValueError('chi has %d entries, need at least nf = %d'
                             % (chi.size, self.nf))
        self.chi = chi[:self.nf].copy()
        return self.chi

    def sign_agreement(self, slip):
        """
        Fraction of patches whose fitted slip deficit points along the
        long-term slip vector (positive projection).

        A locked or partially locked fault must give ~1.0.  A value near 0
        means the slip deficit and s_lt have opposite senses, which usually
        means the block order or the Okada sign calibration is wrong.
        """
        s = np.asarray(slip, dtype=float).reshape(-1, 3)[:self.nf]
        proj = np.einsum('ij,ij->i', s, self.s_lt_vec)
        m = self.mask & (np.hypot(s[:,0], s[:,1]) > 1e-9)
        if not m.any():
            return np.nan, 0
        return float((proj[m] > 0).mean()), int(m.sum())

    # --------------------------------------------------------- presentation
    def depth_profile(self, values=None):
        """Reshape a per-patch quantity to (ndeps, nsegs); the mesh is depth-major."""
        v = self.chi if values is None else np.asarray(values, dtype=float)
        return v.reshape(self.ndeps, self.nsegs)

    def _depth_axis(self):
        return self.fault.geom_grid[:, 1].reshape(self.ndeps, self.nsegs)[:, 0]

    def locking_depth(self, values=None, frac=0.5):
        """
        Depth above which `frac` of the total moment deficit sits.

        Weight is |chi| * patch area; the down-dip width is uniform (1 km) and
        the along-strike length cancels in the ratio, so the sum is taken over
        |chi| alone.  This is the moment-weighted measure.
        """
        v = self.depth_profile(np.abs(self.chi if values is None else values))
        w = np.nan_to_num(v, nan=0.0).sum(axis=1)
        if w.sum() <= 0:
            return np.nan
        return float(np.interp(frac, np.cumsum(w)/w.sum(), self._depth_axis()))

    def half_chi_depth(self, values=None, level=0.5):
        """
        Shallowest depth below the mid-crustal peak where chi falls to `level`.
        """
        v = self.depth_profile(self.chi if values is None else values)
        depth = self._depth_axis()
        prof  = np.nanmedian(v, axis=1)
        if not np.isfinite(prof).any():
            return np.nan
        ipk = int(np.nanargmax(prof))
        below = np.where(np.isfinite(prof[ipk:]) & (prof[ipk:] <= level))[0]
        if below.size == 0:
            return np.nan
        return float(depth[ipk + below[0]])

    def segment_mean(self, nseg=8, values=None):
        """Along-strike means over nseg contiguous column groups -> (nseg,) + edges."""
        v     = self.depth_profile(self.chi if values is None else values)
        edges = np.linspace(0, self.nsegs, nseg+1).astype(int)
        out   = np.full(nseg, np.nan)
        for i in range(nseg):
            blk = v[:, edges[i]:edges[i+1]]
            if np.isfinite(blk).any():
                out[i] = np.nanmean(blk)
        return out, edges

    def report(self, slip=None):
        """Log a one-screen summary of the coupling distribution."""
        if self.chi is None:
            logging.warning('chi not computed yet -- call chi_from_slip() / '
                            'set_chi() first')
            return

        c = self.chi[self.mask]
        if c.size == 0:
            logging.warning('no patch has a defined s_lt -- chi is empty')
            return

        logging.info('chi over %d defined patches:  min %6.2f  median %6.2f  '
                     'max %6.2f', c.size, c.min(), np.median(c), c.max())
        logging.info('  outside [0,1]: %d (%.0f%%)   negative: %d   above 1: %d',
                     int(((c < 0) | (c > 1)).sum()),
                     100.0*((c < 0) | (c > 1)).mean(),
                     int((c < 0).sum()), int((c > 1).sum()))

        if slip is not None:
            agree, n = self.sign_agreement(slip)
            if n:
                lvl = logging.info if agree > 0.9 else logging.warning
                lvl('slip deficit points along s_lt on %.1f%% of %d patches'
                    + ('   <-- check block order / sign calibration'
                       if agree <= 0.9 else ''),
                    100.0*agree, n)

        prof = self.depth_profile()
        logging.info('chi depth profile (median over along-strike):')
        depth = self.fault.geom_grid[:, 1].reshape(self.ndeps, self.nsegs)[:, 0]
        for iz in range(self.ndeps):
            row = prof[iz][np.isfinite(prof[iz])]
            logging.info('   z = %5.2f km   n = %2d   median chi = %s',
                         depth[iz], row.size,
                         '%6.2f' % np.median(row) if row.size else '   n/a')

        ld = self.locking_depth()
        logging.info('depth above which 50%% of the moment deficit sits: %s km',
                     '%.2f' % ld if np.isfinite(ld) else 'n/a')
        hd = self.half_chi_depth()
        logging.info('locking depth (median chi falls back to 0.5): %s km',
                     '%.2f' % hd if np.isfinite(hd) else 'n/a')

    # ------------------------------------------------------------------ IO
    HEADER = """\
# Block backslip inversion -- locking coefficient
#
#   1 width_km    2 depth_km    3 dip_deg     4 lat0       5 lon0       6 lat1       7 lon1
#   8 s_lt_ss     9 s_lt_ds    10 |s_lt|      11 chi       12 blk_left   13 blk_right  14 offset_km
#
# s_lt_ss / s_lt_ds are the long-term strike-slip and dip-slip rates of the
# hanging wall relative to the footwall, projected onto the fault's slip
# vector and sign-calibrated against the Okada Green's function.  |s_lt| is
# sqrt(ss^2 + ds^2).  Units mm/yr.
#
# chi is FITTED in the linked model (s = chi * s_lt) and DERIVED in the
# free-slip model (chi = s . s_lt / |s_lt|^2).  In the linked model it is
# bounded by dict_coupling.chi_range (default [0, 1]).
#
# blk_left / blk_right : block ids on the footwall (-y) and hanging-wall (+y)
#     sides of the strike.  This pair fixes the sign of s_lt and chi.
# offset_km : normal offset at which the two blocks were sampled at the
#     SURFACE trace.  A large value means the trace is drawn away from the
#     digitised block boundary, but the rigid relative rotation is still the
#     long-term rate as long as the pair is the correct one.
#
# UNDEFINED rows (no block pair, or |s_lt| below min_slt) are written as
# s_lt = chi = offset_km = 0.000000 with blk_left = blk_right = -1.
"""

    def save(self, filename='coupling.faultgeom'):
        """
        geom(7) + s_lt_ss + s_lt_ds + |s_lt| + chi + block_left + block_right
        + offset_km.  Undefined entries are written as 0.0 -- see HEADER.
        """
        g    = self.fault.geom_grid
        sv   = self.s_lt_vec if self.s_lt_vec is not None else np.zeros((self.nf, 3))
        slt  = self.s_lt     if self.s_lt     is not None else np.zeros(self.nf)
        chi  = self.chi      if self.chi      is not None else np.full(self.nf, np.nan)
        left = self.blk_left  if self.blk_left  is not None else np.full(self.nf, -1)
        right = self.blk_right if self.blk_right is not None else np.full(self.nf, -1)
        off  = self.offset   if self.offset   is not None else np.full(self.nf, np.nan)
        out = np.column_stack([g[:, 0:7],
                               np.nan_to_num(sv[:, 0]),
                               np.nan_to_num(sv[:, 1]),
                               np.nan_to_num(slt),
                               np.nan_to_num(chi),
                               left, right,
                               np.nan_to_num(off)])
        fmt = 7*['%14.6f'] + 3*['%14.6f'] + ['%14.6f', '%6d', '%6d', '%8.2f']
        with open(filename, 'w') as f:
            f.write(self.HEADER)
            np.savetxt(f, out, fmt=fmt)
        return filename
