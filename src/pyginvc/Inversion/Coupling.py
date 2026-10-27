import logging
import numpy as np

from pyginvc.Blocks.Blocks import Block

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")


class Coupling:
    """
    Locking (coupling) coefficient  chi = slip deficit rate / long-term slip rate.

    Mode A -- free chi
    ------------------
    The inversion estimates an unconstrained slip deficit `s`; chi is derived
    here afterwards and is NOT confined to [0, 1].  Values outside that range are
    a result, not an error: they say the data want more (or less) slip deficit
    than the block model's long-term rate allows.

    How s_lt is obtained
    --------------------
    The long-term rate is the rigid relative rotation of the two blocks flanking
    the fault, evaluated at each patch and projected onto that patch's
    along-strike direction:

        s_lt(x) = e_ss(x) . [ B(x) . (omega_A - omega_B) ]

    `B` is the same rotation matrix and `omega` the same Euler vector the
    inversion used, so no unit conversion is needed -- s_lt comes out in the same
    mm/yr as the fitted slip deficit.  This also means Mode A inherits the
    omega<->s trade-off of the joint inversion: s_lt is not an independent
    measurement, and its amplitude is only as well constrained as omega is.

    Caveat on the dip component
    ---------------------------
    The mesh encodes dip=269 as "vertical, with a reversed sign", not as an
    angle, so the down-dip direction cannot be reconstructed from geom_grid.
    Only the along-strike component is computed.  For a near-vertical
    strike-slip fault this is the component that matters.

    Usage
    -----
        cpl = Coupling(blocks, fault)
        cpl.local_basis()                    # once
        cpl.long_term_slip(sol.euler[-1])    # needs the fitted Euler vectors
        cpl.chi_from_slip(sol.slip[-1])      # -> cpl.chi
        cpl.report()
    """

    def __init__(self, blocks, fault):
        self.blocks = blocks
        self.fault  = fault
        self.nf     = fault.nf
        self.nsegs  = fault.nsegs
        self.ndeps  = fault.ndeps

        self.e_ss   = None    # (nf, 2) unit along-strike vector, (east, north)
        self.strike = None    # (nf,)   azimuth, degrees clockwise from north

        self.v_rel  = None    # (nf, 2) rigid relative velocity at each patch
        self.s_lt   = None    # (nf,)   long-term along-strike rate
        self.mask   = None    # (nf,)   bool, True where chi is defined
        self.blk_left  = None # (nf,)   block id on the left (+y) side of the strike
        self.blk_right = None # (nf,)   block id on the right (-y) side
        self.exact  = None    # (nf,)   bool, sampled inside the polygons
        self.offset = None    # (nf,)   normal offset used, km

        self.chi    = None    # (nf,)   filled by chi_from_slip()
        self.slip   = None    # (nf, 3) the slip deficit chi was built from

    # ------------------------------------------------------------------ basis
    def local_basis(self):
        """
        Along-strike unit vector per patch.

        geom_grid columns: width, depth, dip, lat0, lon0, lat1, lon1.  The top
        edge runs (lat0, lon0) -> (lat1, lon1); both endpoints are at the surface,
        so the strike direction is horizontal and can be expressed directly in
        (east, north) without touching the dip.
        """
        g = self.fault.geom_grid
        lat0, lon0, lat1, lon1 = g[:, 3], g[:, 4], g[:, 5], g[:, 6]

        mlat = 0.5*(lat0 + lat1)
        de   = (lon1 - lon0) * 111.0 * np.cos(np.deg2rad(mlat))
        dn   = (lat1 - lat0) * 111.0
        L    = np.hypot(de, dn)
        if np.any(L <= 0):
            raise ValueError('degenerate patch geometry: zero-length top edge')

        self.e_ss   = np.column_stack([de/L, dn/L])
        self.strike = np.rad2deg(np.arctan2(de, dn)) % 360.0
        return self.e_ss

    # ------------------------------------------------------- long-term rate
    def long_term_slip(self, euler, require_exact=False, min_slt=1.0):
        """
        Long-term along-strike rate for every patch, from the fitted Euler vectors.

        euler         : (nblock, 3) or flat.  Same array the inversion returned.
        require_exact : drop patches whose flanking blocks had to be found with
                        the nearest-polygon fallback.  Off by default, because a
                        large offset is not by itself wrong -- see
                        BlockModel.flanking_blocks.  Turn it on to see the cost.
        min_slt       : |s_lt| below this (mm/yr) is treated as unresolved; the
                        ratio s/s_lt is meaningless there.

        Fills s_lt, mask, blk_a, blk_b, exact, offset.
        """
        if self.e_ss is None:
            self.local_basis()

        euler = np.asarray(euler, dtype=float).reshape(-1, 3)
        g     = self.fault.geom_grid
        lat   = 0.5*(g[:, 3] + g[:, 5])
        lon   = 0.5*(g[:, 4] + g[:, 6])

        left   = np.full(self.nf, -1, dtype=int)
        right  = np.full(self.nf, -1, dtype=int)
        exact  = np.zeros(self.nf, dtype=bool)
        offset = np.full(self.nf, np.nan)

        for k in range(self.nf):
            left[k], right[k], exact[k], offset[k] = self.blocks.flanking_blocks(
                                                lon[k], lat[k], self.strike[k])

        # Sign convention, fixed against the Okada implementation rather than
        # assumed: with dip = 269 (what this mesh uses) a unit positive ss makes
        # the +y side -- the side 90 deg counter-clockwise from the strike -- move
        # in the +strike direction.  So ss = (v_left - v_right) . e_ss, and the
        # relative velocity has to be taken in that order.
        # Verify with Greens/Okada.py: probes at dip=90, dip=270 and dip=269 give
        # opposite senses -- 269 misses the `if dip == 90 or dip == 270` branch,
        # so cos(269) < 0 flips the strike-slip sign relative to dip=90.
        v_rel = np.zeros((self.nf, 2))
        bad   = (left < 0) | (right < 0) | (left == right)
        for k in np.where(~bad)[0]:
            v_rel[k] = Block.buildG(lon[k], lat[k]) @ (euler[left[k]] - euler[right[k]])

        self.blk_left  = left
        self.blk_right = right
        self.exact  = exact
        self.offset = offset
        self.v_rel  = v_rel
        self.s_lt   = np.einsum('ij,ij->i', v_rel, self.e_ss)

        ok        = (~bad) & (np.abs(self.s_lt) >= min_slt)
        if require_exact:
            ok = ok & exact
        self.mask = ok

        n_fail   = int(bad.sum())                                  # no block pair at all
        n_weak   = int(((~bad) & (np.abs(self.s_lt) < min_slt)).sum())   # pair found, s_lt ~ 0
        fin      = offset[np.isfinite(offset)]
        logging.info('s_lt: %d/%d patches defined (%.0f%%)', ok.sum(), self.nf,
                     100.0*ok.sum()/self.nf)
        logging.info('   dropped: %d where no block pair was found, %d where |s_lt| '
                     '< %g mm/yr', n_fail, n_weak, min_slt)
        if fin.size:
            logging.info('   normal offset used to find the pair: median %.1f km, '
                         'max %.1f km  (large values mean the trace sits in a gap '
                         'between the digitised polygons)',
                         np.median(fin), fin.max())
        return self.s_lt

    # ---------------------------------------------------------------- ratio
    def chi_from_slip(self, slip):
        """
        Mode A:  chi = s_ss / s_lt, patch by patch.  Not clipped.

        slip : (nf, 3) or flat (3*nf,), i.e. sol.slip[i] / sol.x[-1][3*nblock:].
        """
        s = np.asarray(slip, dtype=float).reshape(-1, 3)
        if len(s) < self.nf:
            raise ValueError('slip has %d rows, need at least nf = %d' % (len(s), self.nf))
        s = s[:self.nf]

        self.slip = s
        self.chi  = np.full(self.nf, np.nan)
        self.chi[self.mask] = s[self.mask, 0] / self.s_lt[self.mask]
        return self.chi

    def sign_agreement(self, slip):
        """
        Fraction of patches where the fitted s_ss has the same sign as s_lt.

        A locked or partially locked fault must give ~1.0 -- the slip deficit
        accumulates in the same direction as the long-term slip.  A value near 0
        means e_ss disagrees with the sign convention of the Okada ss column, and
        every chi would come out negated: fix the basis before reading anything
        into the numbers.
        """
        s = np.asarray(slip, dtype=float).reshape(-1, 3)[:, 0]
        m = self.mask & (np.abs(s) > 1e-9)
        if not m.any():
            return np.nan, 0
        return float((np.sign(s[m]) == np.sign(self.s_lt[m])).mean()), int(m.sum())

    # --------------------------------------------------------- presentation
    def depth_profile(self, values=None):
        """Reshape a per-patch quantity to (ndeps, nsegs); the mesh is depth-major (>depth, along-strike)."""
        v = self.chi if values is None else np.asarray(values, dtype=float)
        return v.reshape(self.ndeps, self.nsegs)

    def _depth_axis(self):
        return self.fault.geom_grid[:, 1].reshape(self.ndeps, self.nsegs)[:, 0]

    def locking_depth(self, values=None, frac=0.5):
        """
        Depth above which `frac` of the total moment deficit sits.

        Weight is |chi| * patch area; the down-dip width is uniform (1 km) and
        the along-strike length cancels in the ratio, so the sum is taken over
        |chi| alone.  This is the moment-weighted measure -- it answers "how deep
        is the deficit concentrated", which is not the same question as "where
        does locking stop" when chi is not monotonic in depth.  Returns NaN when
        no patch is defined.
        """
        v = self.depth_profile(np.abs(self.chi if values is None else values))
        w = np.nan_to_num(v, nan=0.0).sum(axis=1)
        if w.sum() <= 0:
            return np.nan
        return float(np.interp(frac, np.cumsum(w)/w.sum(), self._depth_axis()))

    def half_chi_depth(self, values=None, level=0.5):
        """
        Shallowest depth below the mid-crustal peak where chi falls to `level`.

        The conventional "locking depth".  Returns NaN if chi never reaches
        `level` below its own peak -- in which case the fault is either nowhere
        locked that strongly, or locked deeper than the mesh bottom.
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
            logging.warning('chi not computed yet -- call chi_from_slip() first')
            return

        c = self.chi[self.mask]
        if c.size == 0:
            logging.warning('no patch has a defined s_lt -- chi is empty')
            return

        logging.info('chi over %d defined patches:  min %6.2f  median %6.2f  max %6.2f',
                     c.size, c.min(), np.median(c), c.max())
        logging.info('  outside [0,1]: %d (%.0f%%)   negative: %d   above 1: %d',
                     int(((c < 0) | (c > 1)).sum()),
                     100.0*((c < 0) | (c > 1)).mean(),
                     int((c < 0).sum()), int((c > 1).sum()))

        if slip is not None:
            agree, n = self.sign_agreement(slip)
            if n:
                lvl = logging.info if agree > 0.9 else logging.warning
                lvl('sign agreement s_ss vs s_lt: %.1f%% over %d patches'
                    + ('   <-- check the e_ss sign convention' if agree <= 0.9 else ''),
                    100.0*agree, n)

        prof = self.depth_profile()
        logging.info('chi depth profile (median over along-strike, %d usable columns):',
                     int(np.isfinite(prof).any(axis=0).sum()))
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
# Block backslip inversion -- Mode A locking coefficient
#
#   1 width_km    2 depth_km    3 dip_deg     4 lat0       5 lon0       6 lat1       7 lon1
#   8 s_lt_mm_yr  9 chi        10 blk_left   11 blk_right 12 offset_km
#
# chi = s_ss / s_lt.  DERIVED, not fitted, and NOT clipped to [0, 1]: a value
# outside that range says the data want more (or less) slip deficit than the
# block model's long-term rate allows.
#
# s_lt is the long-term along-strike rate at the patch, from the rigid relative
# rotation of the two blocks flanking it, using the same B and omega the
# inversion fitted.  It is therefore not an independent measurement -- it
# inherits the omega <-> s trade-off, and is only as well constrained as omega.
#
# blk_left / blk_right : block ids on the counter-clockwise / clockwise side of
#     the strike azimuth.  This pair is what fixes the sign of s_lt and of chi.
# offset_km : normal offset at which the two blocks were sampled -- the smallest
#     of 5/10/20/40 km that put both points in different polygons.  A large
#     value means the trace sits in a gap between the digitised polygons.
#
# UNDEFINED rows (no block pair, or |s_lt| below min_slt) are written as
# s_lt = chi = offset_km = 0.000000 with blk_left = blk_right = -1.  chi = 0 is
# also a legitimate value, so identify an undefined row by the block ids.
"""

    def save(self, filename='coupling.faultgeom'):
        """
        geom(7) + s_lt + chi + block_left + block_right + offset_km.

        Undefined entries are NaN in memory but written as 0.0 -- see HEADER.
        """
        g = self.fault.geom_grid
        out = np.column_stack([g[:, 0:7],
                               np.nan_to_num(self.s_lt),
                               np.nan_to_num(self.chi),
                               self.blk_left, self.blk_right,
                               np.nan_to_num(self.offset)])
        fmt = 7*['%14.6f'] + ['%14.6f', '%14.6f', '%6d', '%6d', '%8.2f']
        with open(filename, 'w') as f:
            f.write(self.HEADER)
            np.savetxt(f, out, fmt=fmt)
        return filename
