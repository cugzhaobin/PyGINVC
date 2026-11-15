#!/usr/bin/env python
"""Plot InSAR LOS data with a cartopy basemap and a diverging colormap.

Supports one or more InSAR LOS files (format: ``Lon Lat LOS [E N U]``), draws
them either as points or on an interpolated grid, and places a horizontal
colorbar directly under each map so that map and colorbar share the same width.
"""
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #
def load_sar(path, mask_zero=False):
    """Read one InSAR file, return (lon, lat, los) with non-finite rows removed."""
    data = np.genfromtxt(path, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 3:
        raise ValueError('%s: need at least 3 columns (Lon Lat LOS)' % path)
    lon, lat, los = data[:, 0], data[:, 1], data[:, 2]
    good = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(los)
    if mask_zero:
        good &= (los != 0)
    lon, lat, los = lon[good], lat[good], los[good]
    if lon.size == 0:
        raise ValueError('%s: no valid points left after filtering' % path)
    return lon, lat, los


def data_extent(lon, lat, margin=0.05):
    """Data extent [lon0, lon1, lat0, lat1] with a margin proportional to the span."""
    dlon = max(np.ptp(lon) * margin, 1e-3)
    dlat = max(np.ptp(lat) * margin, 1e-3)
    return [lon.min() - dlon, lon.max() + dlon,
            lat.min() - dlat, lat.max() + dlat]


def union_extent(extents):
    return [min(e[0] for e in extents), max(e[1] for e in extents),
            min(e[2] for e in extents), max(e[3] for e in extents)]


# --------------------------------------------------------------------------- #
# Layout: one row per file, each map with a horizontal colorbar right below
# --------------------------------------------------------------------------- #
def _layout(extents, max_cols=3, map_h=4.5, cbar_h=0.3, cbar_gap=0.6,
            panel_gap=0.8, row_gap=1.0, ml=0.9, mr=0.35, mt=0.55, mb=0.35):
    """Compute (figsize, [(map_rect, cbar_rect), ...]) in inches/normalised units.

    Each map is sized so its width/height matches its extent's lon/lat span
    (PlateCarree keeps equal degrees), and the colorbar below spans exactly the
    same width as the map. Panels are wrapped into rows of at most ``max_cols``
    columns; each row is centred with ``panel_gap`` inches between columns and
    ``row_gap`` inches between rows.
    """
    n = len(extents)
    ncols = max(1, min(max_cols, n))
    nrows = int(np.ceil(n / ncols))
    widths = [map_h * (e[1] - e[0]) / max(e[3] - e[2], 1e-9) for e in extents]

    rows = []
    for r in range(nrows):
        idx = list(range(r * ncols, min((r + 1) * ncols, n)))
        row_w = sum(widths[i] for i in idx) + (len(idx) - 1) * panel_gap
        rows.append((idx, row_w))
    max_row_w = max(rw for _, rw in rows)

    W = ml + max_row_w + mr
    row_h = map_h + cbar_gap + cbar_h
    H = mt + nrows * row_h + (nrows - 1) * row_gap + mb

    rects = [None] * n
    top = H - mt
    for idx, row_w in rows:
        x = ml + 0.5 * (max_row_w - row_w)          # centre each row
        map_bottom = top - map_h
        cbar_bottom = map_bottom - cbar_gap - cbar_h
        for i in idx:
            w = widths[i]
            rects[i] = ([x / W, map_bottom / H, w / W, map_h / H],
                        [x / W, cbar_bottom / H, w / W, cbar_h / H])
            x += w + panel_gap
        top = cbar_bottom - row_gap
    return (W, H), rects


def _nice_ticks(lo, hi, nbins=6):
    """Return a shared set of tidy tick values spanning [lo, hi]."""
    from matplotlib.ticker import MaxNLocator
    if hi <= lo:
        hi = lo + 1.0
    ticks = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    eps = 1e-9 * max(abs(lo), abs(hi), 1.0)
    ticks = ticks[(ticks >= lo - eps) & (ticks <= hi + eps)]
    return ticks


def _add_basemap(ax, extent, grid_labels=True):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=':')
    gl = ax.gridlines(draw_labels=grid_labels, linewidth=0.2, color='gray',
                      alpha=0.4, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(args):
    datasets = []
    for path in args.sarfile:
        lon, lat, los = load_sar(path, mask_zero=args.mask_zero)
        datasets.append({
            'path': path,
            'lon': lon, 'lat': lat, 'los': los,
            'extent': data_extent(lon, lat, args.margin),
        })

    # Which extent to display
    if args.extent is not None:
        extents = [list(args.extent)] * len(datasets)
    elif args.individual_extent:
        extents = [d['extent'] for d in datasets]
    else:
        shared = union_extent([d['extent'] for d in datasets])
        extents = [shared] * len(datasets)

    # Color limits. By default every panel uses its own data range so that each
    # subplot's colorbar min/max comes from that file's data. --shared-limits
    # forces one common range; --vmin/--vmax override everything.
    def _sym_limits(vals):
        vmax_i = args.vmax if args.vmax is not None else \
            np.percentile(np.abs(vals), args.robust)
        vmin_i = args.vmin if args.vmin is not None else -vmax_i
        if vmax_i <= 0:
            vmax_i = 1.0
        if vmin_i >= vmax_i:
            vmin_i = -vmax_i
        return vmin_i, vmax_i

    def _raw_limits(vals):
        vmin_i = args.vmin if args.vmin is not None else float(vals.min())
        vmax_i = args.vmax if args.vmax is not None else float(vals.max())
        if vmax_i <= vmin_i:
            vmax_i = vmin_i + 1.0
        return vmin_i, vmax_i

    if args.shared_limits:
        allv = np.concatenate([d['los'] for d in datasets])
        fun = _raw_limits if args.data_minmax else _sym_limits
        limits = [fun(allv)] * len(datasets)
    else:
        limits = []
        for d in datasets:
            fun = _raw_limits if args.data_minmax else _sym_limits
            limits.append(fun(d['los']))

    # Per-panel tick sets (shared ticks only when the limits are shared).
    ticks_list = [_nice_ticks(lo, hi, args.cbar_ticks) for lo, hi in limits]

    max_cols = 1 if args.layout == 'vertical' else args.max_cols
    (W, H), rects = _layout(extents, max_cols=max_cols, map_h=args.height,
                            panel_gap=args.panel_gap)
    if W > args.max_width:
        s = args.max_width / W
        W, H = W * s, H * s          # 等比缩小，归一化 rect 不变
    fig = plt.figure(figsize=(W, H))

    use_map = HAS_CARTOPY and not args.no_map
    if not HAS_CARTOPY and not args.no_map:
        print('cartopy is not installed, falling back to plain axes.')

    for i, (d, ext, (map_rect, cbar_rect)) in enumerate(zip(datasets, extents, rects)):
        vmin, vmax = limits[i]
        if use_map:
            ax = fig.add_axes(map_rect, projection=ccrs.PlateCarree())
            _add_basemap(ax, ext, grid_labels=args.grid_labels)
            transform = ccrs.PlateCarree()
        else:
            ax = fig.add_axes(map_rect)
            ax.set_xlim(ext[0], ext[1])
            ax.set_ylim(ext[2], ext[3])
            ax.set_aspect('equal')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            transform = ax.transData

        lon, lat, los = d['lon'], d['lat'], d['los']
        if args.ptype == 'Point':
            mappable = ax.scatter(lon, lat, c=los, cmap=args.cmap,
                                  vmin=vmin, vmax=vmax, s=args.size,
                                  linewidths=0, transform=transform)
        else:  # Grid
            n = args.gridsize
            gx = np.linspace(lon.min(), lon.max(), n)
            gy = np.linspace(lat.min(), lat.max(), n)
            X, Y = np.meshgrid(gx, gy)
            Z = griddata((lon, lat), los, (X, Y), method='cubic')
            levels = np.linspace(vmin, vmax, args.nlevels)
            mappable = ax.contourf(X, Y, Z, levels=levels, cmap=args.cmap,
                                   extend='both', transform=transform)

        ax.set_title(os.path.basename(d['path']))

        cax = fig.add_axes(cbar_rect)
        cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
        cb.set_ticks(ticks_list[i])
        cb.set_label(args.clabel)

    if args.outfile:
        fig.savefig(args.outfile, dpi=args.dpi, bbox_inches='tight')
        print('Figure saved to %s' % args.outfile)
    if args.show:
        plt.show()
    return fig


def cli():
    parser = argparse.ArgumentParser(
        description='Plot one or more InSAR LOS files on a cartopy map '
                    'with a diverging colormap.')
    parser.add_argument('--sarfile', type=str, nargs='+', required=True,
                        help='one or more InSAR files: Lon Lat LOS [E N U]')
    parser.add_argument('--ptype', type=str, default='Point',
                        choices=['Point', 'Grid'], help='Point or Grid (default Point)')
    parser.add_argument('--layout', type=str, default='horizontal',
                        choices=['horizontal', 'vertical'],
                        help='panel arrangement (default horizontal, wrapped)')
    parser.add_argument('--max-cols', type=int, default=3,
                        help='max panels per row before wrapping (default 3)')
    parser.add_argument('--panel-gap', type=float, default=0.8,
                        help='horizontal white space between panels, inches (default 0.8)')
    parser.add_argument('--cmap', type=str, default='RdBu_r',
                        help='diverging colormap (default RdBu_r)')
    parser.add_argument('--vmin', type=float, default=None,
                        help='force color lower limit for all panels')
    parser.add_argument('--vmax', type=float, default=None,
                        help='force color upper limit for all panels')
    parser.add_argument('--robust', type=float, default=99.0,
                        help='percentile of |LOS| used for per-panel vmax (default 99)')
    parser.add_argument('--shared-limits', action='store_true',
                        help='use one common color range for all panels '
                             '(default: each panel uses its own data range)')
    parser.add_argument('--data-minmax', action='store_true',
                        help='use the raw data min/max instead of a symmetric '
                             'range around zero')
    parser.add_argument('--cbar-ticks', type=int, default=6,
                        help='target number of colorbar ticks (default 6)')
    parser.add_argument('--size', type=float, default=4.0,
                        help='scatter marker size for Point mode (default 4)')
    parser.add_argument('--gridsize', type=int, default=200,
                        help='grid resolution for Grid mode (default 200)')
    parser.add_argument('--nlevels', type=int, default=30,
                        help='contour levels for Grid mode (default 30)')
    parser.add_argument('--height', type=float, default=4.5,
                        help='map height in inches (default 4.5)')
    parser.add_argument('--max-width', type=float, default=26.0,
                        help='shrink the whole figure uniformly if it is wider '
                             'than this many inches (default 26)')
    parser.add_argument('--margin', type=float, default=0.05,
                        help='extent margin as a fraction of the span (default 0.05)')
    parser.add_argument('--extent', type=float, nargs=4, default=None,
                        metavar=('LON0', 'LON1', 'LAT0', 'LAT1'),
                        help='force a display extent for all maps')
    parser.add_argument('--individual-extent', action='store_true',
                        help='use each file\'s own extent instead of a shared one')
    parser.add_argument('--mask-zero', action='store_true',
                        help='drop points whose LOS is exactly 0 (nodata)')
    parser.add_argument('--no-map', action='store_true',
                        help='do not draw the cartopy basemap')
    parser.add_argument('--no-grid-labels', dest='grid_labels',
                        action='store_false', help='hide graticule labels')
    parser.add_argument('--clabel', type=str, default='Line of sight changes (mm)',
                        help='colorbar label')
    parser.add_argument('--outfile', type=str, default=None,
                        help='save the figure to this path')
    parser.add_argument('--dpi', type=int, default=300, help='output dpi (default 300)')
    parser.add_argument('--no-show', dest='show', action='store_false',
                        help='do not call plt.show()')
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    cli()
