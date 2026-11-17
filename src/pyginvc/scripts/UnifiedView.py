#!/usr/bin/env python
"""Command-line viewer for unified and legacy inversion result files.

Examples::

    pyginvc unified-view --h5file run/solutions.h5
    pyginvc unified-view --h5file old/solutions.h5 --show
    pyginvc unified-view --h5file run/solutions.h5 --plot slip

The legacy flat HDF5 layout written by ``Export.Output`` is detected and
converted automatically, so old result files remain readable.
"""

import os
import argparse
import logging

from pyginvc.Unified.UnifiedView import UnifiedView

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(filename)s:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")

PLOTS = ('slip', 'obs', 'ramp', 'lcurve', 'all')


def run(h5file, plot='all', outdir=None, show=False, index=None,
        depth_exaggeration=2.0, vecscale=None):
    """Plot the selected products from one solution file."""
    viewer = UnifiedView(h5file)
    if outdir:
        viewer.outdir = outdir
        os.makedirs(outdir, exist_ok=True)

    selected = PLOTS if plot == 'all' else (plot,)
    written = []
    if 'slip' in selected:
        written.append(viewer.plot_slip_3d(
            depth_exaggeration=depth_exaggeration, show=show, index=index))
    if 'obs' in selected:
        written.append(viewer.plot_obs_mod(scale=vecscale, show=show))
    if 'ramp' in selected:
        written.append(viewer.plot_sar_ramp(show=show))
    if 'lcurve' in selected:
        written.append(viewer.plot_lcurve(show=show))
    return [path for item in written if item for path in item]


def main():
    parser = argparse.ArgumentParser(
        description='Plot unified or legacy slip inversion results.')
    parser.add_argument('--h5file', required=True,
                        help='solutions.h5 written by the inversion')
    parser.add_argument('--plot', choices=PLOTS, default='all',
                        help='which figures to produce (default: all)')
    parser.add_argument('--outdir', default=None,
                        help='figure output directory (default: h5 directory)')
    parser.add_argument('--show', action='store_true',
                        help='open interactive matplotlib windows')
    parser.add_argument('--index', type=int, default=None,
                        help='slip model index to plot (default: chosen/last)')
    parser.add_argument('--depth-exaggeration', type=float, default=2.0,
                        help='visual depth exaggeration in the 3-D view')
    parser.add_argument('--vecscale', type=float, default=None,
                        help='GPS quiver scale (default: matplotlib auto)')
    args = parser.parse_args()
    run(args.h5file, plot=args.plot, outdir=args.outdir, show=args.show,
        index=args.index, depth_exaggeration=args.depth_exaggeration,
        vecscale=args.vecscale)


if __name__ == '__main__':
    main()
