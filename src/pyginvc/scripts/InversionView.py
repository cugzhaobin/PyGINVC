#!/usr/bin/env python
from pyginvc.Export import View
import argparse

def main():
    parser = argparse.ArgumentParser(description="Qick view inversion results.")
    parser.add_argument('--solfile', type=str, required=True, help='solution.h5')
    parser.add_argument('--vecscale', type=float, default=None, required=False,
                        help='vector scale; default auto (longest arrow ~10%% of map)')
    args     = parser.parse_args()

    solfile  = args.solfile
    vecscale = args.vecscale

    View.plot_lcurve(solfile)
    View.plot_sar_ramp(solfile)
    View.plot_obs_mod(solfile, scale=vecscale)
    View.plot_slip_3d(solfile)


if __name__ == '__main__':
    main()
