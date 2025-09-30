#!/usr/bin/env python
from pyginvc.Export import View
import argparse

def main(args):
    solfile  = args.solfile
    vecscale = args.vecscale

    View.plot_obs_mod(solfile, scale=vecscale)
    View.plot_slip_3d(solfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Qick view inversion results.")
    parser.add_argument('--solfile', type=str, required=True, help='solution.h5')
    parser.add_argument('--vecscale', type=float, default=1000, required=False, help='')
    args = parser.parse_args()
    main(args)
