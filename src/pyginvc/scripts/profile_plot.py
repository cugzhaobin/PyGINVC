#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Redraw the fit and posterior figures from the saved tables and samples.

    ProfilePlot.py --outdir example/out/XSH_kangding

The raw velocity/InSAR files are not read.
'''
import argparse
import sys

from pyginvc.Profile import run
from pyginvc.Profile.config import ConfigError


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='redraw figures from a finished profile folder')
    ap.add_argument('--outdir', required=True, help='out/PROFILE')
    args = ap.parse_args(argv)
    try:
        return run.plot(args.outdir)
    except (ConfigError, ValueError, OSError, ImportError) as e:
        print('error: %s' % e, file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
