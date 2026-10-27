#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Read the data, cut the profile, save the tables and the A/B figures.

    ProfilePrepare.py --cfgfile example/config_gps.yaml

No inversion is run; use ProfileInversion.py for that.
'''
import argparse
import sys

from pyginvc.Profile import run
from pyginvc.Profile.config import ConfigError


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='prepare a cross-fault profile: data + A/B figures')
    ap.add_argument('--cfgfile', required=True, help='YAML configuration')
    args = ap.parse_args(argv)
    try:
        return run.prepare(args.cfgfile)
    except (ConfigError, ValueError, OSError, ImportError) as e:
        print('error: %s' % e, file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
