#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Prepare and invert one cross-fault profile.

    ProfileInversion.py --cfgfile example/config_gps.yaml

A config with `catalog:` is handed to ProfileBatch.py's workflow.
'''
import argparse
import sys

from pyginvc.Profile import run
from pyginvc.Profile.config import ConfigError


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='prepare and invert one cross-fault profile')
    ap.add_argument('--cfgfile', required=True, help='YAML configuration')
    args = ap.parse_args(argv)
    try:
        return run.invert(args.cfgfile)
    except (ConfigError, ValueError, OSError, ImportError) as e:
        print('error: %s' % e, file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
