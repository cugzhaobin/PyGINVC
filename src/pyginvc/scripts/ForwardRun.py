#!/usr/bin/env python
# Kinematic Fault Source Inversion Program version 1.0
# This program is based on DISMODEL of Roland Burgmann
# Written By Zhao Bin, Institute of Seismology, CEA @ UC Berkeley
#
# Feb. 21 2016
# Mod by Zhao Bin, Mar  7 2016, added green function compulation
# Mod by Zhao Bin, Dec  7 2018, show figures after inversion
# Mod by Zhao Bin, Dec 15 2018, we now use dict as input parameters for inversion
# Mod by Zhao Bin, Jul 16 2020, adopt yaml input file

import yaml, os, sys
import logging, argparse
from pyginvc.Forward.Forward import Forward
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Patch import Fault
from pyginvc.Greens.Okada import Okada


class DispForward:
    '''
    DispForward is a class representing forward calculation of surface displacements
    '''
    def __init__(self, cfgfile, args=None):
        '''
        Input:
            cfgfile = YAML format configure file for forward modelling
            args    = parsed command line arguments; when given, they override
                      the corresponding entries in cfgfile
        '''
        if os.path.isfile(cfgfile) is False:
            logging.warning('Input configure file {} does not exist'.format(cfgfile))
            sys.exit(0)

        with open(cfgfile, 'r') as fid:
            lines = fid.read()
            cfg   = yaml.load(lines, Loader=yaml.FullLoader)
            self.cfg = cfg
        self._apply_overrides(args)
        return

    def _set(self, section, key, value):
        node = self.cfg.get(section)
        if not isinstance(node, dict):
            node = {}
            self.cfg[section] = node
        node[key] = value

    def _apply_overrides(self, args):
        '''Command line options win over the configure file.'''
        if args is None:
            return
        if args.faultfile:
            self._set('dict_fault', 'faultfile', args.faultfile)
        if args.gpsfile:
            self._set('dict_data', 'gpsfile', args.gpsfile)
        if args.sarfile:
            self._set('dict_data', 'sarfile', args.sarfile)

    def run_fwd(self):
        keys = ['dict_data', 'dict_fault', 'dict_green', 'dict_weight', 'dict_bound', 'dict_export']

        for key in keys:
            if key not in self.cfg.keys():
                logging.warning('No key {} in configure file'.format(key))
                sys.exit()

        dict_data  = self.cfg['dict_data']
        dict_fault = self.cfg['dict_fault']
        dict_green = self.cfg['dict_green']
        dict_weight= self.cfg['dict_weight']

        #
        # GeoData
        #
        gpsfile = dict_data['gpsfile']
        sarfile = dict_data['sarfile']
        levfile = dict_data['levfile']
        gfiletype = dict_data['gfiletype']
        data = GeoData(gpsfile, sarfile, levfile, gfiletype)
        data.load_data()
        
        #
        # Fault
        #
        faultfile = dict_fault['faultfile']
        nsegs     = dict_fault['nsegs']
        ndeps     = dict_fault['ndeps']
        doSubFault= dict_fault['doSubFault']
        origin    = dict_fault['origin']
        flt       = Fault(faultfile, nsegs, ndeps, doSubFault, origin=origin)
        flt.load_fault()

        #
        # GreenFunction
        #
        green = Okada(flt, data, dict_green)
        green.build_greens()

        #
        # Weight
        self.wsar = dict_weight['wsar']
        data.wsar = dict_weight['wsar']

        fwd = Forward(flt, data, green)
        fwd.run_forward()

def main():
    parser = argparse.ArgumentParser(description="Forward modelling of surface displacements from rectangular dislocation model.")
    parser.add_argument('--cfgfile', type=str, required=True, help='')
    parser.add_argument('--faultfile', type=str, default='',
                        help='override dict_fault.faultfile')
    parser.add_argument('--gpsfile', type=str, default='',
                        help='override dict_data.gpsfile')
    parser.add_argument('--sarfile', type=str, default='',
                        help='override dict_data.sarfile')

    args    = parser.parse_args()
    cfgfile = args.cfgfile
    slipinv = DispForward(cfgfile, args)
    slipinv.run_fwd()
    
if __name__ == '__main__':
    main()
