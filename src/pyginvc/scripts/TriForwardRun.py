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
from pyginvc.Forward.TriForward import TriForward
from pyginvc.Export.Output import Output


class DispForward:
    '''
    DispForward is a class representing forward calculation of surface displacements
    '''
    def __init__(self, cfgfile):
        '''
        Input:
            cfgfile = YAML format configure file for slip inversion
        '''
        if os.path.isfile(cfgfile) is False:
            logging.warning('Input configure file {} does not exist'.format(cfgfile))
            sys.exit(0)

        with open(cfgfile, 'r') as fid:
            lines = fid.read()
            cfg   = yaml.load(lines, Loader=yaml.FullLoader)
            self.cfg = cfg
        return

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
        dict_export= self.cfg['dict_export']


        fwd = TriForward(dict_fault, dict_data, dict_green, dict_weight)

        #
        # Output
        #
        outpath    = dict_export['outpath']
        out        = Output(fwd.flt, fwd.data, fwd.green, sol=fwd, archdir=outpath)
        out.OutputSolution()
        out.archive_outfile()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Forward modelling of surface displacements from triangular dislocation model.")
    parser.add_argument('--cfgfile', type=str, required=True, help='')
    
    args    = parser.parse_args()
    slipinv = DispForward(cfgfile)
    slipinv.run_fwd()
