#!/usr/bin/env python
from pyginvc.GeomEst.GeomOptimization import GeomOptimization
import yaml, os, sys
import logging

class GeomOptim:
    '''
    GeomOptim is a class representing determination of fault geometry using different methods
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

    def run(self):
        keys = ['dict_data', 'dict_fault', 'dict_green', 'dict_weight', 'dict_param']

        for key in keys:
            if key not in self.cfg.keys():
                logging.warning('No key {} in configure file'.format(key))
                sys.exit()

        dict_data  = self.cfg['dict_data']
        dict_fault = self.cfg['dict_fault']
        dict_green = self.cfg['dict_green']
        dict_weight= self.cfg['dict_weight']
        dict_param = self.cfg['dict_param']


        opt = GeomOptimization(dict_fault, dict_data, dict_green, dict_weight, dict_param)
        if dict_param['method'] == 'emcee':
            opt.bayesian_inversion(method='emcee')
        elif dict_param['method'] == 'diff':
            opt.run_geomest()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        logging.warning('Please input configure file in YAML format')
        sys.exit()

    cfgfile = sys.argv[1]
    geomopt = GeomOptim(cfgfile)
    geomopt.run()
