#!/usr/bin/env python
# Kinematic Fault Source Inversion Program version 1.0
# This program is based on DISMODEL of Roland Burgmann
# Written By Zhao Bin, Institute of Seismology, CEA @ UC Berkeley
#
# Feb. 21 2016
# Mod by Zhao Bin, Mar  7 2016, added green function compulation
# Mod by Zhao Bin, Dec  7 2018, show figures after inversion
# Mod by Zhao Bin, Dec 15 2018, we now use dict as input parameters for inversion

import yaml, os, sys
import logging, argparse
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Triangle import Triangle
from pyginvc.Laplacian.TriLap import TriLap
from pyginvc.Inversion.TriInversion import TriInversion
from pyginvc.Export.Output import Output
from pyginvc.Export import View



class SlipInversion:
    '''
    SlipInversion is a class representing slip inversion on rectangular patches from geodetic data.
    '''
    def __init__(self, args):
        '''
        Input:
            cfgfile = YAML format configure file for slip inversion
        '''
        if os.path.isfile(args.cfgfile) is False:
            logging.warning('Input configure file {} does not exist'.format(args.cfgfile))
            sys.exit(0)

        with open(args.cfgfile, 'r') as fid:
            lines    = fid.read()
            cfg      = yaml.load(lines, Loader=yaml.FullLoader)
            self.cfg = cfg
            self.cfg['dict_data']['sarfile']   = args.sarfile if args.sarfile != "" else ""
            self.cfg['dict_data']['gpsfile']   = args.gpsfile if args.gpsfile != "" else ""
            self.cfg['dict_data']['lapmethod'] = args.lapmethod if args.lapmethod != "" else "2"
        return



    def run_inv(self):
        '''
        '''

        keys = ['dict_data', 'dict_fault', 'dict_green', 'dict_weight', 'dict_bound', 'dict_export']

        for key in keys:
            if key not in self.cfg.keys():
                logging.warning('No key {} in configure file'.format(key))
                sys.exit()

        dict_data  = self.cfg['dict_data']
        dict_fault = self.cfg['dict_fault']
        dict_green = self.cfg['dict_green']
        dict_weight= self.cfg['dict_weight']
        dict_bound = self.cfg['dict_bound']
        dict_export= self.cfg['dict_export']

        #
        # GeoData
        #
        gpsfile   = dict_data['gpsfile']
        sarfile   = dict_data['sarfile']
        levfile   = dict_data['levfile']
        gfiletype = dict_data['gfiletype']
        if 'lapmethod' in dict_data.keys():
            lapmethod = str(dict_data['lapmethod'])
        else:
            lapmethod = '2'
        data      = GeoData(gpsfile, sarfile, levfile, gfiletype)
        data.DumpData()
    
        #
        # Fault
        #
        vertexfile  = dict_fault['vertexfile']
        elementfile = dict_fault['elementfile']
        flt         = Triangle(vertexfile, elementfile, origin=[])
        flt.DumpFault()
        

        #
        # GreenFunction
        if dict_green['grnmethod'] == 'meade':
            from pyginvc.Greens.Meade import Meade
            green = Meade(flt, data, dict_green)
        elif dict_green['grnmethod'] == 'nikkhoo':
            from pyginvc.Greens.Nikkhoo import Nikkhoo
            green = Nikkhoo(flt, data, dict_green)
        elif dict_green['grnmethod'] == 'poly3d':
            from pyginvc.Greens.TriPoly3D import TriPoly3D
            green = TriPoly3D(flt, data, dict_green)
        else:
            from pyginvc.Greens.Nikkhoo import Nikkhoo
            green = Nikkhoo(flt, data, dict_green)

        
        #
        # Laplacian 
        #
        lap        = TriLap(flt.vertex_enu, flt.element, method=lapmethod)
        lap.DumpLap()

        #
        # Inversion
        #
        sol        = TriInversion(flt, data, green, lap, dict_weight, dict_bound)
        
        #
        # Output
        #
        outpath    = dict_export['outpath']
        out        = Output(flt, data, green, sol=sol, archdir=outpath)
        out.OutputSolution()
        out.archive_outfile()
        
        #
        # View observation VS modeling
        #
        if dict_export['view'] == True:
            View.plot_obs_mod(outpath+'/solutions.h5', scale=dict_export['vecscale'])
            View.plot_slip_3d(outpath+'/solutions.h5')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Invert slip distribution on triangular patches from geodetic data.")
    parser.add_argument('--cfgfile', type=str, required=True, help='')
    parser.add_argument('--sarfile', type=str, required=False, help='')
    parser.add_argument('--gpsfile', type=str, required=False, help='')
    parser.add_argument('--lapmethod', type=str, required=False, help='')
    parser.add_argument('--outpath', type=str, required=False, help='')

    args    = parser.parse_args()
    slipinv = SlipInversion(args)
    slipinv.run_inv()
