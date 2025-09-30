#!/usr/bin/env python
# Forward simulation of volcano deformation
# Written By Zhao Bin, Institute of Seismology, CEA @ UC Berkeley

from pyginvc.Forward.VolForward import VolForward
import yaml, os, sys
import logging

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
        keys = ['dict_data', 'dict_volcano', 'dict_green', 'dict_weight', 'dict_export']

        for key in keys:
            if key not in self.cfg.keys():
                logging.warning('No key {} in configure file'.format(key))
                sys.exit()

        dict_data    = self.cfg['dict_data']
        dict_volcano = self.cfg['dict_volcano']
        dict_green   = self.cfg['dict_green']
        dict_weight  = self.cfg['dict_weight']
        dict_export  = self.cfg['dict_export']


        fwd = VolForward(dict_volcano, dict_data, dict_green, dict_weight)
        self.output(fwd)

        #
        # Output
        #
        outpath    = dict_export['outpath']

    def output(self, fwd):
        '''
        '''
        llh_gps = fwd.data.llh_gps
        len_gps = len(fwd.data.d_gps)
        len_sar = len(fwd.data.d_sar)
        if fwd.data.ndim == 2 and len(llh_gps)>0:
            mod = fwd.dhat.reshape((len(llh_gps),2))
            with open('gps_mod.gmtvec', 'w') as fid:
                fid.write("# Lon Lat E  N  Se  Sn  Cor site\n")
                for i in range(len(llh_gps)):
                    if len(fwd.data.station_gps) == 0:
                        fid.write('{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], mod[i,0], mod[i,1], 0.0, 0.0, 0.0))
                    else:
                        fid.write('{:10.3f} {:10.3f} {:10.3d} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], mod[i,0], mod[i,1], 0.0, 0.0, 0.0, fwd.data.station_gps[i]))
        if fwd.data.ndim == 3 and len(llh_gps)>0:
            mod = fwd.dhat.reshape((len(llh_gps),3))
            with open('gps_mod.gmtvec', 'w') as fid:
                fid.write("# Lon Lat E  N  Se  Sn  Cor site\n")
                for i in range(len(llh_gps)):
                    if len(fwd.data.station_gps) == 0:
                        fid.write('{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], mod[i,0], mod[i,1], 0.0, 0.0, 0.0))
                    else:
                        fid.write('{:10.3f} {:10.3f} {:10.3d} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], mod[i,0], mod[i,1], 0.0, 0.0, 0.0, fwd.data.station_gps[i]))
            with open('gps_mod_up.gmtvec', 'w') as fid:
                for i in range(len(llh_gps)):
                    if len(fwd.data.station_gps) == 0:
                        fid.write('{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], 0.0, mod[i,2], 0.0, 0.0, 0.0))
                    else:
                        fid.write('{:10.3f} {:10.3f} {:10.3d} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], 0.0, mod[i,2], 0.0, 0.0, 0.0, fwd.data.station_gps[i]))
            with open('gps_mod_3d.gmtvec', 'w') as fid:
                for i in range(len(llh_gps)):
                    if len(fwd.data.station_gps) == 0:
                        fid.write('{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\
                                   {:10.3f} {:10.3f} {} {:10.3f} {:10.3f}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], mod[i,0], mod[i,1], 0.0, 0.0, 0.0, "NONE", mod[i,2], 0.0))
                    else:
                        fid.write('{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\
                                   {:10.3f} {:10.3f} {} {:10.3f} {:10.3f}\n'.format(
                        llh_gps[i,1], llh_gps[i,0], mod[i,0], mod[i,1], 0.0, 0.0, 0.0, fwd.data.station_gps[i], mod[i,2], 0.0))
                


if __name__ == '__main__':
    if len(sys.argv) != 2:
        logging.warning('Please input configure file in YAML format')
        sys.exit()

    cfgfile = sys.argv[1]
    dspfwd  = DispForward(cfgfile)
    dspfwd.run_fwd()
