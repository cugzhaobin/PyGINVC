#!/usr/bin/env python
# Kinematic Fault Source Inversion Program version 1.0
# This program is based on DISMODEL of Roland Burgmann
# Written By Zhao Bin, Institute of Seismology, CEA @ UC Berkeley
#
# Feb. 21 2016
# Mar 7 2016, added green function compulation


# import models
from pyginvc.Forward.Forward import Forward
from pyginvc.Forward import Checkerboard as chb
from pyginvc.Geometry.Fault import Fault
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Greens.Okada import Okada
from pyginvc.Export.Output import Output
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap
from pyginvc.Inversion.GeoInversion import GeoInversion
import os, yaml, argparse, random, time
import numpy as np

def main(args):
    
    cfgfile = args.cfg
    with open(cfgfile, 'r') as fid:
        lines      = fid.read()
        cfg        = yaml.load(lines, Loader=yaml.FullLoader)
        dict_data  = cfg['dict_data']
        dict_fault = cfg['dict_fault']
        dict_green = cfg['dict_green']
        dict_weight= cfg['dict_weight']
        dict_export= cfg['dict_export']
        dict_bound = cfg['dict_bound']
        dict_export= cfg['dict_export']

    nss       = args.ns
    nds       = args.nd
    faultfile = dict_fault['faultfile']
    nsegs     = dict_fault['nsegs']
    ndeps     = dict_fault['ndeps']
    doSubFault= dict_fault['doSubFault']
    if doSubFault == True:
        flt   = Fault(faultfile, nsegs, ndeps, doSubFault)
        out   = np.column_stack((flt.geom_grid, np.zeros((flt.nf, 3))))
        fmt   = '%10.3f%10.3f%10.3f%10.4f%10.4f%10.4f%10.4f%15.3f%15.3f%15.3f'
        np.savetxt("tmp.faultgeom", out, fmt=fmt)
        faultfile = "tmp.faultgeom"

    # make checker board
    chb.make_checkerborad(faultfile, nsegs, ndeps, nss, nds)
    # new faultgeom file name
    newfaultfile = "FaultGeom_"+str(nss)+"by"+str(nds)
    dict_fault['faultfile'] = newfaultfile
    dict_fault['doSubFault'] = False
    Forward(dict_fault, dict_data, dict_green, dict_weight)
    sige   = np.loadtxt(dict_data['gpsfile'], usecols=[4])
    sign   = np.loadtxt(dict_data['gpsfile'], usecols=[5])
    lonlat = np.loadtxt(dict_data['gpsfile'], usecols=[0,1])
    en     = np.loadtxt('gps_mod.gmtvec', usecols=[2,3])
    for i in range(len(en)):
        en[i,0] = en[i,0] + random.normalvariate(1, 0.1)
        en[i,1] = en[i,1] + random.normalvariate(1, 0.1)
    gps_syn = np.column_stack((lonlat, en, sige, sign, np.zeros(len(en))))
    np.savetxt('gps_syn.gmtvec', gps_syn, fmt="%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f")

    data = GeoData('gps_syn.gmtvec',
            dict_data['sarfile'],
            dict_data['levfile'],
            dict_data['gfiletype'])

    flt  = Fault(newfaultfile,
            nsegs,
            ndeps,
            False)

    green= Okada(flt, data, dict_green)

    lap  = RectPlaneLap(nsegs,
            ndeps,
            dict_green['bcs'],
            method='1')

    sol  = GeoInversion(flt,
            data,
            green,
            lap,
            dict_weight,
            dict_bound)

    outpath = dict_export['outpath']

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)
        time.sleep(1)
    out     = Output(flt, data, green, sol, archdir=outpath)
    out.OutputSolution()
    out.archive_outfile()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run checkerboard test.")
    parser.add_argument('--ns', type=int, required=True, help='number of patches along strike')
    parser.add_argument('--nd', type=int, required=True, help='number of patches along dip')
    parser.add_argument('--std', type=float, required=True, help='uncertainties of obserbation')
    parser.add_argument('--cfg', type=str, required=True, help='yaml file')
    args = parser.parse_args()
    main(args)
