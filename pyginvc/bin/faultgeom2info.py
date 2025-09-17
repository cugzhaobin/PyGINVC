#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to input file of the Coulomb software
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Ocb 4, 2016 
# Rewritten by Zhao Bin, April 6, 2021

from pyginvc.Geometry.Fault import Fault
import sys, argparse

def main(args):
    fltfile = args.faultfile
    display = args.display
    origin  = args.origin
    flt     = Fault(fltfile, 1, 1, False, origin=origin)
    if display == 'complex':
        fmt = 33*"%10.2f"
        print(" Leng, wid, dip, strk, ss, ds, op, ts, rake,\
                top_left_llh, top_right_llh,\
                bottom_right_llh, bottom_left_llh,\
                top_left_neu, top_right_neu,\
                bottom_right_neu, bottom_left_neu")
        for i in range(len(flt.length)):
            print(fmt %(flt.length[i], flt.width[i], flt.dip[i], flt.strike[i],
                    flt.ss[i], flt.ds[i], flt.op[i], flt.ts[i], flt.rake[i],
                    flt.top_left_llh[i,0], flt.top_left_llh[i,1], flt.top_left_llh[i,2],
                    flt.top_right_llh[i,0], flt.top_right_llh[i,1], flt.top_right_llh[i,2],
                    flt.bot_right_llh[i,0], flt.bot_right_llh[i,1], flt.bot_right_llh[i,2],
                    flt.bot_left_llh[i,0], flt.bot_left_llh[i,1], flt.bot_left_llh[i,2],
                    flt.top_left_neu[i,0], flt.top_left_neu[i,1], flt.top_left_neu[i,2],
                    flt.top_right_neu[i,0], flt.top_right_neu[i,1], flt.top_right_neu[i,2],
                    flt.bot_right_neu[i,0], flt.bot_right_neu[i,1], flt.bot_right_neu[i,2],
                    flt.bot_left_neu[i,0], flt.bot_left_neu[i,1], flt.bot_left_neu[i,2],
                    ))
    elif display == 'simple' or display == 'SIMPLE':
        fmt = 12*"%10.3f"
        print(' Leng, wid, dip, strk, ss, ds, op, ts, rake, lon_mid, lat_mid, dep_mid')
        for i in range(len(flt.length)):
            print(fmt %(flt.length[i], flt.width[i], flt.dip[i], flt.strike[i],
                    flt.ss[i], flt.ds[i], flt.op[i], flt.ts[i], flt.rake[i],
                    (flt.top_left_llh[i,0]+flt.top_right_llh[i,0])/2,
                    (flt.top_left_llh[i,1]+flt.top_right_llh[i,1])/2,
                    (flt.top_left_llh[i,2]+flt.top_right_llh[i,2])/2,
                    ))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to Coulomb input file.")
    parser.add_argument('--faultfile', type=str, required=True, 
            help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    parser.add_argument('--display', type=str, required=False, default='simple',
            help='simple/complex')
    parser.add_argument('--origin', default=[], help='[lat, lon]', nargs=2, type=float)
    args = parser.parse_args()
    main(args)
