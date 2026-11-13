#!/usr/bin/env python
# Convert FaultGeom file, output of IDSMODEL/GINVPY to input file of the Coulomb software
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Ocb 4, 2016 
# Rewritten by Zhao Bin, April 6, 2021

from pyginvc.Geometry.Patch import Fault
import argparse
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to Coulomb input file.")
    parser.add_argument('--faultfile', type=str, required=True, 
            help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    parser.add_argument('--display', type=str, required=False, default='simple',
            help='simple/complex')
    parser.add_argument('--origin', default=[], help='[lat, lon]', nargs=2, type=float)
    args    = parser.parse_args()

    fltfile = args.faultfile
    display = args.display
    origin  = args.origin
    nfault = len(np.atleast_2d(np.genfromtxt(fltfile, comments='#')))
    flt    = Fault(fltfile, nfault, 1, False, origin=origin)
    flt.load_fault()
    felem   = flt.FaultGeom2AllVertex()
    if display == 'complex':
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(felem)
    elif display == 'simple' or display == 'SIMPLE':
        """
        fmt = 12*"%10.3f"
        print(' Leng, wid, dip, strk, ss, ds, op, ts, rake, lon_mid, lat_mid, dep_mid')
        for i in range(len(flt.length)):
            print(fmt %(flt.length[i], flt.width[i], flt.dip[i], flt.strike[i],
                    flt.ss[i], flt.ds[i], flt.op[i], flt.ts[i], flt.rake[i],
                    (flt.top_left_llh[i,0]+flt.top_right_llh[i,0])/2,
                    (flt.top_left_llh[i,1]+flt.top_right_llh[i,1])/2,
                    (flt.top_left_llh[i,2]+flt.top_right_llh[i,2])/2,
                    ))"""
        selected_df = felem[["length", "width", "strike", "dip", "rake", "strike_slip", "dip_slip", "tensile_slip", "total_slip",
        "lb_lat", "lb_lon", "lb_dep", "rb_lat","rb_lon","rb_dep"]]
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(selected_df)

if __name__ == '__main__':
    main()