from pyginvc.Geometry.Fault import Fault
from pyginvc.libs import geotools as gt
import numpy as np
import argparse

def get_dim(faultfile):
    dat = np.genfromtxt(faultfile)
    ndeps = int(len(np.unique(dat[:,1])))
    nsegs = int(len(dat)/ndeps)
    return nsegs, ndeps

def main():
    parser = argparse.ArgumentParser(description="Invert slip distribution on rectangular patches from geodetic data.")
    parser.add_argument('--faultfile', type=str, required=True, help='')
    parser.add_argument('--top_dip', type=float, required=True, help='')
    parser.add_argument('--bottom_dip', type=float, required=True, help='')
    parser.add_argument('--nsegs', type=int, default=0, required=False, help='')
    parser.add_argument('--ndeps', type=int, default=0, required=False, help='')
    args      = parser.parse_args()
    faultfile = args.faultfile
    top_dip   = args.top_dip
    bot_dip   = args.bottom_dip
    
    if args.nsegs>0 and args.ndeps>0:
        nsegs = args.nsegs
        ndeps = args.ndeps
    else:
        nsegs, ndeps = get_dim(faultfile)

    flt       = Fault(faultfile, nsegs, ndeps, False)
    new_dip   = np.linspace(top_dip, bot_dip, ndeps)
    for i in range(ndeps):
        for j in range(nsegs):
            k    =  i*flt.nsegs+j
            strk = flt.strike[k]
            dip  = new_dip[i]
            [x1, y1, z1] = flt.top_left_neu[k]
            s     = np.array([np.sin(np.deg2rad(strk)), 
                              np.cos(np.deg2rad(strk)), 0.0])
            d     = np.array([-np.cos(np.deg2rad(strk))*np.cos(np.deg2rad(dip)), 
                              np.sin(np.deg2rad(strk))*np.cos(np.deg2rad(dip)), 
                              np.sin(np.deg2rad(dip))])
#           print(x1, y1, z1) 
            x2    = x1 + flt.length[k]*s[0]
            y2    = y1 + flt.length[k]*s[1]
            z2    = z1

            x3    = x2 + flt.width[k]*d[0]
            y3    = y2 + flt.width[k]*d[1]
            z3    = z2 - flt.width[k]*d[2]

            x4    = x1 + flt.width[k]*d[0]
            y4    = y1 + flt.width[k]*d[1]
            z4    = z1 - flt.width[k]*d[2]

            flt.top_right_neu[k] = [x2, y2, z2]
            flt.bot_right_neu[k] = [x3, y3, z3]
            flt.bot_left_neu[k]  = [x4, y4, z4]

            if i+1<ndeps:
                flt.top_right_neu[k+nsegs] = [x3, y3, z3]
                flt.top_left_neu[k+nsegs]  = [x4, y4, z4]

            llh1   = gt.utm2llh(flt.top_left_neu[k,[1,0]], flt.origin)
            llh2   = gt.utm2llh(flt.top_right_neu[k, [1,0]], flt.origin)

            print(f"{flt.width[k]:.1f} {z1:.1f} {new_dip[i]:.2f} {llh1[0,0]:.4f} {llh1[0,1]:.4f} {llh2[0,0]:.4f} {llh2[0,1]:.4f} {flt.ss[k]:.2f} {flt.ds[k]:.2f} {flt.op[k]:.2f}")


if __name__ == '__main__':
    main()
