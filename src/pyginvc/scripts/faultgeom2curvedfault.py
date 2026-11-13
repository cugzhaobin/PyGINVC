from pyginvc.Geometry.Patch import Fault
from pyginvc.libs import geotools as gt
import numpy as np
import argparse

def get_dim(faultfile):
    dat = np.atleast_2d(np.genfromtxt(faultfile))
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
    flt.load_fault()
    felem     = flt.FaultGeom2AllVertex()
    strike        = felem["strike"].values
    length        = felem["length"].values
    width         = felem["width"].values
    ss            = felem["strike_slip"].values
    ds            = felem["dip_slip"].values
    op            = felem["tensile_slip"].values
    # Old Fault.py named the down-dip (input) edge "top"; Patch names the
    # same corners lb/rb/rt/lt, so map them explicitly to keep the convention.
    top_left_neu  = felem[["lb_e", "lb_n", "lb_u"]].values.copy()
    top_right_neu = felem[["rb_e", "rb_n", "rb_u"]].values.copy()
    bot_right_neu = felem[["rt_e", "rt_n", "rt_u"]].values.copy()
    bot_left_neu  = felem[["lt_e", "lt_n", "lt_u"]].values.copy()
    new_dip   = np.linspace(top_dip, bot_dip, ndeps)
    for i in range(ndeps):
        for j in range(nsegs):
            k    =  i*flt.nsegs+j
            strk = strike[k]
            dip  = new_dip[i]
            [x1, y1, z1] = top_left_neu[k]
            s     = np.array([np.sin(np.deg2rad(strk)), 
                              np.cos(np.deg2rad(strk)), 0.0])
            d     = np.array([-np.cos(np.deg2rad(strk))*np.cos(np.deg2rad(dip)), 
                              np.sin(np.deg2rad(strk))*np.cos(np.deg2rad(dip)), 
                              np.sin(np.deg2rad(dip))])
#           print(x1, y1, z1) 
            x2    = x1 + length[k]*s[0]
            y2    = y1 + length[k]*s[1]
            z2    = z1

            x3    = x2 + width[k]*d[0]
            y3    = y2 + width[k]*d[1]
            z3    = z2 - width[k]*d[2]

            x4    = x1 + width[k]*d[0]
            y4    = y1 + width[k]*d[1]
            z4    = z1 - width[k]*d[2]

            top_right_neu[k] = [x2, y2, z2]
            bot_right_neu[k] = [x3, y3, z3]
            bot_left_neu[k]  = [x4, y4, z4]

            if i+1<ndeps:
                top_right_neu[k+nsegs] = [x3, y3, z3]
                top_left_neu[k+nsegs]  = [x4, y4, z4]

            llh1   = gt.utm2llh(top_left_neu[k,[1,0]], flt.origin)
            llh2   = gt.utm2llh(top_right_neu[k, [1,0]], flt.origin)

            print(f"{width[k]:.1f} {z1:.1f} {new_dip[i]:.2f} {llh1[0,0]:.4f} {llh1[0,1]:.4f} {llh2[0,0]:.4f} {llh2[0,1]:.4f} {ss[k]:.2f} {ds[k]:.2f} {op[k]:.2f}")


if __name__ == '__main__':
    main()
