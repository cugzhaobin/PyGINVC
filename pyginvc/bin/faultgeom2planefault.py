#!/usr/bin/env python
import argparse
import numpy as np
def FaultGeom2slip2D(dis_geom_grid, slip, nseg, ndep):
    '''
    convert FaultGeom fault slip model into 2D: along strike and along dip
    Written by Zhao Bin, Jan 6. 2016

    Input:
    dis_geom_grid : 
    slip          :
    nseg          : number of patches along strike direction
    ndep          : number of patches along dip direction
    Output:

    '''
    import numpy as np
    from numpy import cos, sin

    # get width and length of each patches
    wid  = dis_geom_grid[0,1]
    leng = dis_geom_grid[0,0]
    # get east and north coordinates
    east = dis_geom_grid[:,5]
    north= dis_geom_grid[:,6]
    # E-N 
    mp     = np.column_stack((east, north))
    # strike angle in radian
    strike = dis_geom_grid[0,4]*np.pi/180-np.pi/2
    # rotation matrix
    RS     = np.array([[cos(strike), -sin(strike)],[sin(strike), cos(strike)]])
    # rotate the E-N into along strike and along dip
    mp     = mp.dot(RS.T)

    # get along strike coordinate for the edge of each patch
    x1     = mp[:,0]-(0.5*leng)
    x2     = mp[:,0]+(0.5*leng)
    x      = np.vstack((x1, x2)).T

    # get along dip coordinate for the edge of each patch
    y1     = wid*np.arange(0, nseg)
    y2     = wid*np.arange(1, nseg+1)

    #
    fid1   = open('ss_xy.gmtlin', 'w')
    fid2   = open('ds_xy.gmtlin', 'w')
    fid3   = open('tt_xy.gmtlin', 'w')
    fid4   = open('locxy_slip.gmtlin', 'w')

    slip   = slip.squeeze()

    # for each patch
    for i in range(ndep):
        for j in range(nseg):
            # strike slip
            line = '> -Z%f\n' %(slip[3*i*nseg+3*j])
            fid1.write(line)
            line = '> -Z%f\n' %(slip[3*i*nseg+3*j+1])
            fid2.write(line)
            line = '> -Z%f\n' %(np.sqrt(slip[3*i*nseg+3*j]**2  + slip[3*i*nseg+3*j+1]**2))
            fid3.write(line)
            line = '%f %f\n' %(x[j,0], y1[i])
            fid1.write(line)
            fid2.write(line)
            fid3.write(line)
            line = '%f %f\n' %(x[j,1], y1[i])
            fid1.write(line)
            fid2.write(line)
            fid3.write(line)
            line = '%f %f\n' %(x[j,1], y2[i])
            fid1.write(line)
            fid2.write(line)
            fid3.write(line)
            line = '%f %f\n' %(x[j,0], y2[i])
            fid1.write(line)
            fid2.write(line)
            fid3.write(line)
            center_x = 0.5*(x[j,0]+x[j,1])
            center_y = 0.5*(y1[i]+y2[i])
            line = '%f %f %f %f\n' %(center_x, center_y, slip[3*i*nseg+3*j], slip[3*i*nseg+3*j+1] )
            fid4.write(line)

    fid1.close()
    fid2.close()
    fid3.close()
    fid4.close()

# main program
def main(args):
    from pyginvc.Geometry.Fault import Fault

    faultfile = args.faultfile
    nseg      = args.nseg
    ndep      = args.ndep

    # load the fault file
    flt = Fault(faultfile, nseg, ndep, False)
    slip = np.genfromtxt(faultfile, usecols=[7,8,9], comments='#')
    slip = slip.flatten()
    FaultGeom2slip2D(flt.dis_geom_grid, slip, nseg, ndep)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert FaultGeom file to 2D along-strike/along dip coordinate system.")
    parser.add_argument('--faultfile', type=str, required=True, help='Width depth dip lat0 lon1 lat2 lon2 strike-slip dip-slip tensile')
    parser.add_argument('--nseg', required=True, type=int)
    parser.add_argument('--ndep', required=True, type=int)
    args = parser.parse_args()
    main(args)
