#!/usr/bin/env python
from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Patch import Fault
from pyginvc.Greens.Okada import Okada
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap
from pyginvc.Inversion.BlockBackslipInversion import BlockBackslipInversion
from pyginvc.Export.Output import Output
from pyginvc.Export import View
from pyginvc.Blocks.Blocks import BlockModel
import numpy as np
from scipy.linalg import block_diag
import logging
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Invert for block rotation and back slip distribution on rectangular patches from geodetic data.")
    parser.add_argument('--cfgfile', type=str, required=True, help='')
    
    args    = parser.parse_args()
    inv(args) 


def inv(args):
    #
    # GeoData
    #
    gpsfile    = args.gpsfile
    data       = GeoData(gpsfile, "", "", "GMT2D")
    data.load_data()
    lapmethod  = '1'

    #
    # Blocks
    #
    blkfile    = args.blockfile
    blocks     = BlockModel(blkfile)
    gps_block_id = blocks.assign_point(data.llh_gps[:,[1,0]], "gps")
    data.process_data(gps_block_id, [])

    # Green's
    G_gps_list  = blocks.buildG_block(data.ndim, block_type="gps")
    G_gps_block = block_diag(*G_gps_list)


    #
    # Fault
    #
    origin    = args.origin
    faultfile = args.faultfile
    nsegs     = args.nsegs
    ndeps     = args.ndeps
    flt       = Fault(seg1, 40, 15, False, origin=origin)
    flt1.load_fault()

    dict_green = args.dict_green
    green = Okada(flt1, data, dict_green)
    green.build_greens()

    #
    # Laplacian
    #
    dict_bound = args.dict_bound
    bcs        = dict_green['bcs']
    lap        = RectPlaneLap(flt1.nsegs, flt1.ndeps, bcs, method=lapmethod)
    lap.build_laps()

    #
    # Inversion
    #
    dict_weight = args.dict_wegiht

    nblock       = blocks.nblock
    sol = BlockBackslipInversion(data, flt1, green, lap, dict_weight, dict_bound)
    sol.run_inversion(nblock, G_gps_block)

    plt.subplot(1,2,1)
    ax = plt.gca()
    blocks.gdf.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1)
    data.plot_gps_vector(ax, sol.r, scale=200, color='red')

    plt.subplot(1,2,2)
    ax = plt.gca()
    blocks.gdf.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1)
    data.plot_gps_vector(ax, sol.dhat, scale=200, color='red')
    plt.show()

    if flt1 is not None:
        flt1.plot_faultgeom(value=sol.x[3*nblock:][0::3])
        outdata = np.column_stack([data.llh_gps[:,1], data.llh_gps[:,0], sol.gps_mod_slip.reshape(-1,2)])
        np.savetxt('gps_mod_slip.gmtvec', outdata)

        slip = sol.x[3*nblock:].reshape(-1,3)
        outdata = np.column_stack([flt1.geom_grid[:,0:7], slip])
        np.savetxt('slipmodel.faultgeom', outdata, fmt=10*' %10.3f')

if __name__ == '__main__':
    main()
