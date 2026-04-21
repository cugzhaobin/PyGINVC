import numpy as np
import geopandas as gpd
from pyginvc.GeoData.GeoData import GeoData
from scipy.linalg import block_diag
from scipy import optimize
import matplotlib.pyplot as plt
from pyginvc.Blocks.Blocks import BlockModel
from pyginvc.Inversion.BlockInversion import BlockInversion
import argparse

def parser_args():
    parser = argparse.ArgumentParser(description="Estimate block motion.")
    parser.add_argument('--gpsfile', type=str, required=True, help='GMT2D')
    parser.add_argument('--blockfile', type=str, required=True)
    args   = parser.parse_args()
    return args
    
def main():
    args         = parser_args()
    data         = GeoData(args.gpsfile, "", "", "GMT2D")
    data.load_data()

    blocks       = BlockModel(args.blockfile)
    nblock       = blocks.nblock
    gps_block_id = blocks.assign_point(data.llh_gps[:,[1,0]], "gps")
    data.process_data(gps_block_id, [])

    Grot_list    = blocks.buildG_block(2, "gps")
    G_gps_block  = block_diag(*Grot_list)

    sol          = BlockInversion(data, dict_weight={'wgps':[1], 'wsar':[1]})
    sol.run_inversion(nblock, G_gps_block)

    fig,ax       = plt.subplots()
    blocks.gdf.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1)
    data.plot_gps_vector(ax, data.d_gps, scale=200)
    data.plot_gps_vector(ax, sol.dhat, color='red', scale=200)
    data.plot_gps_vector(ax, sol.r, color='blue', scale=200)
    plt.show()

if __name__ == '__main__':
    main()
