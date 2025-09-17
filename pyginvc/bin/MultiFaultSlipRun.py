#!/usr/bin/env python
# Kinematic Fault Source Inversion Program version 1.0
# This program is based on DISMODEL of Roland Burgmann
# Written By Zhao Bin, Institute of Seismology, CEA @ UC Berkeley
#
# Modified by Zhao Bin, Mar  7, 2016. added green function compulation
# Modified by Zhao Bin, Dec 22, 2016. added slip bound from files
# Modified by Zhao Bin, Apr 13, 2017. multiple fault Lap
# Modified by Zhao Bin, Apr 18, 2018. move bcs to flt_dir


# import models
#import numpy as np
import GreenFunction as gf
import Inverse as inv
import Fault, LoadData, Output
import geotools as gt
import numpy as np

####################################################################################
#                            Geodetic Inversion Python                             #
#                              Set Global Parameters                               #
####################################################################################
#                             Parameters For Fault                                 #
#----------------------------------------------------------------------------------#
# nfltplane          = numbrt of faults                                            #
# faultfile          = absolute path and file name for GPS data. It must exit.     #
# doSubFault         = True/False                                                  #
# nsegs              = numbrt of patches along fault strike angle                  #
# ndeps              = numbrt of patches along fault dip angle                     #
# bcs                = [0,0,0,0] means no slip at the four edges                   #
#                      top, bottom, end-side, start-side                           #
####################################################################################
flt_dict1            = {'faultfile': 'faultgeom_seg1', 'doSubFault': True , 'nsegs': 20, 'ndeps': 10, 'bcs':[1,0,0,0]}
flt_dict2            = {'faultfile': 'faultgeom_seg2', 'doSubFault': True , 'nsegs': 30, 'ndeps': 10, 'bcs':[1,0,0,0]}
flt_dict3            = {'faultfile': 'faultgeom_seg3', 'doSubFault': True , 'nsegs': 20, 'ndeps': 10, 'bcs':[1,0,0,0]}
flt_dict4            = {'faultfile': 'faultgeom_seg4', 'doSubFault': True , 'nsegs': 30, 'ndeps': 10, 'bcs':[1,0,0,0]}
fault_param          = [flt_dict1, flt_dict2, flt_dict3]
fault_param          = [flt_dict1, flt_dict2, flt_dict3, flt_dict4]
####################################################################################
#                          Parameters For Data                                     #
#----------------------------------------------------------------------------------#
# gfiletype          : choose GPS data format/type GMT2D or GMT3D                  #
# gpsfile            : absolute path and file name for GPS data, '' means no data  #
# levfile            : absolute path and file name for SAR data, '' means no data  #
# sarfile            : absolute path and file name for Lev data, '' means no data  #
####################################################################################
gfiletype            = 'GMT2D'    #GMT2D, GMT3D
gpsfile              = 'H_new.gmtvec'
levfile              = ''
sarfile              = ''
####################################################################################
#                          Parameters For Green Function                           #
#----------------------------------------------------------------------------------#
# greentype          : [0,1,0] means compute GF for strike dip open slip           #
# nu                 : Possion ratio                                               #
# greenmethod        : rectokada  --- classic homogenerous model                   #
#                    : triokada   --- triangular model                             #
#                    : layerwang  --- R. J. Wang's layered model                   #
####################################################################################
greentype            = [1, 1, 0]  #strike, dip, open
nu                   = 0.25
greenmethod          = 'rectokada' # not finished
greenfile            = 'SAVE'
####################################################################################
#                          Parameters For Smoothing                                #
#----------------------------------------------------------------------------------#
# smoothfactor       : first S.M.F, end S.M.F and step                             #
####################################################################################
smoothfactor         = [0.005, 0.0051, 0.1]
smoothfactor         = [0.010, 0.0101, 0.1]
####################################################################################

####################################################################################
#                          Parameters For Data Weighting                           #
#----------------------------------------------------------------------------------#
# wsar               : weight for InSAR data                                       #
# wgps               : weight for GPS   data                                       #
# wlev               : weight for Level data                                       #
####################################################################################
wgps                 = [1.0]
wsar                 = [1.0]
wlev                 = [1.0]

####################################################################################
#                          Parameters For bounding                                 #
#----------------------------------------------------------------------------------#
# bound_switch       : True/False                                                  #
# ss_range           : constrain range for strike slip                             #
# ds_range           : constrain range for dip    slip                             #
# op_range           : constrain range for opening slip                            #
# s_ss_range         : surface relative variance range for strike slip             #
# s_ds_range         : surface relative variance range for dip slip                #
# s_op_range         : surface relative variance range for open slip               #
# b_ss_range         : bottom  relative variance range for strike slip             #
# b_ds_range         : bottom  relative variance range for dip slip                #
# b_op_range         : bottom  relative variance range for open slip               #
# surf_slip          : constrained surface slip for nsegs patches                  #
# bot_slip           : constrained bottom  slip for nsegs patches                  #
# sar_plane_range    : constrained SAR parameters                                  #
# slip_lb            : lower bound of slip distribution, same format as FaultGeom  #
# slip_ub            : upper bound of slip distribution, same format as FaultGeom  #
# bound_dict         : python dictionary for bounding                              #
####################################################################################
bound_switch         = True
ss_range             = [-15000, 15000]
ds_range             = [-10000, 10000]
op_range             = [0, 0.01]
s_ss_range           = 0
s_ds_range           = 0
s_op_range           = 0 
b_ss_range           = 0
b_ds_range           = 0
b_op_range           = 0
surf_slip            = []
bot_slip             = []
sar_plane_range      = []
slip_lb              = 'faultgeom.lb'
slip_ub              = 'faultgeom.ub'
slip_lb              = ''
slip_ub              = ''
bound_dict           = {'ss_range': ss_range, 'ds_range': ds_range, 'op_range': op_range, 'slip_lb': slip_lb, 'slip_ub': slip_ub}
####################################################################################

############################################################################
#                    Read in Multiple Faults                               #
############################################################################
# get origin of multiple faults
for i in range(len(fault_param)):
    faultfile = fault_param[i]['faultfile']
    geom      = Fault.LoadFaultGeom(faultfile)
    if len(geom) == 0:
        exit
    if i == 0:
        origin    = np.array([np.mean(geom[:,[3,5]]), np.mean(geom[:,[4,6]])])
    else:
        org       = np.array([np.mean(geom[:,[3,5]]), np.mean(geom[:,[4,6]])])
        origin    = np.array([(origin[0]+org[0])*0.5, (origin[1]+org[1])*0.5])

for i in range(len(fault_param)):
    faultfile = fault_param[i]['faultfile']
    nsegs     = fault_param[i]['nsegs']
    ndeps     = fault_param[i]['ndeps']
    doSubFault= fault_param[i]['doSubFault']
    if i == 0:
        [dis_geom_grid, geom_grid] = Fault.LoadFault_origin(faultfile, nsegs, ndeps, doSubFault, origin)
    else:
        [a, b]                     = Fault.LoadFault_origin(faultfile, nsegs, ndeps, doSubFault, origin)
        dis_geom_grid              = np.vstack((dis_geom_grid, a))
        geom_grid                  = np.vstack((geom_grid, b))


############################################################################
#                    Load data                                             #
############################################################################
[llh_gps, llh_lev, llh_sar, d_gps, d_lev, d_sar, unit, W] = \
        LoadData.LoadAllData(gpsfile, sarfile, levfile, gfiletype, origin)

############################################################################
#                    Compute Green Function                                #
############################################################################
if gfiletype == 'GMT2D':
    ndim = 2
if gfiletype == 'GMT3D':
    ndim = 3
[G, G_sar] = gf.GenGreen(llh_gps, llh_lev, llh_sar, unit, origin, dis_geom_grid, greentype, nu, ndim)

nplanes = 1
# make Lap
for i in range(len(fault_param)):
    nsegs     = fault_param[i]['nsegs']
    ndeps     = fault_param[i]['ndeps']
    bcs       = fault_param[i]['bcs']
    if i == 0:
        G_lap     = gf.GenLap(nplanes, nsegs, ndeps, bcs)
    else:
        a         = gf.GenLap(nplanes, nsegs, ndeps, bcs)
        G_lap     = gt.mdiag(G_lap, a)


############################################################################
#                     Run the Inserion                                     #
############################################################################
[ruff, misfit, slip, sig_slip, r, dhat] = \
        inv.multi_inversion(d_gps, d_lev, d_sar, W, wsar, G, G_sar, G_lap, geom_grid, dis_geom_grid, fault_param, smoothfactor, bound_switch, bound_dict)


############################################################################
#                     Output Data and Figure                               #
############################################################################
Output.WriteGMTSlip(geom_grid, dis_geom_grid, slip[0,:], 'ss_slip_llh.gmtlin', 'ds_slip_llh.gmtlin', 'tt_slip_llh.gmtlin')
#Output.WriteGMTSlip(geom_grid, dis_geom_grid, sig_slip[0,:], 'sig_ss_slip_llh.gmtlin', 'sig_ds_slip_llh.gmtlin', 'sig_tt_slip_llh.gmtlin')
Output.WriteData(gpsfile, levfile, sarfile, gfiletype, dhat, r, wsar)
