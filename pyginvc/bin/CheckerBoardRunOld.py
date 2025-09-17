#!/usr/bin/env python
# Kinematic Fault Source Inversion Program version 1.0
# This program is based on DISMODEL of Roland Burgmann
# Written By Zhao Bin, Institute of Seismology, CEA @ UC Berkeley
#
# Feb. 21 2016
# Mar 7 2016, added green function compulation


# import models
#import numpy as np
import Inverse as inv
import Forward as fwd
import Checkerboard as chb
import os

####################################################################################
#                            Geodetic CheckerBroad Python                          #
#                              Set Global Parameters                               #
####################################################################################
#                             Parameters For Fault                                 #
#----------------------------------------------------------------------------------#
# nflts              = numbrt of faults                                            #
# faultfile          = absolute path and file name for GPS data. It must exit.     #
#                    : format: wid dep dip lat1 lon1 lat2 lon2 ss ds op in mm      #
#                    : ss>0 right-lateral slip | ss<0 left-lateral slip            #
#                    : ds>0 inverse slip       | ds<0 normal slip                  #
# nsegs              = numbrt of patches along fault strike angle                  #
# ndeps              = numbrt of patches along fault dip angle                     #
####################################################################################
nflts                = 1
faultfile            = '/Users/zhao/earthquake/nepal/postseismic/afterslip/dismodel2D/FaultGeom'
doSubFault           = True
nsegs                = 10
ndeps                = 10
####################################################################################
#                          Parameters For Data                                     #
#----------------------------------------------------------------------------------#
# gpsfile            : absolute path and file name for GPS data, '' means no data  #
#                    : GMT2D format: lon lat ve vn se sn cor site in mm            #
#                    : GMT3D format: lon lat ve vn vu se sn su site in mm          #
# levfile            : absolute path and file name for SAR data, '' means no data  #
# sarfile            : absolute path and file name for Lev data, '' means no data  #
#                    : INSAR format: lon lat los unit_e unit_n unit_u in mm        #
####################################################################################
gfiletype            = 'GMT2D'    #GMT2D, GMT3D
gpsfile              = '/Users/zhao/earthquake/nepal/postseismic/afterslip/dismodel2D/nepal_post_90_new.gmtvec'
levfile              = ''
sarfile              = ''

####################################################################################
#                          Parameters For Green Function                           #
#----------------------------------------------------------------------------------#
# greentype          : [0,1,0] means compute GF for strike dip open slip           #
# nu                 : Possion ratio                                               #
# bcs                : [0,0,0,0] means no slip at the four edges                   #
# greenmethod        : rectokada  --- classic homogenerous model                   #
#                    : triokada   --- triangular model                             #
#                    : layerwang  --- R. J. Wang's layered model                   #
####################################################################################
greentype            = [0, 1, 0]  #strike, dip, open
nu                   = 0.25
bcs                  = [0,0,0,0]
greenmethod          = 'rectokada' # not finished

####################################################################################
#                          Parameters For Smoothing                                #
#----------------------------------------------------------------------------------#
# smoothfactor       : first S.M.F, end S.M.F and step                             #
####################################################################################
smoothfactor         = [0.1, 0.11, 0.1]
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
# bound_dict         : python dictionary for bounding                              #
####################################################################################
bound_switch         = True
ss_range             = [0, 0.01]
ds_range             = [0,  200]
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
bound_dict = {'ss_range': ss_range, 'ds_range': ds_range, 'op_range': op_range}
####################################################################################
####################################################################################
#                          Parameters Checkerbroad                                 #
#----------------------------------------------------------------------------------#
# nss                : number of small patch along strike direction                #
# nds                : number of small patch along dip direction                   #
####################################################################################
nss                  = 3
nds                  = 3
mean                 = 2
std                  = 0.3

# make checker board
chb.make_checkerborad(faultfile, nsegs, ndeps, nss, nds)

# new faultgeom file name
newfaultfile = "FaultGeom_"+str(nss)+"by"+str(nds)

# if the new faultgeom exist
if os.path.exists(newfaultfile) == True:
    # run forward model
    fwd.run_forward(newfaultfile, gpsfile, sarfile, levfile, gfiletype, greentype, nu, wsar)
    # add uncertainty to forward gps data, for GMT2D format
    chb.syn_gps("gps_mod.gmtvec", mean, std)
    newgpsfile = "sync_gps.gmtvec"
    # if new gps file exist
    if os.path.exists(newgpsfile) == True:
        # inversion
        inv.run_inversion(faultfile, nsegs, ndeps, doSubFault, newgpsfile, sarfile, levfile, gfiletype, wsar, greentype, nu, bcs, smoothfactor, bound_switch, bound_dict)
