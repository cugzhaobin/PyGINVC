#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Determine the optimal fault geometry using Simulated Annealing(SA) algorithm.
This code depends on functions of PyGINV program and simanneal, which can be 
found in https://github.com/perrygeo/simanneal.

Created on Tue Aug 28 21:07:08 2018

@author: Zhao Bin, Institute of Seismology, CEA.
"""

import numpy as np
import random
import Inverse as inv
import Fault
import LoadData
from simanneal import Annealer
import math, time


class GeomOptimization(Annealer):
    '''
    A simulated annealing class succeed from Annealer.
    '''
    
    def __init__(self, dict_param):
        '''
        Construction the class.
        Input:
            dict_param = dictionary type contains parameters
        '''

        # location of observations and geodetic data
        self.llh_gps    = dict_param['llh_gps']
        self.llh_lev    = dict_param['llh_lev']
        self.llh_sar    = dict_param['llh_sar']
        self.d_gps      = dict_param['d_gps']
        self.d_lev      = dict_param['d_lev']
        self.d_sar      = dict_param['d_sar']
        self.unit       = dict_param['unit']
        self.ndim       = dict_param['ndim']
        # initial fault geometry and origin
        self.nf         = dict_param['nf']
        self.x          = dict_param['x']
        self.origin     = dict_param['origin']
        # Bounds of fault geometry
        self.xlb        = dict_param['xlb']
        self.xub        = dict_param['xub']
        # Green function parameter
        self.greentype  = dict_param['greentype']
        self.nu         = dict_param['nu']
        # data weighting
        self.W          = dict_param['W']
        self.wsar       = dict_param['wsar']
        
        # The initial state is random
        state           = np.zeros((self.nf, 7))
        for i in range(self.nf):
            for j in range(7):
                state[i,j] = random.uniform(self.xlb[i,j], self.xub[i,j])
        # import
        super(GeomOptimization, self).__init__(state)
        
    
    def move(self):
        '''
        Change the state for iteration
        '''
        for i in range(self.nf):
            for j in range(7):
                delt = random.randint(-1,1)*(self.xub[i,j]-self.xlb[i,j])/50.0
                self.state[i,j] += delt
                if self.state[i,j] > self.xub[i,j]:
                    self.state[i,j] = self.xub[i,j]
                if self.state[i,j] < self.xlb[i,j]:
                    self.state[i,j] = self.xlb[i,j]
#                self.state[i,j] = random.uniform(self.xlb[i,j], self.xub[i,j])
                
        
    def energy(self):
        '''
        Compute the cost for each iteration
        '''
        e = inv.disl_wr(self.state, self.llh_gps, self.llh_lev, self.llh_sar, self.unit, self.origin,
                           self.greentype, self.nu, self.ndim, self.d_gps, self.d_lev, self.d_sar, self.W, self.wsar, self.nf)
        return e

    def anneal2(self):
        """Minimizes the energy of a system by simulated annealing.

        Parameters
        state : an initial arrangement of the system

        Returns
        (state, energy): the best state and energy found.
        """
        step = 0
        self.start = time.time()

        # Precompute factor for exponential cooling from Tmax to Tmin
        if self.Tmin <= 0.0:
            raise Exception('Exponential cooling requires a minimum "\
                "temperature greater than zero.')
        Tfactor = 0.9

        # Note initial state
        T = self.Tmax
        E = self.energy()
        prevState = self.copy_state(self.state)
        prevEnergy = E
        self.best_state = self.copy_state(self.state)
        self.best_energy = E
        trials, accepts, improves = 0, 0, 0
        if self.updates > 0:
            updateWavelength = self.steps / self.updates
            self.update(step, T, E, None, None)

        # Attempt moves to new states
        while T>self.Tmin and step < self.steps and not self.user_exit:
            step += 1
            self.move()
            E = self.energy()
            dE = E - prevEnergy
            trials += 1
            if dE > 0.0 and math.exp(-dE / T) < random.random():
                # Restore previous state
                self.state = self.copy_state(prevState)
                E = prevEnergy
            else:
                # Accept new state and compare to best state
                accepts += 1
                if dE < 0.0:
                    improves += 1
                prevState = self.copy_state(self.state)
                prevEnergy = E
                if E < self.best_energy:
                    self.best_state = self.copy_state(self.state)
                    self.best_energy = E
                print(T, dE, math.exp(-dE / T))
            if self.updates > 1:
                if (step // updateWavelength) > ((step - 1) // updateWavelength):
                    self.update(
                        step, T, E, accepts / trials, improves / trials)
                    trials, accepts, improves = 0, 0, 0
            T    *= Tfactor
        self.state = self.copy_state(self.best_state)
        if self.save_state_on_exit:
            self.save_state()

        # Return best state and energy
        return self.best_state, self.best_energy    

    def anneal3(self):
        """Minimizes the energy of a system by simulated annealing.

        Parameters
        state : an initial arrangement of the system

        Returns
        (state, energy): the best state and energy found.
        """
        step = 0
        self.start = time.time()

        # Precompute factor for exponential cooling from Tmax to Tmin
        if self.Tmin <= 0.0:
            raise Exception('Exponential cooling requires a minimum "\
                "temperature greater than zero.')
        Tfactor = 0.98

        # Note initial state
        T = self.Tmax
        E = self.energy()
        prevState = self.copy_state(self.state)
        prevEnergy = E
        self.best_state = self.copy_state(self.state)
        self.best_energy = E
        trials, accepts, improves = 0, 0, 0
        if self.updates > 0:
            updateWavelength = self.steps / self.updates
            self.update(step, T, E, None, None)

        # Attempt moves to new states
        while T > self.Tmin:
            for i in range(5):
                step += 1
                self.move()
                E = self.energy()
                dE = E - prevEnergy
                trials += 1
                if dE > 0.0 and math.exp(-dE / T) < random.random():
                    # Restore previous state
                    self.state = self.copy_state(prevState)
                    E = prevEnergy
                else:
                    # Accept new state and compare to best state
                    accepts += 1
                    if dE < 0.0:
                        improves += 1
                    prevState = self.copy_state(self.state)
                    prevEnergy = E
                    if E < self.best_energy:
                        self.best_state = self.copy_state(self.state)
                        self.best_energy = E
                if self.updates > 1:
                    if (step // updateWavelength) > ((step - 1) // updateWavelength):
                        self.update(
                            step, T, E, accepts / trials, improves / trials)
                        trials, accepts, improves = 0, 0, 0
            self.state = self.copy_state(self.best_state)
            if self.save_state_on_exit:
                self.save_state()
            T *= Tfactor
            print(T)
        # Return best state and energy
        return self.best_state, self.best_energy  

if __name__ == '__main__':

    faultfile  = 'faultgeom_all'
    faultbound = 'faultbounds_all_2'
#   faultfile  = 'FaultGeom'
#   faultbound = 'FaultBounds'
    gpsfile    = './comb_3d.gmtvec'
    sarfile    = ''
    levfile    = ''
    gfiletype  = 'GMT3D'
    nu         = 0.25
    greentype  = [1,0,0]
    
    ############################################################################
    #                    Load Fault from files                                 #
    ############################################################################
    [dis_geom_grid, geom_grid, origin] = Fault.LoadFault(faultfile, 1, 1, False)
    x0 = Fault.geom2x(geom_grid, origin)
    nf = len(geom_grid)

    ############################################################################
    #                    Load Fault from files                                 #
    ############################################################################
    [xlb, xub] = Fault.GetFaultBounds(faultbound, nf)
    print(xlb)
    print(xub)

    ############################################################################
    #                    Load Data from files                                  #
    ############################################################################
    [llh_gps, llh_lev, llh_sar, d_gps, d_lev, d_sar, unit, W] = \
    LoadData.LoadAllData(gpsfile, sarfile, levfile, gfiletype, origin)

    ndim = int(gfiletype[3])
    dict_param = {'x'       : x0,
                  'llh_gps' : llh_gps,
                  'llh_sar' : llh_sar,
                  'llh_lev' : llh_lev,
                  'd_gps'   : d_gps,
                  'd_sar'   : d_sar,
                  'd_lev'   : d_lev,
                  'nu'      : nu,
                  'nf'      : len(geom_grid),
                  'origin'  : origin,
                  'ndim'    : ndim,
                  'unit'    : unit,
                  'W'       : W,
                  'wsar'    : [1],
                  'xlb'     : xlb,
                  'xub'     : xub,
                  'greentype': greentype}
    
    geomopt = GeomOptimization(dict_param)
    geomopt.steps = 1000
    geomopt.copy_strategy ='deepcopy'
    geomopt.Tmax = 15000
    geomopt.Tmin = 1e-1
    state, wrms = geomopt.anneal3()
    print(state)
    print(wrms)
    a,b  = Fault.x2geom(state, origin)
    slip = np.zeros((4,3))
    a    = np.hstack((a, slip))
    np.savetxt('faultgeom_optim', a, fmt="%f %f %f %f %f %f %f %f %f %f")
#    print b
