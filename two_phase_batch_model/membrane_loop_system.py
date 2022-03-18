#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 20:19:21 2022

This is a inital implement of membrane loop system calculations
based on our PBR30 reactor set up with this particular membrane unit

@author: yli5
"""
import numpy as np

class MembraneLoopSystem():
    def __init__(self,permeate_mass_dot,permeate_rou,fis):
        """
        Calculate power requirement for the membrane loop
        Assumptions: each membrane-power relationship scale linearly based on
                     number of membrane units used

        use get_power() to obtain pump power required ( in kW)
        Parameters
        ----------
        permeate_mass_dot : float
            mass flow rate of permeate (sugar stream)
        permeate_rou : float
            density of permeate
        fis : float (range 0 ~ 1)
            insoluble fraction
        """
        # assume membrane area per membrane unit is
        self.area_per_unit = 0.036483011 #m2    # this is the unit on PBR30
        # based on our experiment observation, let's say a permeate flux is 20 LMH
        # this can be manually requested in CEH operation but need to be a
        #realistic number.
        self.p_flux = 20 # permeate flux in LMH
        self.permeate_mass_dot = permeate_mass_dot # permeate mass flow rate in kg/h
        self.permeate_rou = permeate_rou # permeate density g/L
        self.fis = fis # insoluble solids
        self._calc_n_units()
        self._calc_slurry_velocity()
        self._calc_pump_power()

    def _calc_n_units(self):
        permeat_V_dot = self.permeate_mass_dot/self.permeate_rou
        self.total_area = permeat_V_dot/self.p_flux
        self.total_units = int(np.ceil(self.total_area/self.area_per_unit))

    def _calc_slurry_velocity(self):
        self.slurry_velocity = (self.p_flux
                                + 39.1134
                                - 0.078742*pow(self.fis+0.737,-25)
                                )/9.8
    def _calc_pump_power(self):
        if self.slurry_velocity > 3.1:
            self.pump_power_per_unit = 2.1/(1+np.exp(-3.2*self.slurry_velocity+12.1))
        else:
            self.pump_power_per_unit = (0.213278/3.1)*self.slurry_velocity
        self.total_pump_power = self.pump_power_per_unit*self.total_units

    @property
    def pump_power(self):
        return self.total_pump_power

    @property
    def membrane_units(self):
        return self.total_units

    @property
    def permeate_flux(self):
        return self.p_flux

    @permeate_flux.setter
    def permeate_flux(self,value):
        self.p_flux = value
        self._calc_n_units()
        self._calc_slurry_velocity()
        self._calc_pump_power()


    """
    def get_pump_power(self):
        return self.total_pump_power

    def get_membrane_units(self):
        return self.total_units
    """
