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
    def __init__(self,permeate_mass_dot,permeate_rou,fis,system_name:str = "PBR30", n_units: int=20,
                tube_id=0.5*0.0254,tube_l=36*0.0254,p_flux=20):
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
        n_units: int
            number of membrane units operated in serial.
        tube_id and tube_l are membrane module dimensions. default value is for CEH_PBR30
        p_flux: float
            membrane permeate flux , unit (LMH)
        """
        self.n_units = n_units
        self.p_flux = p_flux # permeate flux in LMH
        self.permeate_mass_dot = permeate_mass_dot # permeate mass flow rate in kg/h
        self.permeate_rou = permeate_rou # permeate density g/L
        self.fis = fis # insoluble solids

        if system_name == "PBR30":
            self.tube_id = 0.5*0.0254
            self.tube_l = 36*0.0254
        elif system_name == "Pall":
            self.tube_id = 0.5*0.0254
            self.tube_l = 6*12*0.0254
        else:
            print(f"System name: {system_name} not found! Choose from: PBR30 or Pall. Otherwise, provide tube_id and tube_l")
            print(f"Using cutomized membrane system with ID: {tube_id} and Length: {tube_l}.")
            self.tube_id = tube_id
            self.tube_l = tube_l

        self.permeance = 50 # LMH/bar
        self.area_per_unit = np.pi * self.tube_id * self.tube_l #m2
        self.area_n_units = self.n_units * self.area_per_unit
        # based on our experiment observation, let's say a permeate flux is 20 LMH
        # this can be manually requested in CEH operation but need to be a
        #realistic number.
        self._calc_n_units()
        self._calc_slurry_velocity()
        self._calc_retentate_permeate_ratio()
        self._calc_pump_power()

    def _calc_n_units(self):
        permeat_V_dot = self.permeate_mass_dot/self.permeate_rou
        self.total_area = permeat_V_dot/self.p_flux
        self.total_units = int(np.ceil(self.total_area/self.area_per_unit))

    def _calc_slurry_velocity(self):
        self.slurry_velocity = (self.p_flux
                                + 26.692
                                - 0.166*pow(self.fis+0.737,-25)
                                )/9.8

    def _calc_retentate_permeate_ratio(self):
        self.retentate_flow_rate = np.pi * pow(self.tube_id / 2, 2) * self.slurry_velocity * 1000 * 3600 # kg/h
        self.permeate_flow_rate = self.p_flux * self.area_n_units
        self.retentate_permeate_ratio = self.retentate_flow_rate/self.permeate_flow_rate
        self.total_retentate_flow_rate = self.retentate_flow_rate * self.total_units / self.n_units

    def _calc_pump_power(self):
        self.pump_power_per_unit = 2.1/(1+np.exp(-3.2*self.slurry_velocity+12.1))
        self.total_pump_power = self.pump_power_per_unit*self.total_units

    @property
    def total_membrane_area(self):
        return self.total_area

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
        self._calc_retentate_permeate_ratio()
        self._calc_pump_power()


    """
    def get_pump_power(self):
        return self.total_pump_power

    def get_membrane_units(self):
        return self.total_units
    """
