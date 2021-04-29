#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Simulation Paramter Class
# =============================================================================

import numpy as np
from scipy.constants import physical_constants

class Sim_param:
    
    
    
    def __init__(self, N, t_max):
        self.N = N
        self.t_max = t_max
        
        self.fs2au = 1e-15/physical_constants['atomic unit of time'][0]
        self.dt = self.fs2au*(2*self.t_max)/self.N
        self.t = (np.arange(-1*self.t_max, self.t_max, self.dt))*self.fs2au
        
    def get_Precision(self):
        return self.N
    
    def set_Precision(self, N):
        self.N = N
        return self.N
    
    def get_tmax(self):
        return self.t_max
    
    def set_tmax(self, t_max):
        self.t_max = t_max
        return self.t_max
    
    def get_t(self):
        return self.t
    
    def get_dt(self):
        return self.dt
    
    def get_fs2au(self):
        return self.fs2au


class Pulse(Sim_param):

    def __init__(self, Wavelength, Intensity, Duration, Delay, simparam):
        
        self.N = simparam.N
        self.t_max = simparam.t_max
        self.fs2au = simparam.fs2au
        self.dt = simparam.dt
        self.t = simparam.t


        self.Wavelength = Wavelength
        self.Intensity = Intensity
        self.Duration = Duration*self.fs2au
        self.Delay = Delay*self.fs2au
        
        self.Omega = (0.057*800)/self.Wavelength
        self.E_field = np.sqrt(self.Intensity/3.509e16)

    def get_Wavelength(self):
        return self.Wavelength

    def set_Wavelength(self, Wavelength):
        self.Wavelength = Wavelength
        self.Omega = (0.057*800)/self.Wavelength
        return self.Wavelength

    def get_Omega(self):
        return self.Omega

    def set_Omega(self, Omega):
        self.Omega = Omega
        self.Wavelength = (0.057*800)/self.Omega
        return self.Omega

    def get_t_fs(self):
        return self.t/self.fs2au

    def get_Delay_fs(self):
        return self.Delay/self.fs2au

    def get_Delay_au(self):
        return self.Delay

    def set_Delay(self, Delay):
        self.Delay = Delay*self.fs2au
        return self.Delay



    def get_E_field_t(self):
        self.E_field_t = self.E_field*np.exp(-2*np.log(2)*(self.t - self.Delay)**(2)/self.Duration**2)*np.exp(1j*self.Omega*(self.t - self.Delay))
        return self.E_field_t

