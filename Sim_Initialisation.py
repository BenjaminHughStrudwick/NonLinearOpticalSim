#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Simulation Paramter Class
# =============================================================================

import numpy as np
from scipy.constants import physical_constants

class Sim_param:
    
    fs2au = 1e-15/physical_constants['atomic unit of time'][0]
    
    def __init__(self, N, t_max):
        self.N = N
        self.t_max = t_max
        

        self.dt = Sim_param.fs2au*(2*self.t_max)/self.N
        self.t = (np.arange(-1*self.t_max, self.t_max, self.dt))*Sim_param.fs2au

        self.Omega = 2*np.pi*np.arange(-1*self.N/2 +1, self.N/2 + 1)/self.N/self.dt
        self.omega_shifted = np.fft.fftshift(self.Omega)
        self.lam = np.divide((0.057*800), self.Omega, out=np.zeros_like(self.Omega), where=self.Omega!=0)
        self.lam[int(np.where(self.lam == 0)[0])] = np.inf #Corrects previous line x / 0 = inf
        self.d_Omega = self.Omega[1] - self.Omega[0]

        
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
        return Sim_param.fs2au
    
    def get_d_Omega(self):
        return self.d_Omega


class Pulse(Sim_param):

    
    def __init__(self, Intensity, Duration, Delay, simparam, Wavelength=None):
        
        self.N = simparam.N
        self.t_max = simparam.t_max
        self.dt = simparam.dt
        self.t = simparam.t
        self.Intensity = Intensity
        self.Duration = Duration*Sim_param.fs2au
        self.Delay = Delay*self.fs2au
        self.Wavelength = Wavelength if Wavelength is not None else 1
        
        self.Pulse_Omega = (0.057*800)/self.Wavelength
        self.E_field = np.sqrt(self.Intensity/3.509e16)

    def get_Wavelength(self):
        return self.Wavelength

    def set_Wavelength(self, Wavelength):
        self.Wavelength = Wavelength
        self.Pulse_Omega = (0.057*800)/self.Wavelength
        return self.Wavelength

    def get_Pulse_Omega(self):
        return self.Pulse_Omega

    def set_Pulse_Omega(self, Pulse_Omega):
        self.Pulse_Omega = Pulse_Omega
        self.Wavelength = (0.057*800)/self.Pulse_Omega
        return self.Pulse_Omega

    def get_t_fs(self):
        return self.t/Sim_param.fs2au

    def get_Delay_fs(self):
        return self.Delay/Sim_param.fs2au

    def get_Delay_au(self):
        return self.Delay

    def set_Delay(self, Delay):
        self.Delay = Delay*Sim_param.fs2au
        return self.Delay

    def get_lamda(self):
        return self.lam
    
    def get_dOmega(self):
        return self.d_Omega


    def get_E_field_t(self):
        self.E_field_t = self.E_field*np.exp(-2*np.log(2)*(self.t - self.Delay)**(2)/self.Duration**2)*np.exp(1j*self.Pulse_Omega*(self.t - self.Delay))
        return self.E_field_t

