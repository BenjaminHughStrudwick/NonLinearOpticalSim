#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Crystal Paramters and Functions
# =============================================================================

import scipy.constants as cnst
import numpy as np



class Crystal:

    au2mm = 1e3*cnst.physical_constants['atomic unit of length'][0] #5.292e-8  conversion of mm length of crystal to au
    c_au = cnst.c/cnst.physical_constants['atomic unit of velocity'][0] # speed of light in atomic units
    fs2au = 1e-15/cnst.physical_constants['atomic unit of time'][0]
    pmpv2au = 0.51422
    theta = 52.9
    
    def __init__(self, medium, step_size):

        self.medium = medium
        self.step_size = step_size
        
        
    def Sellmeier_cnst(self):
        if self.medium == 'AgGaSe2':
            self.n_Ord = {'a':6.8507,'b':0.4297,'c':0.1584,'d':0.00125}
            self.n_exO = {'a':6.6792,'b':0.4598,'c':0.2122,'d':0.00126}
            

        return (self.n_Ord, self.n_exO)

    def get_deff(self):
        self.deff = 33.3*Crystal.pmpv2au
        return self.deff

        
             
    def refindexoe(self, Wavelength):
        
        self.Wavelength = Wavelength*1e-3
        self.Sell_const_Or = self.Sellmeier_cnst()[0]
        self.Sell_const_Ex = self.Sellmeier_cnst()[1]

        self.no = np.sqrt(self.Sell_const_Or['a'] + self.Sell_const_Or['b']/(self.Wavelength**2 - self.Sell_const_Or['c']) - self.Sell_const_Or['d']*self.Wavelength**2)
        self.ne = np.sqrt(self.Sell_const_Ex['a'] + self.Sell_const_Ex['b']/(self.Wavelength**2 - self.Sell_const_Ex['c']) - self.Sell_const_Ex['d']*self.Wavelength**2)
        
        return (self.no, self.ne)

    def refindex_angle(self, nOnEx, theta):
        self.ne = nOnEx[1]
        self.no = nOnEx[0]
        self.angle_rad = np.deg2rad(Crystal.theta)
        self.na = 1/np.sqrt((np.sin(self.angle_rad)/self.ne)**2 + (np.cos(self.angle_rad)/self.no)**2)
        return self.na

    def get_theta(self):
        return self.theta

    def get_dispersion(self, n_pump_omega, n_seed_omega, n_idler_omega, d_omega):
        self.n_omega_list = [n_pump_omega, n_seed_omega, n_idler_omega]
        self.n_omega_keys = ['pump', 'seed', 'idler']
        self.d_omega = d_omega
        self.dispersion_dic = {}
        d = 0
        for i in self.n_omega_list:
            
            self.dndo = np.diff(i)/self.d_omega
            self.dndo = np.insert(self.dndo, 0, self.dndo[0], axis=0)
            self.dn2do = np.diff(i, 2)/self.d_omega**2
            self.dn2do = np.insert(self.dn2do, 0, [self.dn2do[0], self.dn2do[0]] , axis=0)

            self.dispersion_dic[self.n_omega_keys[d]] = (self.dndo, self.dn2do)
            d+=1

        return self.dispersion_dic

    def get_gvd(self, n_pump_omega, n_seed_omega, n_idler_omega, dispersion):
        self.n_omega_list = [n_pump_omega, n_seed_omega, n_idler_omega]
        self.n_omega_keys = ['pump', 'seed', 'idler']
        self.dispersion = dispersion
        self.gvd_dic = {}

        d = 0
        for i in self.n_omega_list:
            
            
            self.gvd_dic[self.n_omega_keys[d]] = (2/Crystal.c_au*dispersion[self.n_omega_keys[d]][0]+(i/Crystal.c_au)*dispersion[self.n_omega_keys[d]][1])
            d+=1

        return self.gvd_dic

    def get_gvd_fs2pmm(self, gvd):
        self.n_omega_keys = ['pump', 'seed', 'idler']
        self.gvd = gvd
        self.gvd_fs2pmm_dic = {}
        for i in self.n_omega_keys:
            self.gvd_fs2pmm_dic[i] = (self.gvd[i]/self.fs2au**2)/Crystal.au2mm
        
        return self.gvd_fs2pmm_dic

    def get_dk(self, Wp, Np, Ws, Ns, Wi, Ni):
        self.Wp = Wp
        self.Np = Np
        self.Ws = Ws
        self.Ns = Ns
        self.Wi = Wi
        self.Ni = Ni

        self.dk = 2*np.pi*(self.Np/self.Wp - self.Ns/self.Ws - self.Ni/self.Wi)
        self.dk=self.dk*1e6*Crystal.au2mm #Change to au

        return self.dk

    def get_c_au(self):
        return Crystal.c_au

    def get_stepsize_au(self):
        self.step_size = self.step_size*1e3*self.au2mm
        return self.step_size
    
    def get_au2mm(self):
        return Crystal.au2mm

        

        



