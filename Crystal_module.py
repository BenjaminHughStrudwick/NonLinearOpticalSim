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
    pmpv2au = 0.51422
    
    def __init__(self, medium, step_size):

        self.medium = medium
        self.step_size = step_size
        
        
    def Sellmeier_cnst(self):
        if self.medium == 'AgGaSe2':
            self.n_Ord = {'a':6.8507,'b':0.4297,'c':0.1584,'d':0.00125}
            self.n_exO = {'a':6.6792,'b':0.4598,'c':0.2122,'d':0.00126}
            self.deff = 33.3*Crystal.pmpv2au

        return (self.n_Ord, self.n_exO, self.deff)

        
             
    def refindexoe(self, Wavelength):
        
        self.Wavelength = Wavelength*1e-3
        self.Sell_const_Or = self.Sellmeier_cnst()[0]
        self.Sell_const_Ex = self.Sellmeier_cnst()[1]

        self.no = np.sqrt(self.Sell_const_Or['a'] + self.Sell_const_Or['b']/(self.Wavelength**2 - self.Sell_const_Or['c']) - self.Sell_const_Or['d']*self.Wavelength**2)
        self.ne = np.sqrt(self.Sell_const_Ex['a'] + self.Sell_const_Ex['b']/(self.Wavelength**2 - self.Sell_const_Ex['c']) - self.Sell_const_Ex['d']*self.Wavelength**2)
        
        return (self.no, self.ne, self.Wavelength*1e3)

    def refindex_angle(self, ne, no, angle_rad):
        self.ne = ne
        self.no = no
        self.angle_rad = angle_rad
        self.na = 1/np.sqrt((np.sin(self.angle_rad)/self.ne)**2 + (np.cos(self.angle_rad)/self.no)**2)
        return self.na

    

