#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Crystal Paramters and Functions
# =============================================================================

import scipy.constants as cnst


class Crystal:
    
    def __init__(self, medium, step_size):

        self.medium = medium
        self.step_size = step_size
        self.au2mm = 1e3*cnst.physical_constants['atomic unit of length'][0] #5.292e-8  conversion of mm length of crystal to au
        self.c_au = cnst.c/cnst.physical_constants['atomic unit of velocity'][0] # speed of light in atomic units
        
    def Sellmeier_cnst(self):
        if self.medium is 'AgGaSe2':
            self.n_Ord = {'a':6.8507,'b':0.4297,'c':0.1584,'d':0.00125}
            self.n_exO = {'a':6.6792,'b':0.4598,'c':0.2122,'d':0.00126}
            self.pmpv2au=0.51422
            self.deff = 33.3*self.pmpv2au
        else:
            print('Please Input a valid Crystal code')
        return (self.n_Ord, self.n_exO, self.deff)
             
