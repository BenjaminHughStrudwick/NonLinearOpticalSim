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
