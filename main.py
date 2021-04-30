#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

import Sim_Initialisation as Si
import Crystal_module as Cm

Init = Si.Sim_param(2048**2, 5000)

Pump = Si.Pulse(1e10, 100, 0, Init, 1370)
Seed = Si.Pulse(1e9, 100, 0, Init, 1450)
Idler = Si.Pulse(1e10 , 100, 0, Init)
Idler.set_Pulse_Omega(Pump.get_Pulse_Omega() - Seed.get_Pulse_Omega())



medium = Cm.Crystal('AgGaSe2', 1e-6)



#plt.plot(Idler.get_t_fs(), Idler.get_E_field_t())
#plt.xlim(-400,400)
#plt.show()

print(medium.Sellmeier_cnst())
