#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

import Sim_Initialisation as si
import Crystal_module as Cm

Init = si.Sim_param(2048**2, 5000)

Pump = si.Pulse(1370, 1e10, 100, 0, Init)
Seed = si.Pulse(1500, 1e9, 100,0, Init)
Idler = si.Pulse(10000, 1e10 , 100, 0, Init)
Idler.set_Omega(Pump.get_Omega() - Seed.get_Omega())



medium = Cm.Crystal('AgGaSe2', 1e-6)


plt.plot(Idler.get_t_fs(), Idler.get_E_field_t())
plt.xlim(-400,400)
plt.show()