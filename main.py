#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

import Sim_Initialisation as Si
import Crystal_module as Cm

Init = Si.Sim_param(2048**2, 5000)

Pump = Si.Pulse(1e10, 100, 0, Init, 1370)
Seed = Si.Pulse(1e9, 100, 0, Init, 1500)
Idler = Si.Pulse(1e10 , 100, 0, Init)
Idler.set_Pulse_Omega(Pump.get_Pulse_Omega() - Seed.get_Pulse_Omega())



medium = Cm.Crystal('AgGaSe2', 1e-6)
theta = medium.theta
n_pump = medium.refindex_angle(medium.refindexoe(Pump.get_Wavelength()), theta)
n_seed = medium.refindexoe(Seed.get_Wavelength())[0]
n_idler = medium.refindex_angle(medium.refindexoe(Idler.get_Wavelength()), theta)

n_pump_omega = medium.refindex_angle(medium.refindexoe(Init.lam), theta)
n_seed_omega = medium.refindexoe(Init.lam)[0]
n_idler_omega = n_pump_omega

dispersion = medium.get_dispersion(n_pump_omega, n_seed_omega, n_idler_omega, Init.get_d_Omega())




plt.plot(Idler.get_t_fs(), Idler.get_E_field_t())

plt.xlim(-400,400)
plt.show()


