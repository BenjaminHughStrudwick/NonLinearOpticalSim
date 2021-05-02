#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import Sim_Initialisation as Si
import Crystal_module as Cm

# =============================================================================
# Initilising the Simulation and creating the Pulses
# =============================================================================

Init = Si.Sim_param(2048**2, 5000)

Pump = Si.Pulse(1e10, 100, 0, Init, 1370)
Seed = Si.Pulse(1e9, 100, 0, Init, 1500)
Idler = Si.Pulse(0 , 100, 0, Init)
Idler.set_Pulse_Omega(Pump.get_Pulse_Omega() - Seed.get_Pulse_Omega())


# =============================================================================
# Creating the Nonlinear Crystal medium and parameters associated with each pulse
# =============================================================================

medium = Cm.Crystal('AgGaSe2', 1e-6)
theta = medium.theta
n_pump = medium.refindex_angle(medium.refindexoe(Pump.get_Wavelength()), theta)
n_seed = medium.refindexoe(Seed.get_Wavelength())[0]
n_idler = medium.refindex_angle(medium.refindexoe(Idler.get_Wavelength()), theta)

n_pump_omega = medium.refindex_angle(medium.refindexoe(Init.lam), theta)
n_seed_omega = medium.refindexoe(Init.lam)[0]
n_idler_omega = n_pump_omega

dk = medium.get_dk(Pump.get_Wavelength(), n_pump, Seed.get_Wavelength(), n_seed, Idler.get_Wavelength(), n_idler)

dispersion = medium.get_dispersion(n_pump_omega, n_seed_omega, n_idler_omega, Init.get_d_Omega())
gvd = medium.get_gvd(n_pump_omega, n_seed_omega, n_idler_omega,dispersion)
gvd_fs2pmm = medium.get_gvd_fs2pmm(gvd)


# =============================================================================
# Engine - numerically solves nonlinear partial differential eqn - Mixes the pulses 
# =============================================================================

I_peak_ef1 = [0,0]
I_peak_ef_dfg = [0,0]
I_peak_ef3 = [0,0]
z_au2mm = [0,0]

a_s = (8*Seed.get_Pulse_Omega()*medium.get_deff())/(medium.get_c_au()*n_seed)
a_i = (8*Idler.get_Pulse_Omega()*medium.get_deff())/(medium.get_c_au()*n_idler)
a_p = (8*Pump.get_Pulse_Omega()*medium.get_deff())/(medium.get_c_au()*n_pump)

dz = medium.get_stepsize_au()
n = 0

ef1 = Seed.get_E_field_t()
ef3 = Pump.get_E_field_t()
ef_dfg = Idler.get_E_field_t()



while(np.diff(I_peak_ef1)[-1] >= 0):
    n += 1
    z=n*dz
# =============================================================================
# Step 1, half step size, linear part in the frequyency domain
# =============================================================================
    fe1 = np.fft.fftshift(np.fft.fft(ef1))
    fe1 = fe1*np.exp(-1j*0.5*dz*dispersion['seed'][0][Seed.get_idx()]*(Init.Omega - Seed.get_Pulse_Omega())-0.5j*gvd['seed'][Seed.get_idx()]*0.5*dz*(Init.Omega-Seed.get_Pulse_Omega())**2) 
    ef1 = np.fft.ifft(np.fft.fftshift(fe1))

    fe3 = np.fft.fftshift(np.fft.fft(ef3))
    fe3 = fe3*np.exp(-1j*0.5*dz*dispersion['pump'][0][Pump.get_idx()]*(Init.Omega-Pump.get_Pulse_Omega())-0.5j*gvd['pump'][Pump.get_idx()]*0.5*dz*(Init.Omega-Pump.get_Pulse_Omega())**2) 
    ef3 = np.fft.ifft(np.fft.fftshift(fe3))

    fe_dfg = np.fft.fftshift(np.fft.fft(ef_dfg))
    fe_dfg = fe_dfg*np.exp(-1j*0.5*dz*dispersion['idler'][0][Idler.get_idx()]*(Init.Omega-Idler.get_Pulse_Omega())-0.5j*gvd['idler'][Idler.get_idx()]*0.5*dz*(Init.Omega-Idler.get_Pulse_Omega())**2) 
    ef_dfg = np.fft.ifft(np.fft.fftshift(fe_dfg))
    
# =============================================================================
#     Step 2, full step size, nonlinear part, fouth-order runge-kutta method
# =============================================================================
    z1 = (n+0.5)*dz
    z2 = (n+1)*dz

    k1 = 1j*a_s*ef3*np.conj(ef_dfg)*np.exp(1j*dk*z)
    h1 = 1j*a_i*ef3*np.conj(ef1)*np.exp(1j*dk*z)
    g1 = 1j*a_p*ef1*ef_dfg*np.exp(-1j*dk*z)
    
    k2 = 1j*a_s*(ef3+0.5*g1*dz)*np.conj(ef_dfg+0.5*h1*dz)*np.exp(1j*dk*z1)
    h2 = 1j*a_i*(ef3+0.5*g1*dz)*np.conj(ef1+0.5*k1*dz)*np.exp(1j*dk*z1)
    g2 = 1j*a_p*(ef1+0.5*k1*dz)*(ef_dfg+0.5*h1*dz)*np.exp(-1j*dk*z1)
    
    k3 = 1j*a_s*(ef3+0.5*g2*dz)*np.conj(ef_dfg+0.5*h2*dz)*np.exp(1j*dk*z1)
    h3 = 1j*a_i*(ef3+0.5*g2*dz)*np.conj(ef1+0.5*k2*dz)*np.exp(1j*dk*z1)
    g3 = 1j*a_p*(ef1+0.5*k2*dz)*(ef_dfg+0.5*h2*dz)*np.exp(-1j*dk*z1)
    
    k4 = 1j*a_s*(ef3+g3*dz)*np.conj(ef_dfg+h3*dz)*np.exp(1j*dk*z2)
    h4 = 1j*a_i*(ef3+g3*dz)*np.conj(ef1+k3*dz)*np.exp(1j*dk*z2)
    g4 = 1j*a_p*(ef1+k3*dz)*(ef_dfg+h3*dz)*np.exp(-1j*dk*z2)
    
    ef1 = ef1+dz*(k1+2*k2+2*k3+k4)/6
    ef_dfg = ef_dfg+dz*(h1+2*h2+2*h3+h4)/6
    ef3 = ef3+dz*(g1+2*g2+2*g3+g4)/6
    
# =============================================================================
#    Step 3, half step size, linear part in the frequyency domain 
# =============================================================================
    fe1 = np.fft.fftshift(np.fft.fft(ef1))
    fe1 = fe1*np.exp(-1j*0.5*dz*dispersion['seed'][0][Seed.get_idx()]*(Init.Omega-Seed.get_Pulse_Omega())-0.5j*gvd['seed'][Seed.get_idx()]*0.5*dz*(Init.Omega-Seed.get_Pulse_Omega())**2) 
    ef1 = np.fft.ifft(np.fft.fftshift(fe1))
    
    fe3 = np.fft.fftshift(np.fft.fft(ef3))
    fe3 = fe3*np.exp(-1j*0.5*dz*dispersion['pump'][0][Pump.get_idx()]*(Init.Omega-Pump.get_Pulse_Omega())-0.5j*gvd['pump'][Pump.get_idx()]*0.5*dz*(Init.Omega-Pump.get_Pulse_Omega())**2) 
    ef3 = np.fft.ifft(np.fft.fftshift(fe3))
    
    fe_dfg = np.fft.fftshift(np.fft.fft(ef_dfg))
    fe_dfg = fe_dfg*np.exp(-1j*0.5*dz*dispersion['idler'][0][Idler.get_idx()]*(Init.Omega-Idler.get_Pulse_Omega())-0.5j*gvd['idler'][Idler.get_idx()]*0.5*dz*(Init.Omega-Pump.get_Pulse_Omega())**2) 
    ef_dfg = np.fft.ifft(np.fft.fftshift(fe_dfg)) 
    
# =============================================================================
#    Appends the Energy to the intensty list         
# =============================================================================

    
    I_peak_ef1.append(float(sum(np.real(ef1)**2)))
    I_peak_ef_dfg.append(float(sum(np.real(ef_dfg)**2)))
    I_peak_ef3.append(float(sum(np.real(ef3)**2)))
#    z_au2mm.append(z*medium.get_au2mm)



fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Electric Field of Idler before and after Mixing Process')

ax1.plot(Idler.get_t_fs(), Idler.get_E_field_t(), c='g')
#ax1.plot(Pump.get_t_fs(), Pump.get_E_field_t(), c='r')
#ax1.plot(Seed.get_t_fs(), Seed.get_E_field_t(),c='b')

ax2.plot(Idler.get_t_fs(), ef_dfg, c='g')
#ax2.plot(Pump.get_t_fs(), ef3, c='r')
#ax2.plot(Seed.get_t_fs(), ef1,c='b')


ax1.set_xlim(-400,400)
ax2.set_xlim(-400,400)
plt.show()


