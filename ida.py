#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:53:07 2020

@author: john
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font = {'size':10}
matplotlib.rc('font', **font)

G = 6.674e-11                               # gravitational constanst
au = 1.496e11                               # astronomical unit
year = 365.25*24.*60.*60.                   # year
mstar = 2.0**30                             # mass of star
mplanet = 1.0**24                           # mass of orbiting planet
r = 10*au                                   # radial position of planet
omegak = np.sqrt(G*(mstar+mplanet)/r**3)    # Keplerian frequency
h = 0.05                                    # scale height of disc
sigma = 17000*(r/au)**(-1)*mstar        # surface density of disc
e = np.arange(0, 0.5, 0.0001)                 # eccentricites
ehat = e/h                                  # eccentricity divided by scale height of disc

# equation 7 from Ida 2020
t_wave = 1/((mplanet/mstar)*(sigma*r**2/mstar)*h**(-4)*omegak)

# equations for tau_e, tau_a, and tau_m by Cresswell & Nelson 2008 as taken from Ida 2020
tau_e_CN = 1/(0.78*(1-0.14*ehat**2+0.06*ehat**3)**(-1)/t_wave)                       # equation 17 from Ida 2020
tau_m_CN = 1/((2.7+1.1*0.5)/(2)*h**2*(1-(ehat/2.02)**4)/(1+(ehat/2.25)**0.5+(ehat/2.84)**6)/t_wave) # equation 18 from Ida 2020
tau_a_CN = 1/(2/tau_m_CN+2*e**2/tau_e_CN)                                   # equation 16 from Ida 2020




# p = -dln(sigma)/dln(r), q = -dln(T)/dln(r)
p_red, q_red = 1.0, 0.5
# equations for Cm and Ct taken from table 1 in Ida 2020
Ct_red = 2.73+1.08*p_red+0.87*q_red
Cm_red = 6*(2*p_red-q_red+2)

# magenta line in plot using different values of p and q
p_mag, q_mag = 0.5, 1.0
Ct_mag = 2.73+1.08*p_mag+0.87*q_mag
Cm_mag = 6*(2*p_mag-q_mag+2)

# equations for tau_e, tau_a, and tau_m by Ida 2020
tau_e = 1/(0.780*(1+(1/15)*ehat**3)**(-1)/t_wave)               # equation 34 from Ida 2020
tau_a_red = 1/(Ct_red*h**2*(1+(Ct_red/Cm_red)*ehat)**(-1)/t_wave)   # equation 38 from Ida 2020
tau_m_red = 1/(0.5/tau_a_red - e**2/tau_e)           # equation 39 from Ida 2020
# magenta lines in plot
tau_a_mag = 1/(Ct_mag*h**2*(1+(Ct_mag/Cm_mag)*ehat)**(-1)/t_wave)   # equation 38 from Ida 2020
tau_m_mag = 1/(0.5/tau_a_mag - e**2/tau_e)           # equation 39 from Ida 2020



fig, ax = plt.subplots(3, figsize=(6,12))

# figure 1a from Ida 2020 - tau_e as a function of e/h
ax[0].plot(ehat, 1/tau_e_CN*t_wave, label='CN08', c='blue')
ax[0].plot(ehat, 1/tau_e*t_wave, label='IDA20', c='red')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel(r'$\tau^{-1}_e/t_{wave}^{-1}$')
ax[0].set_xlim(0.1, 10)
ax[0].set_ylim(1e-3, 1e1)
ax[0].tick_params(which='both', direction="in", top=True, right=True)
ax[0].legend()

# figure 1b from Ida 2020 - tau_m as a function of e/h
ax[1].plot(ehat, 1/tau_m_CN*t_wave, label='CN08', c='blue')
ax[1].plot(ehat, -1/tau_m_CN*t_wave, c='blue', ls='--')
ax[1].plot(ehat, 1/tau_m_red*t_wave, label='IDA20 p=1 q=0.5', c='red')
ax[1].plot(ehat, -1/tau_m_red*t_wave, c='red', ls='--')
ax[1].plot(ehat, 1/tau_m_mag*t_wave, label='IDA20 p=0.5 q=1', c='magenta')
ax[1].plot(ehat, -1/tau_m_mag*t_wave, c='magenta', ls='--')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_ylabel(r'$\tau^{-1}_m/t_{wave}^{-1}$')
ax[1].set_xlim(0.1, 10)
ax[1].set_ylim(1e-4, 1e-1)
ax[1].tick_params(which='both', direction="in", top=True, right=True)
ax[1].legend()

# figure 1c from Ida 2020 - tau_a as a function of e/h
ax[2].plot(ehat, 1/tau_a_CN*t_wave, label='CN08', c='blue')
ax[2].plot(ehat, 1/tau_a_red*t_wave, label='IDA20 p=1 q=0.5', c='red')
ax[2].plot(ehat, 1/tau_a_mag*t_wave, label='IDA20 p=0.5 q=1', c='magenta')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
ax[2].set_xlabel('e/h')
ax[2].set_ylabel(r'$\tau^{-1}_a/t_{wave}^{-1}$')
ax[2].set_xlim(0.1, 10)
ax[2].set_ylim(1e-4, 1e-1)
ax[2].tick_params(which='both', direction="in", top=True, right=True)
ax[2].legend()

# fig.savefig('/home/john/Desktop/summerproject/img/ida2020_fig1.png', bbox_inches='tight')

# %%

fig, ax = plt.subplots(3, figsize=(6,12))

# tau_e as a function of e/h
ax[0].plot(ehat, tau_e_CN*year, label='CN08', c='blue')
ax[0].plot(ehat, tau_e*year, label='IDA20', c='red')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel(r'$\tau_e$ (yr)')
ax[0].set_xlim(0.1, 10)
ax[0].set_ylim(1e2, 1e5)
ax[0].tick_params(which='both', direction="in", top=True, right=True)
ax[0].legend()

# tau_m as a function of e/h
ax[1].plot(ehat, tau_m_CN*year, label='CN08', c='blue')
ax[1].plot(ehat, -tau_m_CN*year, c='blue', ls='--')
ax[1].plot(ehat, tau_m_red*year, label='IDA20 p=1 q=0.5', c='red')
ax[1].plot(ehat, -tau_m_red*year, c='red', ls='--')
ax[1].plot(ehat, tau_m_mag*year, label='IDA20 p=0.5 q=1', c='magenta')
ax[1].plot(ehat, -tau_m_mag*year, c='magenta', ls='--')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_ylabel(r'$\tau_m$ (yr)')
ax[1].set_xlim(0.1, 10)
# ax[1].set_ylim(1e-4, 1e-1)
ax[1].tick_params(which='both', direction="in", top=True, right=True)
ax[1].legend()

# tau_a as a function of e/h
ax[2].plot(ehat, tau_a_CN*year, label='CN08', c='blue')
ax[2].plot(ehat, tau_a_red*year, label='IDA20 p=1 q=0.5', c='red')
ax[2].plot(ehat, tau_a_mag*year, label='IDA20 p=0.5 q=1', c='magenta')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
ax[2].set_xlabel('e/h')
ax[2].set_ylabel(r'$\tau_a$ (yr)')
ax[2].set_xlim(0.1, 10)
ax[2].set_ylim(1e4, 1e6)
ax[2].tick_params(which='both', direction="in", top=True, right=True)
ax[2].legend()

# fig.savefig('/home/john/Desktop/summerproject/img/tau_in_years.png', bbox_inches='tight')