#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:53:07 2020

@author: johngillan
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font = {'size':10}
matplotlib.rc('font', **font)

G = 6.674e-11                               # gravitational constanst
au = 1.496e11                               # astronomical unit
year = 365.25*24.*60.*60.                   # year
mstar = 1.989e30                              # mass of star
mplanet = 6e23                            # mass of orbiting planet
r = 1*au   # radial position of planet

h = 0.05                                    # scale height of disc

# switch these:
i = 0
e = np.arange(0, 0.5, 0.005)                # eccentricites
# e = 0
# i = np.arange(0, h*10, 0.0001)              # inclinations

ehat = e/h                                  # eccentricity divided by scale height of disc
ihat = i/h                                  # inclinations divided by scale height of disc

omegak = np.sqrt(G*(mstar+mplanet)/r**3)                # Keplerian frequency
sigma = 17000*(r/au)**(-3/2)                              # surface density of disc
t_wave = mstar/mplanet*mstar/sigma/r**2*h**(4)/omegak   # equation 7 from Ida 2020

# equations for tau_e, tau_a, and tau_m by Cresswell & Nelson 2008 as taken from Ida 2020
t_e_CN = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3+0.18*ehat*ihat**2)
t_i_CN = t_wave/0.544*(1-0.30*ihat**2+0.24*ihat**3+0.14*ihat*ehat**2)
Pe = (1+(ehat/2.25)**1.2+(ehat/2.84)**6)/(1-(ehat/2.02)**4)
t_m_CN = t_wave*2/h**2/(2.7+1.1*0.5)*(Pe+Pe/np.abs(Pe)*(0.07*ihat+0.085*ihat**4-0.08*ehat*ihat**2))
t_a_CN = (2/t_m_CN+2*e**2/t_e_CN)**(-1)                                                          # equation 16 from Ida 2020    

Cm = 21.
Ct = 4.25

# equations for tau_e, tau_a, tau_m and tau_i by Ida 2020
t_e = t_wave/0.780*(1+1/15*(ehat**2+ihat**2)**(3/2))
t_i = t_wave/0.544*(1+1/21.5*(ehat**2+ihat**2)**(3/2))
t_a = t_wave/Ct/h**2*(1+Ct/Cm*(ehat**2+ihat**2)**(1/2))

fig, ax = plt.subplots(2, figsize=(7,7))

# switch these:
x = ehat
# x = ihat

ax[0].plot(x, t_e/year, label='IDA20')
ax[0].plot(x, t_e_CN/year, label='CN08')
# ax[0].plot(x, tau_e/year)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel(r'$t_e$ (years)')
# ax[0].set_ylabel(r'$\tau^{-1}_e/t_{wave}^{-1}$')
ax[0].set_xlim(0.1, 10)
# ax[0].set_ylim(0.001, 10)
ax[0].tick_params(which='both', direction="in", top=True, right=True)
ax[0].legend()

ax[1].plot(x, t_a/year, label='IDA20')
ax[1].plot(x, t_a_CN/year, label='CN08')
# ax[1].plot(x, t_a/year)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel('e/h')
ax[1].set_ylabel(r'$\tau_a$ (years)')
# ax[0].set_ylabel(r'$\tau^{-1}_a/t_{wave}^{-1}$')
ax[1].set_xlim(0.1, 10)
# ax[1].set_ylim(10e5)
ax[1].tick_params(which='both', direction="in", top=True, right=True)
ax[1].legend()

# fig.savefig('/home/john/Desktop/summerproject/img/i_timescales10.png', bbox_inches='tight')