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
mstar = 2.0e30                              # mass of star
mplanet = 5.0e24                            # mass of orbiting planet
r = 1*au   # radial position of planet

h = 0.05                                    # scale height of disc
e = np.arange(0, 0.5, 0.005)                # eccentricites
# e = 0.5
ehat = e/h                                  # eccentricity divided by scale height of disc
i = np.arange(0, h*10, 0.0001)              # inclinations
i = 0
ihat = i/h                                  # inclinations divided by scale height of disc

omegak = np.sqrt(G*(mstar+mplanet)/r**3)                # Keplerian frequency
sigma = 17000*(r/au)**(-3/2)                              # surface density of disc
t_wave = mstar/mplanet*mstar/sigma/r**2*h**(4)/omegak   # equation 7 from Ida 2020

# equations for tau_e, tau_a, and tau_m by Cresswell & Nelson 2008 as taken from Ida 2020
tau_e_CN = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)                                                    # equation 17 from Ida 2020
tau_m_CN = 1/((2.7+1.1*0.5)/(2)*h**2*(1-(ehat/2.02)**4)/(1+(ehat/2.25)**0.5+(ehat/2.84)**6)/t_wave)     # equation 18 from Ida 2020
tau_a_CN = 1/(2/tau_m_CN+2*e**2/tau_e_CN)                                                               # equation 16 from Ida 2020    

# p = -dln(sigma)/dln(r), q = -dln(T)/dln(r)
p, q = 1.0, 0.5
# equations for Cm and Ct taken from table 1 in Ida 2020
Ct = 2.73+1.08*p+0.87*q
Cm = 6*(2*p-q+2)

# equations for tau_e, tau_a, tau_m and tau_i by Ida 2020
tau_e = 1/(0.780/t_wave*(1+1/15*(ehat**2+ihat**2)**(3/2))**(-1))            # equation D1 in Ida 2020
tau_i = 1/(0.544/t_wave*(1+1/21.5*(ehat**2+ihat**2)**(3/2))**(-1))          # equation D2 in Ida 2020
tau_a = 1/(Ct*h**2*(1+Ct/Cm*(ehat**2+ihat**2)**(1/2))**(-1)/t_wave)         # equation D4 in Ida 2020
tau_m = 1/(0.5/tau_a-e**2/tau_e-i**2/tau_i)                                 # equation D3 in Ida 2020

# %%
fig, ax = plt.subplots(2, figsize=(7,7))

x = ehat

ax[0].plot(x, tau_e/year, label='IDA20')
ax[0].plot(x, tau_e_CN/year, label='CN08')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel(r'$\tau_e$ (years)')
# ax[0].set_ylabel(r'$\tau^{-1}_e/t_{wave}^{-1}$')
ax[0].set_xlim(0.1, 10)
ax[0].tick_params(which='both', direction="in", top=True, right=True)
ax[0].legend()

ax[1].plot(x, tau_a/year, label='IDA20')
ax[1].plot(x, tau_a_CN/year, label='CN08')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel('e/h')
ax[1].set_ylabel(r'$\tau_a$ (years)')
# ax[0].set_ylabel(r'$\tau^{-1}_a/t_{wave}^{-1}$')
ax[1].set_xlim(0.1, 10)
ax[1].tick_params(which='both', direction="in", top=True, right=True)
ax[1].legend()

# fig.savefig('/home/john/Desktop/summerproject/img/eccentricity_timescales.png', bbox_inches='tight')