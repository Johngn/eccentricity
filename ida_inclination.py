#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:51:19 2020

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
r = 10*au                                  # radial position of planet
omegak = np.sqrt(G*(mstar+mplanet)/r**3)    # Keplerian frequency
h = 0.05                                    # scale height of disc
sigma = 17000*(r/au)**(-1)*mstar          # surface density of disc
i = np.arange(0, h, 0.0001)                 # inclinations
ihat = i/h                                  # inclinations divided by scale height of disc
e = [0,0.01,0.02,0.03,0.04,0.05]
colors = ['red','orange','yellow','green','blue','black']


fig, ax = plt.subplots(4, figsize=(6,12))

for j in range(len(e)):    
    ehat = e[j]/h    
    
    # inverse of equation 7 from Ida 2020
    t_wave = 1/((mplanet/mstar)*(sigma*r**2/mstar)*h**(-4)*omegak)    
    
    # equation from appendix of Ida 2020 for inclination damping
    tau_e = 1/(0.780/t_wave*(1+1/15*(ehat**2+ihat**2)**(3/2))**(-1)) # equation D1 in Ida 2020
    tau_i = 1/(0.544/t_wave*(1+1/21.5*(ehat**2+ihat**2)**(3/2))**(-1)) # equation D2 in Ida 2020    
    
    # p = -dln(sigma)/dln(r), q = -dln(T)/dln(r)
    p, q = 1.0, 0.5
    # equations for Cm and Ct taken from table 1 in Ida 2020
    Ct = 2.73+1.08*p+0.87*q
    Cm = 6*(2*p-q+2)
    
    tau_a = 1/(Ct*h**2*(1+Ct/Cm*(ehat**2+ihat**2)**(1/2))**(-1)/t_wave) # equation D4 in Ida 2020    
    tau_m = 1/(0.5/tau_a-e[j]**2/tau_e-i**2/tau_i)
    
    
    # tau_e as a function of i/h
    ax[0].plot(ihat, 1/tau_e*t_wave, label=f'e = {e[j]}', c=colors[j])
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$\tau^{-1}_e/t_{wave}^{-1}$')
    ax[0].set_xlim(0.01, 1)
    # ax[0].set_ylim(1e-3, 1e1)
    ax[0].tick_params(which='both', direction="in", top=True, right=True)
    ax[0].legend()
    
    # tau_m as a function of i/h
    ax[1].plot(ihat, 1/tau_m*t_wave, label=f'e = {e[j]}', c=colors[j])
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$\tau^{-1}_m/t_{wave}^{-1}$')
    ax[1].set_xlim(0.01, 1)
    # ax[1].set_ylim(1e-4, 1e-1)
    ax[1].tick_params(which='both', direction="in", top=True, right=True)
    ax[1].legend()
    
    # tau_a as a function of i/h
    ax[2].plot(ihat, 1/tau_a*t_wave, label=f'e = {e[j]}', c=colors[j])
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].set_ylabel(r'$\tau^{-1}_a/t_{wave}^{-1}$')
    ax[2].set_xlim(0.01, 1)
    # ax[2].set_ylim(1e-4, 1e-1)
    ax[2].tick_params(which='both', direction="in", top=True, right=True)
    ax[2].legend()
    
    # tau_i as a function of i/h
    ax[3].plot(ihat, 1/tau_i*t_wave, label=f'e = {e[j]}', c=colors[j])
    ax[3].set_xscale('log')
    ax[3].set_yscale('log')
    ax[3].set_xlabel('i/h')
    ax[3].set_ylabel(r'$\tau^{-1}_i/t_{wave}^{-1}$')
    ax[3].set_xlim(0.01, 1)
    # ax[3].set_ylim(1e-4, 1e-1)
    ax[3].tick_params(which='both', direction="in", top=True, right=True)
    ax[3].legend()

fig.savefig('/home/john/Desktop/summerproject/img/inclination_plot.png', bbox_inches='tight')

# %%

fig, ax = plt.subplots(4, figsize=(6,12))

for j in range(len(e)):    
    ehat = e[j]/h    
    
    # inverse of equation 7 from Ida 2020
    t_wave = 1/((mplanet/mstar)*(sigma*r**2/mstar)*h**(-4)*omegak)
    
    # equation from appendix of Ida 2020 for inclination damping
    tau_e = 1/(0.780/t_wave*(1+1/15*(ehat**2+ihat**2)**(3/2))**(-1)) # equation D1 in Ida 2020
    tau_i = 1/(0.544/t_wave*(1+1/21.5*(ehat**2+ihat**2)**(3/2))**(-1)) # equation D2 in Ida 2020    
    
    # p = -dln(sigma)/dln(r), q = -dln(T)/dln(r)
    p, q = 1.0, 0.5
    # equations for Cm and Ct taken from table 1 in Ida 2020
    Ct = 2.73+1.08*p+0.87*q
    Cm = 6*(2*p-q+2)
    
    tau_a = 1/(Ct*h**2*(1+Ct/Cm*(ehat**2+ihat**2)**(1/2))**(-1)/t_wave) # equation D4 in Ida 2020    
    tau_m = 1/(0.5/tau_a-e[j]**2/tau_e-i**2/tau_i)
    

    # tau_e as a function of i/h
    ax[0].plot(ihat, tau_e*year, label=f'e = {e[j]}', c=colors[j])
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$\tau_e$ (yr)')
    ax[0].set_xlim(0.01, 1)
    # ax[0].set_ylim(1e-2, 1e3)
    ax[0].tick_params(which='both', direction="in", top=True, right=True)
    ax[0].legend()
    
    # tau_m as a function of i/h
    ax[1].plot(ihat, tau_m*year, label=f'e = {e[j]}', c=colors[j])
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$\tau_m$ (yr)')
    ax[1].set_xlim(0.01, 1)
    # ax[1].set_ylim(1e-4, 1e-1)
    ax[1].tick_params(which='both', direction="in", top=True, right=True)
    ax[1].legend()
    
    # tau_a as a function of i/h
    ax[2].plot(ihat, tau_a*year, label=f'e = {e[j]}', c=colors[j])
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].set_ylabel(r'$\tau_a$ (yr)')
    ax[2].set_xlim(0.01, 1)
    # ax[2].set_ylim(1e-4, 1e-1)
    ax[2].tick_params(which='both', direction="in", top=True, right=True)
    ax[2].legend()
    
    # tau_i as a function of i/h
    ax[3].plot(ihat, tau_i*year, label=f'e = {e[j]}', c=colors[j])
    ax[3].set_xscale('log')
    ax[3].set_yscale('log')
    ax[3].set_xlabel('i/h')
    ax[3].set_ylabel(r'$\tau_i$ (yr)')
    ax[3].set_xlim(0.01, 1)
    # ax[3].set_ylim(1e-4, 1e-1)
    ax[3].tick_params(which='both', direction="in", top=True, right=True)
    ax[3].legend()

fig.savefig('/home/john/Desktop/summerproject/img/inclination_timescale.png', bbox_inches='tight')