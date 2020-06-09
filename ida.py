#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:53:07 2020

@author: john
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.misc import derivative

font = {'size':12}

matplotlib.rc('font', **font)


# %%
# constants
G = 6.674e-11
au = 1.496e11
mstar = 2.0**30
mplanet = 1.0**24   
h = 0.05            # scale height of disc
r = 0.1*au          # radial position of planet

epsilon = 5.68e-3*(r/au)**(-2.168)*mstar    # surface density of disc (from Crida 2009)
omegak = np.sqrt(G*(mstar+mplanet)/r**3)    # Keplerian frequency

# eccentricites
e = np.arange(0, 1, 0.0001)
ehat = e/h

# equation 7 from Ida 2020
twave = (mplanet/mstar)*(epsilon*r**2/mstar)*h**(-4)*omegak


# equation 34 from Ida 2020
invtau_e = 0.780*(1+(1/15)*ehat**3)**(-1)*twave





p1, p2 = 1.0, 0.5
q1, q2 = 0.5, 1.0
ct1 = 2.73+1.08*p1+0.87*q1
cm1 = 6*(2*p1-q1+2)
ct2 = 2.73+1.08*p2+0.87*q2
cm2 = 6*(2*p2-q2+2)

tauared = ct1*h**2*(1+(ct1/cm1)*ehat)**(-1)*twave
tauamag = ct2*h**2*(1+(ct2/cm2)*ehat)**(-1)*twave




taumred = 0.5*tauared - e**2*invtau_e
taummag = 0.5*tauamag - e**2*invtau_e

# figure 1a from Ida 2020 - red line
fig, ax = plt.subplots(3, figsize=(6,10))
ax[0].plot(ehat, invtau_e, label='Ida 2020', c='red')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel(r'$\tau^{-1}_e$')
ax[0].set_xlim(0.1, 10)
ax[0].set_ylim(1e-3, 1e1)
ax[0].tick_params(which='both',direction="in", top=True, right=True)
ax[0].legend()

# figure 1b from Ida 2020 - red line
ax[1].plot(ehat, taumred, label='Ida 2020', c='red')
ax[1].plot(ehat, -taumred, c='red', ls='--')
ax[1].plot(ehat, taummag, label='Ida 2020', c='magenta')
ax[1].plot(ehat, -taummag, c='magenta', ls='--')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_ylabel(r'$\tau^{-1}_m$')
ax[1].set_xlim(0.1, 10)
ax[1].set_ylim(1e-4, 1e-1)
ax[1].tick_params(which='both',direction="in", top=True, right=True)
ax[1].legend()

# figure 1b from Ida 2020 - red line
ax[2].plot(ehat, tauared, label='Ida 2020', c='red')
ax[2].plot(ehat, tauamag, label='Ida 2020', c='magenta')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
ax[2].set_xlabel('e/h')
ax[2].set_ylabel(r'$\tau^{-1}_a$')
ax[2].set_xlim(0.1, 10)
ax[2].set_ylim(1e-4, 1e-1)
ax[2].tick_params(which='both',direction="in", top=True, right=True)
ax[2].legend()