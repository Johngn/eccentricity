#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:51:19 2020

@author: john
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:53:07 2020

@author: john
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font = {'size':12}
matplotlib.rc('font', **font)

G = 6.674e-11                               # gravitational constanst
au = 1.496e11                               # astronomical unit
mstar = 2.0**30                             # mass of star
mplanet = 1.0**24                           # mass of orbiting planet
r = 10*au                                  # radial position of planet
omegak = np.sqrt(G*(mstar+mplanet)/r**3)    # Keplerian frequency
h = 0.05                                    # scale height of disc
sigma = 17000*(r/au)**(-2./2.)*mstar          # surface density of disc
i = np.arange(0, h, 0.0001)                 # inclinations
ihat = i/h                                  # inclinations divided by scale height of disc
e = 0
ehat = e/h


# equation 7 from Ida 2020
t_wave_inverse = (mplanet/mstar)*(sigma*r**2/mstar)*h**(-4)*omegak


# equation from appendix of Ida 2020 for inclination damping
tau_e_inverse = 0.780*(1+1/15*(ehat**2+ihat**2)**(3/2))**(-1) # equation D1 in Ida 2020
tau_i_inverse = 0.544*(1+1/21.5*(ehat**2+ihat**2)**(3/2))**(-1) # equation D2 in Ida 2020


# p = -dln(sigma)/dln(r), q = -dln(T)/dln(r)
p_red, q_red = 1.0, 0.5
# equations for Cm and Ct taken from table 1 in Ida 2020
Ct_red = 2.73+1.08*p_red+0.87*q_red
Cm_red = 6*(2*p_red-q_red+2)

# magenta line in plot using different values of p and q
p_mag, q_mag = 0.5, 1.0
Ct_mag = 2.73+1.08*p_mag+0.87*q_mag
Cm_mag = 6*(2*p_mag-q_mag+2)

tau_a_inverse_red = Ct_red*h**2*(1+Ct_red/Cm_red*(ehat**2+ihat**2)**(1/2))**(-1) # equation D4 in Ida 2020
tau_a_inverse_mag = Ct_mag*h**2*(1+Ct_mag/Cm_mag*(ehat**2+ihat**2)**(1/2))**(-1) # equation D4 in Ida 2020

tau_m_inverse_red = 0.5*tau_a_inverse_red-e**2*tau_e_inverse-i**2*tau_i_inverse
tau_m_inverse_mag = 0.5*tau_a_inverse_mag-e**2*tau_e_inverse-i**2*tau_i_inverse

fig, ax = plt.subplots(4, figsize=(6,12))

# tau_e as a function of i/h
ax[0].plot(ihat, tau_e_inverse, label='IDA20', c='red')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel(r'$\tau^{-1}_e$')
ax[0].set_xlim(0.01, 1)
# ax[0].set_ylim(1e-3, 1e1)
ax[0].tick_params(which='both', direction="in", top=True, right=True)
ax[0].legend()

# tau_m as a function of i/h
ax[1].plot(ihat, tau_m_inverse_mag, label='IDA20 p=0.5 q=1', c='magenta')
ax[1].plot(ihat, tau_m_inverse_red, label='IDA20 p=1 q=0.5', c='red')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel('e/h')
ax[1].set_ylabel(r'$\tau^{-1}_m$')
ax[1].set_xlim(0.01, 1)
# ax[1].set_ylim(1e-4, 1e-1)
ax[1].tick_params(which='both', direction="in", top=True, right=True)
ax[1].legend()

# tau_a as a function of i/h
ax[2].plot(ihat, tau_a_inverse_mag, label='IDA20 p=0.5 q=1', c='magenta')
ax[2].plot(ihat, tau_a_inverse_red, label='IDA20 p=1 q=0.5', c='red')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
ax[2].set_xlabel('e/h')
ax[2].set_ylabel(r'$\tau^{-1}_a$')
ax[2].set_xlim(0.01, 1)
# ax[2].set_ylim(1e-4, 1e-1)
ax[2].tick_params(which='both', direction="in", top=True, right=True)
ax[2].legend()

# tau_i as a function of i/h
ax[3].plot(ihat, tau_i_inverse, label='IDA20', c='red')
ax[3].set_xscale('log')
ax[3].set_yscale('log')
ax[3].set_ylabel(r'$\tau^{-1}_i$')
ax[3].set_xlim(0.01, 1)
# ax[3].set_ylim(1e-4, 1e-1)
ax[3].tick_params(which='both', direction="in", top=True, right=True)
ax[3].legend()

# fig.savefig('/home/john/Desktop/summerproject/img/inclination_plot_1.png', bbox_inches='tight')
