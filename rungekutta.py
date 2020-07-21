#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:55:33 2020

@author: johngillan
"""

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from timeit import default_timer as timed


year = 365.25*24.*60.*60.
au = 1.496e11
m = 1.989e30

G = 4*np.pi**2
mstar = 1
mplanet = 3e-5
dt = 0.01
r = 1
sigma = 17000*r**(-1/2)*au**2/m
lim = 1.5

# G = 6.674e-11
# mstar = 1.989e30
# mplanet = 5.972e25
# dt = year*0.02
# r = au*5
# sigma = 17000*(r/au)**(-3/2)
# lim = au*1.5

mu = G*mstar
noutputs = 1000
totaltime = noutputs*dt
h = 0.05
i = 0
ihat = i/h

e = 0.15
ehat = e/h
r_p = r*(1-e) # perihelion
r_a = r*(1+e) # aphelion

Cm = 21.
Ct = 4.25

vorb0 = np.sqrt(mu*(2/r_a-1/r)) # instantaneous orbital speed at aphelion

W0 = np.array([0,r_a,0,-vorb0,0,0])
# %%
def acceleration(W0):
    R = W0[0:3]                             # position vector
    V = W0[3:6]                             # velocity vector
    r = np.linalg.norm(R)                   # magnitude of position vector
    v = np.linalg.norm(V)                   # magnitude of velocity vector
    vk = np.sqrt(mu/r)                      # instantaneous keplerian velocity
    
    u_r = R/r                               # unit vector in radial direction
    vr = np.dot(V, u_r)                     # radial velocity
    Va = V-(vr*u_r)                         # azimuthal velocity vector
    va = np.linalg.norm(Va)                 # azimuthal velocity
    u_a = Va/va                             # unit vector in azimuthal direction
    
    dvdtG = -mu*R/r**3                      # acceleration due to gravity
    
    e = np.linalg.norm((v**2/mu - 1/r)*R - np.dot(R,V)/mu*V) # eccentricity
    ehat = e/h                   
    
    omegak = np.sqrt(mu/r**3)               # new omegak every timestep
    t_wave = mstar/mplanet*mstar/sigma/r**2*h**4/omegak # new t_wave every timestep
    
    # Ida prescription
    # tau_e = t_wave/0.780*(1+1/15*(ehat**2+ihat**2)**(3/2))
    # tau_i = t_wave/0.544*(1+1/21.5*(ehat**2+ihat**2)**(3/2))
    # tau_a = t_wave/Ct/h**2*(1+Ct/Cm*(ehat**2+ihat**2)**(1/2))
    # tau_m = 1/(0.5/tau_a-e**2/tau_e-i**2/tau_i)
    
    # CN prescription
    tau_e = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)
    tau_m = 1/((2.7+1.1*0.5)/2*h**2*(1-(ehat/2.02)**4)/(1+(ehat/2.25)**0.5+(ehat/2.84)**6)/t_wave)
    tau_a = 1/(2/tau_m+2*e**2/tau_e)
    
    # dvdt2 = -vk/2/tau_a*u_a-vr/tau_e*u_r-(va-vk)/tau_e*u_a # equation 46 in Ida 2020
    dvdt2 = -V/tau_m-2*np.dot(V,R)*R/r**2/tau_e # equation 15 in Creswell+Nelson 2008
    
    return np.hstack((V, dvdtG+dvdt2))

def rungekutta(W0):
    
    W = np.zeros((noutputs,6))
    
    for i in range(noutputs):
        fa = acceleration(W0)
        Wb = W0 + dt/2*fa
        fb = acceleration(Wb)
        Wc = W0 + dt/2*fb
        fc = acceleration(Wc)
        Wd = W0 + dt*fc
        fd = acceleration(Wd)    
        W0 = W0 + dt/6*fa + dt/3*fb + dt/3*fc + dt/6*fd
        
        W[i] = W0
        
    return W

timer = timed()
W = rungekutta(W0)
print(timed()-timer)

omegak = np.sqrt(mu/r**3)
t_wave = mstar/mplanet*mstar/sigma/r**2*h**4/omegak

R = W[:,0:3]
V = W[:,3:6]
r = np.linalg.norm(R, axis=1)
v = np.linalg.norm(V, axis=1)

# e_total_low = np.linalg.norm([ (v[i]**2/mu - 1/r[i])*R[i] - (np.dot(R[i],V[i])*V[i])/mu for i in range(len(W)) ], axis=1)
# e_total_mid = np.linalg.norm([ (v[i]**2/mu - 1/r[i])*R[i] - (np.dot(R[i],V[i])*V[i])/mu for i in range(len(W)) ], axis=1)
e_total_high = np.linalg.norm([ (v[i]**2/mu - 1/r[i])*R[i] - (np.dot(R[i],V[i])*V[i])/mu for i in range(len(W)) ], axis=1)

E_total = v**2/2-mu/r

a_total = -mu/2/E_total

# e_results_low = np.zeros(noutputs)
# e_results_mid = np.zeros(noutputs)
e_results_high = np.zeros(noutputs)

a_results = np.zeros(noutputs)

for j in range(noutputs):    
    ehat = e/h
    
    # tau_e = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    tau_e = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
    
    de = -e/tau_e*dt
    e = e+de
    # e_results_low[j] = e
    # e_results_mid[j] = e
    e_results_high[j] = e
    
    # tau_m = 1/((2.7+1.1*0.5)/(2)*h**2*(1-(ehat/2.02)**4)/(1+(ehat/2.25)**0.5+(ehat/2.84)**6)/t_wave)
    # tau_a  = 1/(2/tau_m+2*e**2/tau_e)
    
    # da = -a/t_a*dt
    # a = a+da
    # a_results[j] = a
# %%
fig, ax = plt.subplots(1, figsize=(7,5))

# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_results_low, label='analytical', c='tab:blue')
# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_total_low, linestyle='--', label='numerical', c='tab:orange')
# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_results_mid, c='tab:blue')
# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_total_mid, linestyle='--', c='tab:orange')
ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_results_high, c='tab:blue')
ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_total_high, linestyle='--', c='tab:orange')
ax.axhline(h, c='black', label='e = H/r')
ax.set_xlabel('time (years)')
ax.set_ylabel('e')
ax.set_xlim(0, noutputs*dt/2/np.pi*omegak)
ax.set_ylim(0)
ax.tick_params(which='both', direction="in", top=True, right=True)
ax.grid()
ax.legend()

fig.savefig('/home/john/Desktop/summerproject/img/num+analyticIDA.pdf', bbox_inches='tight')
# %%
times = np.linspace(0,totaltime,noutputs) # for displaying time elapsed in plot
fig, ax = plt.subplots(1, figsize=(9, 9))

ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
ax.set_xlabel('Distance (AU)')
ax.set_ylabel('Distance (AU)')
ax.grid()

star, = ax.plot([], [], 'o')
planet, = ax.plot([], [], 'o', c='tab:blue')
planetline, = ax.plot([], [], lw=0.1, c='tab:blue')
text = ax.text(-lim+lim/4, +lim-lim/4, s='', fontsize=15)

def animate(i):
    star.set_data(0,0)
    planet.set_data(W[i,0],W[i,1])
    planetline.set_data(W[0:i,0],W[0:i,1])
    text.set_text('{:.1f} years'.format(times[i]))
    return star, planet, planetline, text
    
anim = animation.FuncAnimation(fig, animate, frames=noutputs, interval=1, blit=True)
# %%
fig, ax = plt.subplots(1, figsize=(8, 8))
ax.set_title('{:.0f} years'.format(totaltime))
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
ax.set_xlabel('Distance (AU)')
ax.set_ylabel('Distance (AU)')
ax.grid()
ax.scatter(0,0, c='black')
ax.scatter(W[-1,0], W[-1,1], c='steelblue')
ax.plot(W[:,0], W[:,1], linewidth=0.01, c='steelblue')

# plt.savefig('/home/john/Desktop/summerproject/img/5000years.png', bbox_inches='tight')
