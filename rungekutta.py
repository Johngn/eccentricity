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
a = 1
sigma = 17000*a**(-3/2)*au**2/m
lim = 1.5

# G = 6.674e-11
# mstar = 1.989e30
# mplanet = 5.972e25
# dt = year*0.01
# a = au*5
# sigma = 17000*(a/au)**(-3/2)
# lim = au*1.5

mu = G*mstar

noutputs = 150000
totaltime = noutputs*dt

h = 0.05

e = 0.25
ehat = e/h
r_p = a*(1-e) # perihelion
r_a = a*(1+e) # aphelion

Cm = 21.
Ct = 4.25

vorb0 = np.sqrt(mu*(2/r_a-1/a)) # instantaneous orbital speed at aphelion

W0 = np.array([0,r_a,0, -vorb0,0,0])

def acceleration(W0):
    x = W0[0:3]
    v = W0[3:6]
    r = np.linalg.norm(x)
    vk = np.linalg.norm(v)
    
    u_r = x/r                               # unit vector in radial direction 
    vr = np.dot(v, u_r)                     # radial velocity
    vtheta = np.linalg.norm(v-(vr*u_r))     # azimuthal velocity
    u_a = (v-(vr*u_r))/vtheta               # unit vector in azimuthal direction
    
    dvdtG = -mu*x/r**3      # acceleration due to gravity
    
    e = np.linalg.norm( (vk**2/mu - 1/r)*x - (np.dot(x,v)*v)/mu ) # eccentricity vector
    ehat = e/h
    E = vk**2/2-mu/r
    a = -mu/2/E
    
    omegak = np.sqrt(mu/r**3)

    t_wave = mstar/mplanet*mstar/sigma/a**2*h**4/omegak    
    
    # Ida prescription
    # tau_e = t_wave/0.780*(1+1/15*ehat**3)
    # tau_a = t_wave/Ct/h**2*(1+Ct/Cm*ehat)
    # tau_m = (0.5/tau_a-e**2/tau_e)**-1
    
    # CN prescription
    tau_e = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)
    tau_m = 1/((2.7+1.1*0.5)/(2)*h**2*(1-(ehat/2.02)**4)/(1+(ehat/2.25)**0.5+(ehat/2.84)**6)/t_wave)
    tau_a = 1/(2/tau_m+2*e**2/tau_e)
    
    # dvdt2 = -vk/2/tau_a*u_a-vr/tau_e*u_r-(vtheta-vk)/tau_e*u_a # equation 46 in Ida 2020
    dvdt2 = -v/tau_m-2*(np.dot(v,x)*x)/r**2/tau_e # equation 15 in Creswell+Nelson 2008
    
    return np.hstack((v, dvdtG+dvdt2))

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

omegak = np.sqrt(mu/a**3)
t_wave = mstar/mplanet*mstar/sigma/a**2*h**4/omegak

x = W[:,0:3]
v = W[:,3:6]
r = np.linalg.norm(x, axis=1)
vk = np.linalg.norm(v, axis=1)

e_total = np.linalg.norm([ (vk[i]**2/mu - 1/r[i])*x[i] - (np.dot(x[i],v[i])*v[i])/mu for i in range(len(W)) ], axis=1)

E_total = vk**2/2-mu/r

a_total = -mu/2/E_total

e_results = np.zeros(noutputs)
a_results = np.zeros(noutputs)

for i in range(noutputs):    
    ehat = e/h
    # tau_e = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    tau_e = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
    # tau_m = 1/((2.7+1.1*0.5)/(2)*h**2*(1-(ehat/2.02)**4)/(1+(ehat/2.25)**0.5+(ehat/2.84)**6)/t_wave)
    # tau_a  = 1/(2/tau_m+2*e**2/tau_e)   
    de = -e/tau_e*dt
    e = e+de
    e_results[i] = e
    # da = -a/t_a*dt
    # a = a+da
    # a_results[i] = a
# %%
fig, ax = plt.subplots(1, figsize=(7,5))

ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_results, label='analytical')
ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_total, label='numerical')
# ax.plot(np.arange(0,noutputs,1)*dt/year, e_results2)
# ax.plot(np.arange(0,noutputs,1)*dt/year, e_results3)
# ax.plot(np.arange(0,noutputs,1)*dt/year, e_results4)
ax.axhline(h, c='black', label='e = H/r')
ax.set_xlabel('time (years)')
ax.set_ylabel('e')
ax.set_xlim(0, noutputs*dt/2/np.pi*omegak)
ax.set_ylim(0)
ax.tick_params(which='both', direction="in", top=True, right=True)
ax.grid()
ax.legend()

fig.savefig('/home/john/Desktop/summerproject/img/high_e_num+analyticCN.png', bbox_inches='tight')
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
planet, = ax.plot([], [], 'o', c='steelblue')
planetline, = ax.plot([], [], lw=0.1, c='steelblue')
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
ax.plot(W[:,0], W[:,1], linewidth=0.005, c='steelblue')

# plt.savefig('/home/john/Desktop/summerproject/img/5000years.png', bbox_inches='tight')
