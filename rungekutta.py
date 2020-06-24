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

# G = 4*np.pi**2
# mstar = 1
# mplanet = 3e-6
# dt = 0.001
# a = 1
# sigma = 17000*a**(-3/2)*au**2/m
# lim = 2

G = 6.674e-11
mstar = 1.989e30
mplanet = 5.972e24
dt = year*0.001
a = au
sigma = 17000*(a/au)**(-3/2)
lim = au*1.5

noutputs = 10000000
totaltime = noutputs*dt/year

h = 0.05

e = 0.3
ehat = e/h
r_p = a*(1-e)
r_a = a*(1+e)

omegak = np.sqrt(G*(mstar)/a**3)

t_wave = mstar/mplanet*mstar/sigma/a**2*h**4/omegak

p, q = 1.0, 0.5
Ct = 2.73+1.08*p+0.87*q
Cm = 6*(2*p-q+2)

tau_e = t_wave/0.780*(1+1/15*ehat**3)
tau_a = t_wave/Ct/h**2*(1+Ct/Cm*ehat)
tau_m = (0.5/tau_a-e**2/tau_e)**-1

vorb = np.sqrt(G*mstar*(2/r_a-1/a))

W0 = np.array([0,r_a,-vorb,0])
# %%
def acceleration(W0):
    x = W0[0:2]
    v = W0[2:5]
    r = np.linalg.norm(x)
    vk = np.linalg.norm(v)
    
    dvdtG = -G*mstar*x/r**3
    
    uv_r = x/r # unit vector in radial direction
    uv_a = np.array([-uv_r[1], uv_r[0]]) # unit vector in azimuthal direction
    vr = np.dot(uv_r, v) # radial velocity
    vtheta = np.dot(uv_a, v) # azimuthal velocity
    
    dvdt1 = -2*(np.dot(v,x)*x)/r**2/tau_e # equation 15 in Creswell+Nelson 2008
    # print('1: ', dvdt1)
    # dvdt1 = -vk/2/tau_a*uv_a-vr/tau_e*uv_r-(vtheta-vk)/tau_e*uv_a # equation 46 in Ida 2020
    # print('3: ', dvdt2)
    
    # print(dvdt1)
    
    return np.hstack((v, dvdtG+dvdt1))

def rungekutta(W0):
    W = np.zeros((noutputs,4))
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

df = pd.DataFrame(W)
df.to_csv('/home/john/Desktop/summerproject/data/10thousandyears.csv')
# %%
times = np.linspace(0,noutputs*dt,noutputs) # for displaying time elapsed in plot
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
    text.set_text('{:.1f} years'.format(times[i]/year*1000))
    return star, planet, planetline, text
    
anim = animation.FuncAnimation(fig, animate, frames=noutputs, interval=1, blit=True)
# %%
Wf = pd.read_csv('/home/john/Desktop/summerproject/data/10thousandyears.csv').values[:,1:5]

W = Wf[0::1000]

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

# plt.savefig('/home/john/Desktop/summerproject/img/10thousandyears.png', bbox_inches='tight')