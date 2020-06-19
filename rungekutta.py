#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:55:33 2020

@author: john
"""

# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from timeit import default_timer as timed

G = 4*np.pi**2                              # gravitational constsant
au = 1.496e11
m = 1.99e30
mstar = 1                                   # mass of star
mplanet = 5.97e-6                           # mass of orbiting planet
a = 1

noutputs = 100000                           # number of timesteps
dt = 0.01                                    # timestep
totaltime = noutputs*dt
h = 0.05
e = 0.3
ehat = e/h
r_p = a*(1-e)
r_a = a*(1+e)
sigma = 17000*a**(-3/2)*au**2/m
omegak = np.sqrt(G*(mstar)/a**3)                # Keplerian frequency
t_wave = mstar/mplanet*mstar/sigma/a**2*h**(4)/omegak
tau_e = t_wave/0.780*(1+1/15*(ehat**3))

vorb = np.sqrt(G*mstar*(2/r_a-1/a))

# initial position and velocity of star and planet [x,y,z,vx,vy,vz]
W0 = np.array([[0,r_a,0,-vorb,0,0]])

# %%
# force equation
def force(W0):
    r = np.linalg.norm(W0[0,0:3])
    acc_g = -G*mstar*(W0[0,0:3])/r**3
    acc_e = -2*(np.dot(W0[0,3:6],W0[0,0:3])*W0[0,0:3])/r**2/tau_e
    # print(acc_e)
    W = np.hstack((W0[0,3:6], acc_g))
    return W

# runge kutta integrator
def rungekutta(W0):
    W = np.zeros((noutputs, 1, 6))
    for i in range(noutputs): 
        fa = force(W0)
        Wb = W0 + dt/2*fa
        fb = force(Wb)
        Wc = W0 + dt/2*fb
        fc = force(Wc)
        Wd = W0 + dt*fc
        fd = force(Wd)    
        W0 = W0 + dt/6*fa + dt/3*fb + dt/3*fc + dt/6*fd
        W[i] = W0
    return W    

timer = timed()
W = rungekutta(W0) # 6 dimensional phase space for each particle at each timestep
print(timed()-timer)
# %%
times = np.linspace(0,noutputs*dt,noutputs) # for displaying time elapsed in plot
lim = 2
fig, ax = plt.subplots(1, figsize=(9, 9))

ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
ax.set_xlabel('Distance (AU)')
ax.set_ylabel('Distance (AU)')
ax.grid()

star, = ax.plot([], [], 'o')
planet, = ax.plot([], [], 'o')
planetline, = ax.plot([], [], lw=0.1, color='orange')
text = ax.text(-lim+0.1, lim-0.2, s='', fontsize=15)

def animate(i):
    star.set_data(0,0)
    planet.set_data(W[i,0,0],W[i,0,1])
    planetline.set_data(W[0:i,0,0],W[0:i,0,1])
    text.set_text('{:.1f} years'.format(times[i]))
    return star, planet, planetline, text
    
anim = animation.FuncAnimation(fig, animate, frames=noutputs, interval=1, blit=True)
# anim.save('solar.gif')
# %%
lim = 2
fig, ax = plt.subplots(1, figsize=(8, 8))
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
ax.set_xlabel('Distance (AU)')
ax.set_ylabel('Distance (AU)')
ax.grid()
ax.scatter(0,0)
ax.scatter(W[-1,0,0], W[-1,0,1], c='orange')
ax.plot(W[:,0,0], W[:,0,1], linewidth=0.1, c='orange')
# plt.savefig('solar0001.png', bbox_inches='tight')
# %%