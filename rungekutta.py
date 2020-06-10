#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:55:33 2020

@author: john
"""

# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from timeit import default_timer as timed

G = 6.674e-11                               # gravitational constanst
au = 1.496e11                               # astronomical unit
mstar = 2.0**30                             # mass of star
mplanet = 1.0**24                           # mass of orbiting planet
r = 0.1*au                                  # radial position of planet   
omegak = np.sqrt(G*(mstar+mplanet)/r**3)    # Keplerian frequency      
vorb = np.sqrt(G*(mstar+mplanet)/r)
noutputs = 100000
h = 1e14

W = np.array([[0,0,0,0,0,0],
      [0,r,0,-vorb*1.1,0,0]])

# %%
# force equation
def force(W):        
    sum_f_star = - G * mplanet * (W[0,0:3] - W[1,0:3]) / (np.linalg.norm(W[0,0:3] - W[1,0:3]) ** 3)
    sum_f_planet = - G * mstar * (W[1,0:3] - W[0,0:3]) / (np.linalg.norm(W[1,0:3] - W[0,0:3]) ** 3)
    g_star = np.hstack((W[0,3:6], sum_f_star))
    g_planet = np.hstack((W[1,3:6], sum_f_planet))
    g = np.vstack((g_star,g_planet))
    return g

# runge kutta integrator
def rungekutta(W):
    W_total = np.zeros((noutputs, 2, 6))
    for i in range(noutputs): 
        fa = force(W)
        Wb = W + h/2*fa
        fb = force(Wb)
        Wc = W + h/2*fb
        fc = force(Wc)
        Wd = W + h*fc
        fd = force(Wd)    
        W = W + h/6*fa + h/3*fb + h/3*fc + h/6*fd
        W_total[i] = W
    return W_total    

timer = timed()
W_total = rungekutta(W) # 6 dimensional phase space for each particle at each timestep
print(timed()-timer)

# %%
lim = 0.15
fig, axes = plt.subplots(1, figsize=(8, 8))
axes.set_xlim(-lim, lim)
axes.set_ylim(-lim, lim)
axes.set_aspect('equal')
axes.set_xlabel('Distance (AU)')
axes.set_ylabel('Distance (AU)')
axes.scatter(W_total[-1,0,0]/au, W_total[-1,0,1]/au, color='gold', label='Star')
axes.scatter(W_total[-1,1,0]/au, W_total[-1,1,1]/au, label='Planet')
axes.plot(W_total[:,1,0]/au, W_total[:,1,1]/au, linewidth=0.2)
axes.legend()
# plt.savefig('solar0001.png', bbox_inches='tight')

# %%
fig, axes = plt.subplots(1, figsize=(9, 9))

axes.set_xlim(-lim, lim)
axes.set_ylim(-lim, lim)
axes.set_aspect('equal')
axes.set_xlabel('Distance (AU)')
axes.set_ylabel('Distance (AU)')

star, = axes.plot([], [], color='gold', marker='o', label='Star')
planet, = axes.plot([], [], 'o', color='blue')
planetline, = axes.plot([], [], lw=0.2, color='blue')

def animate(i):
    star.set_data(W_total[i,0,0]/au, W_total[i,0,1]/au)
    planet.set_data(W_total[i,1,0]/au, W_total[i,1,1]/au)
    planetline.set_data(W_total[0:i,1,0]/au, W_total[0:i,1,1]/au)
    return star, planet, planetline

    
anim = FuncAnimation(fig, animate, frames=noutputs, interval=1, blit=True)
# anim.save('solar.gif')

# %%