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


year = 365.25*24.*60.*60.
au = 1.496e11
m = 1.99e30

# G = 4*np.pi**2
# mstar = 1
# mplanet = 3e-6
# dt = 0.001
# a = 1
# lim = 2

G = 6.674e-11
mstar = 1.99e30
mplanet = 6.0e24
dt = year*0.001
a = au
lim = au*2

noutputs = 3000000
totaltime = noutputs*dt/year

h = 0.05

# initial position and velocity of star and planet [x,y,z,vx,vy,vz]
W0 = np.array([0,r_a,0,-vorb,0,0])

e = 0.3
ehat = e/h
r_p = a*(1-e)
r_a = a*(1+e)
sigma = 17000*(a/au)**(-3/2)
omegak = np.sqrt(G*(mstar)/a**3)

t_wave = mstar/mplanet*mstar/sigma/a**2*h**4/omegak

p, q = 1.0, 0.5
Ct = 2.73+1.08*p+0.87*q
Cm = 6*(2*p-q+2)

tau_e = t_wave/0.780*(1+1/15*ehat**3)
tau_a = t_wave/Ct/h**2*(1+Ct/Cm*ehat)
tau_m = (0.5/tau_a-e**2/tau_e)**-1

vorb = np.sqrt(G*mstar*(2/r_a-1/a))



# %%
# force equation
def force(W0):
    x = W0[0:3]
    v = W0[3:6]
    r = np.linalg.norm(x)
    vmag = np.linalg.norm(v)
    
    dvdt1 = -G*mstar*x/r**3
    
    u_norm = x/r
    u_tang = np.array([-u_norm[1], u_norm[0], 0])
    v_norm = np.dot(u_norm, v)
    v_tang = np.dot(u_tang, v)
    
    dvdt2 = -2*(np.dot(v,x)*x)/r**2/tau_e
    # print(dvdt2)
    dvdt2 = -v/tau_m-2*v_norm/tau_e*u_norm
    # print(dvdt2)
    dvdt2 = -vmag/2/tau_a*u_tang-v_norm/tau_e*u_norm-(v_tang-vmag)/tau_e*u_tang    
    # print(dvdt2)
    
    W = np.hstack((v, dvdt1+dvdt2))
    print(dvdt1)
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
fig, ax = plt.subplots(1, figsize=(9, 9))

ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect('equal')
ax.set_xlabel('Distance (AU)')
ax.set_ylabel('Distance (AU)')
ax.grid()

star, = ax.plot([], [], 'o')
planet, = ax.plot([], [], 'o')
planetline, = ax.plot([], [], lw=0.5, color='orange')
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