#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 09:41:39 2020

@author: john
"""

# %%
import rebound
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from scipy.spatial.distance import pdist
from timeit import default_timer as timed


# constants
G = 6.67428e-11                     # gravitational constanct in SI units
au = 1.496e11                       # astronomical unit    
Msun = 1.9891e30                    # mass of sun
year = 365.25*24.*60.*60.           # number of seconds in a year
m1 = 5.972e25
rsun = 1.*au                       # distance of centre of mass of binary from the sun 
mu = G*Msun
OmegaK = np.sqrt(mu/rsun**3) # keplerian frequency at this distance

h = 0.05

e = 0.25
ehat = e/h
r_p = rsun*(1-e) # perihelion
r_a = rsun*(1+e) # aphelion

Cm = 21.
Ct = 4.25

vorb0 = np.sqrt(mu*(2/r_a-1/rsun)) # instantaneous orbital speed at apoapsis

sim = rebound.Simulation()
sim.G = G
sim.dt = 0.01*year 
sim.integrator = "ias15"
sim.gravity    = "basic"
sim.collision  = "none"
sim.ri_ias15.epsilon=0

sim.add(m=Msun, hash="sun")
sim.add(m=m1, y=r_a, vx=-vorb0, hash="primary")

Noutputs = 200000
totaltime = Noutputs*dt

times = np.linspace(0.,totaltime, Noutputs)

p, sun = np.zeros((Noutputs, 3)), np.zeros((Noutputs, 3))
vp, vsun = np.zeros((Noutputs, 3)), np.zeros((Noutputs, 3))

ps = sim.particles

# Cd = 2.
# rho_g = 1e-20
# drag1 = 0.5*Cd*np.pi*s1**2*rho_g
# drag2 = 0.5*Cd*np.pi*s2**2*rho_g
# drag3 = 0.5*Cd*np.pi*simp**2*rho_g

# def dragForce(reb_sim):
#     ps["primary"].ax -= (ps["primary"].vx)**2*drag1
#     ps["primary"].ay -= (ps["primary"].vy)**2*drag1
#     ps["primary"].az -= (ps["primary"].vz)**2*drag1
#     ps["secondary"].ax -= (ps["secondary"].vx)**2*drag2
#     ps["secondary"].ay -= (ps["secondary"].vy)**2*drag2
#     ps["secondary"].az -= (ps["secondary"].vz)**2*drag2
#     ps["impactor"].ax -= (ps["impactor"].vx)**2*drag3
#     ps["impactor"].ay -= (ps["impactor"].vy)**2*drag3
#     ps["impactor"].az -= (ps["impactor"].vz)**2*drag3
    
# sim.additional_forces = dragForce
# sim.force_is_velocity_dependent = 1


timer = timed()
for i, time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    p[i] = [ps["primary"].x, ps["primary"].y, ps["primary"].z]
    sun[i] = [ps["sun"].x, ps["sun"].y, ps["sun"].z]
print(timed()-timer)
# %%
lim = 2*au
fig, axes = plt.subplots(1, figsize=(9, 9))
axes.set_ylim(-lim,lim)
axes.set_xlim(-lim,lim)
primaryline, = axes.plot([], [], label="primary", c="tab:orange", lw=1.2)
primarydot, = axes.plot([], [], marker="o", ms=7, c="tab:orange")
axes.grid()
axes.legend()

def animate(i):
    # primaryline.set_data(p[0:i], p[0:i])
    primarydot.set_data(p[i], p[i])
    text.set_text('{} Years'.format(np.round(times[i]/(year), 1)))
    return primarydot, primaryline

anim = animation.FuncAnimation(fig, animate, frames=Noutputs, interval=1,blit=True)
# anim.save(f'{path}/videos/2D.mp4')