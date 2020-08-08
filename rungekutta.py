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
from mpl_toolkits.mplot3d import axes3d
from timeit import default_timer as timed
from glob import glob
# %%
year = 365.25*24.*60.*60.
au = 1.496e11
mstarkg = 1.989e30
G = 4*np.pi**2
mstar = 1
mplanet = 1e28/1.989e33
dt = 0.02
r = 1
sigma = 17000*r**(-1/2)*au**2/mstarkg
lim = 1.2
mu = G*mstar
noutputs = 1000000
totaltime = noutputs*dt
h = 0.05
e = 0.25
i = 0.0
ehat = e/h
ihat = i/h
r_p = r*(1-e) # perihelion
r_a = r*(1+e) # aphelion
Cm = 21.
Ct = 4.25
y = r_a*np.cos(i)
z = r_a*np.sin(i)
vorb = np.sqrt(mu*(2/r_a-1/r)) # instantaneous orbital speed at aphelion
W0 = np.array([0,y,z,-vorb,0,0])
omegak = np.sqrt(mu/r**3)
t_wave = mstar/mplanet*mstar/sigma/r**2*h**4/omegak

ida = True # select prescription
# %%
def acceleration(W0):
    R = W0[0:3]                             # position vector
    V = W0[3:6]                             # velocity vector
    r = np.linalg.norm(R)                   # magnitude of position vector
    v = np.linalg.norm(V)                   # magnitude of velocity vector    
    dvdtG = -mu*R/r**3                      # acceleration due to gravity
    vk = np.sqrt(mu/r)                      # instantaneous keplerian velocity
    omegak = np.sqrt(mu/r**3)               
    t_wave = mstar/mplanet*mstar/sigma/r**2*h**4/omegak
    
    theta = np.arctan2(R[1],R[0])                           # angle of rotation about z axis
    u_r = np.array([np.cos(theta), np.sin(theta), 0])       # radial unit vector
    u_a = np.array([-np.sin(theta), np.cos(theta), 0])      # azimuthal unit vector
    u_z = np.array([0,0,1])

    vr = np.dot(V,u_r)                      # radial velocity
    va = np.dot(V,u_a)                      # azimuthal velocity
    vz = np.dot(V,u_z)
     
    L = np.cross(R,V)  # angular momentum
    i = np.arccos(L[2]/np.linalg.norm(L)) # inclination
    ihat = i/h    
    e = np.linalg.norm((v**2/mu-1/r)*R-np.dot(R,V)/mu*V) # eccentricity
    ehat = e/h    
    
    if ida:
        # Ida prescription
        t_a = t_wave/Ct/h**2*(1+Ct/Cm*(ehat**2+ihat**2)**(1/2))
        t_e = t_wave/0.78*(1+1/15*(ehat**2+ihat**2)**(3/2))
        t_i = t_wave/0.544*(1+1/21.5*(ehat**2+ihat**2)**(3/2))
        dvdt2 = - vk/2/t_a*u_a - vr/t_e*u_r - (va-vk)/t_e*u_a - vz/t_i*u_z
    else:
        # CN prescription
        Pe = (1+(ehat/2.25)**1.2+(ehat/2.84)**6)/(1-(ehat/2.02)**4)
        t_a = t_wave*2/h**2/(2.7+1.1*0.5)*(Pe+Pe/np.abs(Pe)*(0.07*ihat+0.085*ihat**4-0.08*ehat*ihat**2))
        t_e = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3+0.18*ehat*ihat**2)
        t_i = t_wave/0.544*(1-0.30*ihat**2+0.24*ihat**3+0.14*ihat*ehat**2)    
        dvdt2 = - V/t_a - 2*np.dot(V,R)*R/r**2/t_e - vz/t_i*u_z
        
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
# %%
R = W[:,0:3]
V = W[:,3:6]
r = np.linalg.norm(R, axis=1)
v = np.linalg.norm(V, axis=1)
e_rk = np.linalg.norm([(v[i]**2/mu-1/r[i])*R[i]-(np.dot(R[i],V[i])*V[i])/mu for i in range(len(W))], axis=1)
L = np.cross(R,V)  # angular momentum
i_rk = np.arccos(L[:,2]/np.linalg.norm(L, axis=1))
E = v**2/2 - mu/r
a = -mu/2/E
# %%
ida1 = np.array(pd.read_csv('./data/datatestida.csv'))
cn1 = np.array(pd.read_csv('./data/datatestcn.csv'))
ida2 = np.array(pd.read_csv('./data/datatestida2.csv'))
cn2 = np.array(pd.read_csv('./data/datatestcn2.csv'))
idal = np.array(pd.read_csv('./data/datatestidalarge.csv'))
cnl = np.array(pd.read_csv('./data/datatestcnlarge.csv'))
# genga_p = genga_data[:,0]
# genga_t = genga_data[:,1]
# genga_a = genga_data[:,2]
# genga_e = genga_data[:,3]
# genga_i = genga_data[:,4]
total_time = int(ida1[0,1])

fig, ax = plt.subplots(1, figsize=(9,7))
ax.scatter(ida1[:,2], ida1[:,3], s=20, alpha=0.7, label='IDA')
ax.scatter(cn1[:,2], cn1[:,3], s=20, alpha=0.7, label='CN')
# ax.scatter(ida2[:,2], ida2[:,3], s=20, alpha=0.7, label='IDA')
# ax.scatter(cn2[:,2], cn2[:,3], s=20, alpha=0.7, label='CN')
ax.scatter(idal[:,2], idal[:,3], s=20, alpha=0.7, label='IDA')
ax.scatter(cnl[:,2], cnl[:,3], s=20, alpha=0.7, label='CN')
ax.set_xlabel('semi-major axis (AU)')
ax.set_ylabel('eccentricity')
ax.set_ylim(0)
ax.set_title(f'Eccentricity and semi-major axis after {total_time} years')
ax.legend()

fig, ax = plt.subplots(1, figsize=(9,7))
ax.scatter(ida1[:,2], ida1[:,4], s=20, alpha=0.7, label='IDA')
ax.scatter(cn1[:,2], cn1[:,4], s=20, alpha=0.7, label='CN')
# ax.scatter(ida2[:,2], ida2[:,4], s=20, alpha=0.7, label='IDA')
# ax.scatter(cn2[:,2], cn2[:,4], s=20, alpha=0.7, label='CN')
ax.scatter(idal[:,2], idal[:,4], s=20, alpha=0.7, label='IDA')
ax.scatter(cnl[:,2], cnl[:,4], s=20, alpha=0.7, label='CN')
ax.set_xlabel('semi-major axis (AU)')
ax.set_ylabel('inclination')
ax.set_ylim(0)
ax.set_title(f'Inclination and semi-major axis after {total_time} years')
ax.legend()
# %%
idaR = np.loadtxt('./data/Outtestidasmalli_000000000000.dat')
x = idaR[:,4]
y = idaR[:,5]
z = idaR[:,6]

fig = plt.figure(figsize=(11,11))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z, s=1, alpha=0.5)
# ax.plot(W[0::1,0], W[0::1,1], W[0::1,2], linewidth=0.005, c='steelblue')
# ax.view_init(elev=2, azim=340)
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
lim = 0.8
ax.set_xlim([-lim,lim])
ax.set_ylim([-lim,lim])
ax.set_zlim([-lim,lim])
# %%

fig, ax = plt.subplots(1, figsize=(9,6))
ax.hist(genga_ida[:,2], 20)
# %%
e_euler = np.zeros(noutputs)
i_euler = np.zeros(noutputs)

for j in range(noutputs):    
    ehat = e/h    
    ihat = i/h
    if ida:
        t_e = t_wave/0.78*(1+1/15*(ehat**2+ihat**2)**(3/2))                     # equation 34 from Ida 20
        t_i = t_wave/0.544*(1+1/21.5*(ehat**2+ihat**2)**(3/2))                  # Ida
    else:
        t_e = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3+0.18*ehat*ihat**2)     # equation 11 from Creswell+Nelson 08
        t_i = t_wave/0.544*(1-0.30*ihat**2+0.24*ihat**3+0.14*ihat*ehat**2)    # CN    
    de = -e/t_e*dt
    e = e+de
    e_euler[j] = e    
    di = -i/t_i*dt
    i = i+di
    i_euler[j] = i    
# %%
fig, ax = plt.subplots(1, figsize=(9,7))
# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_euler, label=r"$t_{ecc}$")
# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, e_rk, linestyle='--', label="RK")
# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, i_euler, c='tab:blue', label=r"$t_{inc}$")
# ax.plot(np.arange(0,noutputs,1)*dt/2/np.pi*omegak, i_rk, linestyle='-', c='tab:orange', label="RK")
ax.plot(genga_t, genga_e, linestyle='--', c='tab:red', label='genga')
ax.axhline(h, c='black', label='H/r')
ax.set_xlabel('time (years)')
ax.set_ylabel('e')
ax.set_xlim(0)
ax.set_ylim(0)
ax.set_title("Ida" if ida else "CN")
ax.tick_params(which='both', direction="in", top=True, right=True)
ax.grid()
ax.legend()

# fig.savefig('/home/john/summerproject/img/test1.png', bbox_inches='tight')
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
ax.scatter(W[-1,1], W[-1,2], c='steelblue')
ax.plot(W[:,1], W[:,2], linewidth=0.05, c='steelblue')

# plt.savefig('/home/john/Desktop/summerproject/img/5000years.png', bbox_inches='tight')
# %%
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
lim = 0.8
ax.scatter(0,0,0, c='black')
ax.scatter(W[-1,0], W[-1,1], W[-1,2], c='steelblue')
ax.plot(W[0::1,0], W[0::1,1], W[0::1,2], linewidth=0.005, c='steelblue')
ax.view_init(elev=2, azim=340)
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_xlim([-lim,lim])
ax.set_ylim([-lim,lim])
ax.set_zlim([-lim,lim])

# plt.savefig(f'close-{particles}-{runtime}-{ts}.png', bbox_inches='tight')
# %%
files = glob(f'./Out*.dat')
files.sort()
data = np.concatenate([np.loadtxt(i) for i in files])

planet_i = 1

planet1 = data[data[:,1] == planet_i]

