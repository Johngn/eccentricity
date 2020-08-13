#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:44:56 2020

@author: johngillan
"""

# %%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

k = 0.0172020989
G = k**2
days_per_year = 365.25
mu = G
solar_to_earth_mass = 1.99e30/5.97e24

def orbitalelements(data):
    m = [data[i][2]*solar_to_earth_mass for i, item in enumerate(data)]
    t = [data[i][0] for i, item in enumerate(data)]
    R = [data[i][4:7] for i, item in enumerate(data)]
    V = [data[i][7:10]*k for i, item in enumerate(data)]
    r = np.linalg.norm(R, axis=1)
    v = np.linalg.norm(V, axis=1)
    L = np.cross(R,V)
    i = np.arccos(L[:,2]/np.linalg.norm(L, axis=1))
    e = np.linalg.norm([(v[i]**2/mu - 1/r[i])*R[i] - np.dot(R[i],V[i])*V[i]/mu for i, item in enumerate(data)], axis=1)
    E = v**2/2 - mu/r
    a = -mu/2/E
    
    a = np.reshape(a, (len(a), 1))
    e = np.reshape(e, (len(e), 1))
    i = np.reshape(i, (len(i), 1))
    m = np.reshape(m, (len(m), 1))
    
    return np.hstack((a,e,i,m))

idal = orbitalelements(np.loadtxt('./data/Outidal_000010000000.dat'))
cnl = orbitalelements(np.loadtxt('./data/Outcnl_000010000000.dat'))

ida_300 = orbitalelements(np.loadtxt('./data/Outida_300_000010000000.dat'))
cn_300 = orbitalelements(np.loadtxt('./data/Outcn_300_000010000000.dat'))

ida_300_longer = orbitalelements(np.loadtxt('./data/Outida_300_longer_000030000000.dat'))
cn_300_longer = orbitalelements(np.loadtxt('./data/Outcn_300_longer_000030000000.dat'))

ida_mars = orbitalelements(np.loadtxt('./data/Outida_mars_000100000000.dat'))
cn_mars = orbitalelements(np.loadtxt('./data/Outcn_mars_000100000000.dat'))

ida_mars_2 = orbitalelements(np.loadtxt('./data/Outida_mars_2_000100000000.dat'))
cn_mars_2 = orbitalelements(np.loadtxt('./data/Outcn_mars_2_000100000000.dat'))

ida_mars_3 = orbitalelements(np.loadtxt('./data/Outida_mars_3_000100000000.dat'))
cn_mars_3 = orbitalelements(np.loadtxt('./data/Outcn_mars_3_000100000000.dat'))

ida_mars_4 = orbitalelements(np.loadtxt('./data/Outida_mars_4_000100000000.dat'))
cn_mars_4 = orbitalelements(np.loadtxt('./data/Outcn_mars_4_000100000000.dat'))
# %%
ida_data = ida_mars_4
cn_data = cn_mars_4

fig, ax = plt.subplots(2, figsize=(9,9))

ax[0].scatter(ida_data[:,0], ida_data[:,1], s=ida_data[:,3]*200, alpha=0.7, label='IDA')
ax[0].scatter(cn_data[:,0], cn_data[:,1], s=cn_data[:,3]*200, alpha=0.7, label='CN')
ax[0].set_ylabel('eccentricity')
ax[0].set_ylim(0)
ax[0].set_title(r'${10^6}$ years')
ax[0].legend()

ax[1].scatter(ida_data[:,0], ida_data[:,2], s=ida_data[:,3]*200, alpha=0.7, label='IDA')
ax[1].scatter(cn_data[:,0], cn_data[:,2], s=cn_data[:,3]*200, alpha=0.7, label='CN')
ax[1].set_ylabel('inclination')
ax[1].set_ylim(0)
ax[1].legend()

# ax[2].scatter(ida_data[:,0], ida_data[:,3], s=ida_data[:,3]*200, alpha=0.7, label='IDA')
# ax[2].scatter(cn_data[:,0], cn_data[:,3], s=cn_data[:,3]*200, alpha=0.7, label='CN')
# ax[2].set_xlabel('semi-major axis (AU)')
# ax[2].set_ylabel(r'mass (${M_e}$)')
# ax[2].set_ylim(0)
# ax[2].legend()

# fig.savefig('/home/john/summerproject/img/genga_mars_3.png', bbox_inches='tight')
# %%
fig, ax = plt.subplots(3, figsize=(8,12))
fig.subplots_adjust(hspace=0.4)

ax[0].hist([ida_data[:,0], cn_data[:,0]], 20, label=['IDA', 'CN'])
ax[0].set_xlabel('semi-major axis')
ax[0].legend(loc='upper right')

ax[1].hist([ida_data[:,1], cn_data[:,1]], 20, label=['IDA', 'CN'])
ax[1].set_xlabel('eccentricity')
ax[1].legend(loc='upper right')

ax[2].hist([ida_data[:,2], cn_data[:,2]], 20, label=['IDA', 'CN'])
ax[2].set_xlabel('inclination')
ax[2].legend(loc='upper right')

# fig.savefig('/home/john/summerproject/img/hist_mars.png', bbox_inches='tight')
# %%
data = np.loadtxt('./data/Outida_300_000010000000.dat')

fig = plt.figure(figsize=(9,9))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(0,0,0, s=200, color="gold")
ax.scatter(data[:,4], data[:,5], data[:,6], s=20)
# ax.view_init(elev=2, azim=340)
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
lim = 1.5
ax.set_xlim([-lim,lim])
ax.set_ylim([-lim,lim])
ax.set_zlim([-lim,lim])