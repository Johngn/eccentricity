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

def orbitalelements(data):
    p = [data[i][1] for i, item in enumerate(data)]
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
    
    return np.hstack((a,e,i))

ida = orbitalelements(np.loadtxt('./data/Outida_000001000000.dat'))
idal = orbitalelements(np.loadtxt('./data/Outidal_000010000000.dat'))
ida_300 = orbitalelements(np.loadtxt('./data/Outida_300_000010000000.dat'))
cn_300 = orbitalelements(np.loadtxt('./data/Outcn_300_000010000000.dat'))
cnl = orbitalelements(np.loadtxt('./data/Outcnl_000010000000.dat'))
# %%
fig, ax = plt.subplots(1, figsize=(9,7))
# ax.scatter(ida[:,0], ida[:,1], s=20, alpha=0.7, label='IDA')
ax.scatter(ida_300[:,0], ida_300[:,1], s=40, alpha=0.8, label='IDA')
ax.scatter(cn_300[:,0], cn_300[:,1], s=40, alpha=0.8, label='CN')
ax.set_xlabel('semi-major axis (AU)')
ax.set_ylabel('eccentricity')
ax.set_ylim(0)
ax.legend()
fig.savefig('/home/john/summerproject/img/genga_300_e.png', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, figsize=(9,7))
# ax.scatter(ida[:,0], ida[:,2], s=20, alpha=0.7, label='IDA')
ax.scatter(ida_300[:,0], ida_300[:,2], s=40, alpha=0.8, label='IDA')
ax.scatter(cn_300[:,0], cn_300[:,2], s=40, alpha=0.8, label='CN')
ax.set_xlabel('semi-major axis (AU)')
ax.set_ylabel('inclination')
ax.set_ylim(0)
ax.legend()
fig.savefig('/home/john/summerproject/img/genga_300_i.png', bbox_inches='tight')
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
# %%
fig, ax = plt.subplots(1, figsize=(9,6))
ax.hist([ida_300[:,1], cn_300[:,1]], 20, label=['IDA', 'CN'])
ax.set_xlabel('eccentricity')
plt.legend(loc='upper right')
fig.savefig('/home/john/summerproject/img/e_hist_300.png', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, figsize=(9,6))
ax.hist([ida_300[:,2], cn_300[:,2]], 20, label=['IDA', 'CN'])
ax.set_xlabel('inclination')
plt.legend(loc='upper right')
fig.savefig('/home/john/summerproject/img/i_hist_300.png', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, figsize=(9,6))
ax.hist([ida_300[:,0], cn_300[:,0]], 20, label=['IDA', 'CN'])
ax.set_xlabel('semi-major axis')
plt.legend(loc='upper right')
fig.savefig('/home/john/summerproject/img/a_hist_300.png', bbox_inches='tight')