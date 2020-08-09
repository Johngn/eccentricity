#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:44:56 2020

@author: johngillan
"""

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import axes3d
# %%
data = np.loadtxt('./data/Outida_000001000000.dat')

k = 0.0172020989
G = k**2
days_per_year = 365.25

mass = data[0][2]
mu = G

p = [data[i][1] for i in range(len(data))]
t = [data[i][0] for i in range(len(data))]
R = [data[i][4:7] for i in range(len(data))]
V = [data[i][7:10]*k for i in range(len(data))]
r = np.linalg.norm(R, axis=1)
v = np.linalg.norm(V, axis=1)
L = np.cross(R,V)
i = np.arccos(L[:,2]/np.linalg.norm(L, axis=1))
e = np.linalg.norm([(v[i]**2/mu - 1/r[i])*R[i] - np.dot(R[i],V[i])*V[i]/mu for i in range(len(data))], axis=1)
E = v**2/2 - mu/r
a = -mu/2/E

fig, ax = plt.subplots(1, figsize=(9,7))
ax.scatter(a, e, s=20, alpha=0.7, label='IDA')
# ax.scatter(cn1[:,2], cn1[:,3], s=20, alpha=0.7, label='CN')
# ax.scatter(ida2[:,2], ida2[:,3], s=20, alpha=0.7, label='IDA')
# ax.scatter(cn2[:,2], cn2[:,3], s=20, alpha=0.7, label='CN')
# ax.scatter(idal[:,2], idal[:,3], s=20, alpha=0.7, label='IDA')
# ax.scatter(cnl[:,2], cnl[:,3], s=20, alpha=0.7, label='CN')
ax.set_xlabel('semi-major axis (AU)')
ax.set_ylabel('eccentricity')
ax.set_ylim(0)
ax.set_title(f'Eccentricity and semi-major axis after {np.round(t[-1],0)} years')
ax.legend()

fig, ax = plt.subplots(1, figsize=(9,7))
ax.scatter(a, i, s=20, alpha=0.7, label='IDA')
# ax.scatter(cn1[:,2], cn1[:,4], s=20, alpha=0.7, label='CN')
# ax.scatter(ida2[:,2], ida2[:,4], s=20, alpha=0.7, label='IDA')
# ax.scatter(cn2[:,2], cn2[:,4], s=20, alpha=0.7, label='CN')
# ax.scatter(idal[:,2], idal[:,4], s=20, alpha=0.7, label='IDA')
# ax.scatter(cnl[:,2], cnl[:,4], s=20, alpha=0.7, label='CN')
ax.set_xlabel('semi-major axis (AU)')
ax.set_ylabel('inclination')
ax.set_ylim(0)
ax.set_title(f'Inclination and semi-major axis after {np.round(t[-1],0)} years')
ax.legend()
# %%
idaR = np.loadtxt('./data/Outida2_000000000000.dat')
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
lim = 1
ax.set_xlim([-lim,lim])
ax.set_ylim([-lim,lim])
ax.set_zlim([-lim,lim])
# %%

fig, ax = plt.subplots(1, figsize=(9,6))
ax.hist(i, 20)