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
from timeit import default_timer as timed
from glob import glob
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