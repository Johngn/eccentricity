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
import seaborn as sns
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

ida_1 = orbitalelements(np.loadtxt('./data/Outida1_000100000000.dat'))
ida_2 = orbitalelements(np.loadtxt('./data/Outida2_000100000000.dat'))
ida_3 = orbitalelements(np.loadtxt('./data/Outida3_000100000000.dat'))
# ida_4 = orbitalelements(np.loadtxt('./data/Outida4_000100000000.dat'))
# ida_5 = orbitalelements(np.loadtxt('./data/Outida5_000100000000.dat'))

cn_1 = orbitalelements(np.loadtxt('./data/Outcn_mars_4_000100000000.dat'))
cn_2 = orbitalelements(np.loadtxt('./data/Outcn_mars_5_000100000000.dat'))
cn_3 = orbitalelements(np.loadtxt('./data/Outcn_mars_6_000100000000.dat'))
cn_4 = orbitalelements(np.loadtxt('./data/Outcn_7_000100000000.dat'))
cn_5 = orbitalelements(np.loadtxt('./data/Outcn_8_000100000000.dat'))

no_damping_1 = orbitalelements(np.loadtxt('./data/Outno_damping_2_000100000000.dat'))

# convert to log values for total kde plot below
# no_damping_1 = np.log10(no_damping_1)
# %%
ida_data = ida_1
cn_data = cn_1

fig, ax = plt.subplots(3, figsize=(9,10))

ax[0].scatter(ida_data[:,0], ida_data[:,1], s=ida_data[:,3]*100, alpha=0.7, label='IDA')
ax[0].scatter(cn_data[:,0], cn_data[:,1], s=cn_data[:,3]*100, alpha=0.7, label='CN')
ax[0].scatter(no_damping_1[:,0], no_damping_1[:,1], s=no_damping_1[:,3]*100, alpha=0.7, label='no damping')
ax[0].set_ylabel('eccentricity')
ax[0].set_ylim(0)
ax[0].legend()

ax[1].scatter(ida_data[:,0], ida_data[:,2], s=ida_data[:,3]*100, alpha=0.7, label='IDA')
ax[1].scatter(cn_data[:,0], cn_data[:,2], s=cn_data[:,3]*100, alpha=0.7, label='CN')
ax[1].scatter(no_damping_1[:,0], no_damping_1[:,2], s=no_damping_1[:,3]*100, alpha=0.7, label='no damping')
ax[1].set_ylabel('inclination')
ax[1].set_ylim(0)
ax[1].legend()

ax[2].scatter(ida_data[:,0], ida_data[:,3], s=ida_data[:,3]*100, alpha=0.7, label='IDA')
ax[2].scatter(cn_data[:,0], cn_data[:,3], s=cn_data[:,3]*100, alpha=0.7, label='CN')
ax[2].scatter(no_damping_1[:,0], no_damping_1[:,3], s=no_damping_1[:,3]*100, alpha=0.7, label='no damping')
ax[2].set_xlabel('semi-major axis (AU)')
ax[2].set_ylabel(r'mass (${M_e}$)')
ax[2].set_ylim(0) 
ax[2].legend()

# fig.savefig('/home/john/summerproject/img/genga_scatter.pdf', bbox_inches='tight')
# %%
ida_data = [ida_1, ida_2, ida_3, ida_3, ida_3]
cn_data = [cn_1, cn_2, cn_3, cn_4, cn_5]

for i, item in enumerate(ida_data):
    ida_data[i] = np.log10(item)
    
for i, item in enumerate(cn_data):
    cn_data[i] = np.log10(item)
    
ida_data_x = np.concatenate(ida_data)
cn_data_x = np.concatenate(cn_data)

fig, ax = plt.subplots(3, figsize=(9,10))
fig.subplots_adjust(hspace=0.3)

bins1 = np.linspace(0,3,20)
bins2 = np.linspace(0,0.1,20)
bins3 = np.linspace(0,0.02,20)

kde = True
hist = False

# multiple kde plots
for i, item in enumerate(ida_data):
    sns.distplot(item[:,0], hist=hist, kde=kde, color='tab:blue', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[0])
    sns.distplot(item[:,1], hist=hist, kde=kde, color='tab:blue', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[1])
    sns.distplot(item[:,2], hist=hist, kde=kde, color='tab:blue', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[2])
for i, item in enumerate(cn_data):
    sns.distplot(item[:,0], hist=hist, kde=kde, color='tab:orange', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[0])
    sns.distplot(item[:,1], hist=hist, kde=kde, color='tab:orange', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[1])
    sns.distplot(item[:,2], hist=hist, kde=kde, color='tab:orange', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[2])

# total kde and histograms
sns.distplot(ida_data_x[:,0], color='tab:blue', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[0], label=['IDA'])
sns.distplot(ida_data_x[:,1], color='tab:blue', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[1])
sns.distplot(ida_data_x[:,2], color='tab:blue', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[2])
sns.distplot(cn_data_x[:,0], color='tab:orange', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[0], label=['CN'])
sns.distplot(cn_data_x[:,1], color='tab:orange', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[1])
sns.distplot(cn_data_x[:,2], color='tab:orange', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[2])
sns.distplot(no_damping_1[:,0], color='tab:green', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[0], label='no damping')
sns.distplot(no_damping_1[:,1], color='tab:green', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[1])
sns.distplot(no_damping_1[:,2], color='tab:green', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[2])

ax[0].set_xlabel('log(a/AU)')
ax[1].set_xlabel('log(e)')
ax[1].axvline(np.log10(0.05), linestyle='--', c='black')
ax[2].set_xlabel('log(i)')
ax[2].axvline(np.log10(0.05), linestyle='--', c='black')

# fig.savefig('/home/john/summerproject/img/dist_sum.pdf', bbox_inches='tight')
# %%
# 3d plot of system
data = np.loadtxt('./data/Outida1_000010000000.dat')

fig = plt.figure(figsize=(9,9))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(0,0,0, s=200, color="gold")
ax.scatter(data[:,4], data[:,5], data[:,6], s=20)
# ax.view_init(elev=2, azim=340)
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
lim = 3
ax.set_xlim([-lim,lim])
ax.set_ylim([-lim,lim])
ax.set_zlim([-lim,lim])
# %%
# difference between t_a for different simulations
data_ida = np.array(pd.read_csv('./data/dataida1.csv'))
data_cn = np.array(pd.read_csv('./data/datacn_mars_3.csv'))

fig, axes = plt.subplots(1, figsize=(8, 3))
axes.plot(data_ida[:,1], data_cn[:,3] - data_ida[:,3])
# axes.plot(data_cn[:,1], data_cn[:,3])
axes.set_xlabel('years')
axes.set_ylabel('AU')