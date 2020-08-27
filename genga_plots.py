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

ida_4 = orbitalelements(np.loadtxt('./data/Outida_mars_4_000100000000.dat'))
cn_4 = orbitalelements(np.loadtxt('./data/Outcn_mars_4_000100000000.dat'))

ida_5 = orbitalelements(np.loadtxt('./data/Outida_mars_5_000100000000.dat'))
cn_5 = orbitalelements(np.loadtxt('./data/Outcn_mars_5_000100000000.dat'))

# ida_6 = orbitalelements(np.loadtxt('./data/Outida_mars_6_000100000000.dat'))
# cn_6 = orbitalelements(np.loadtxt('./data/Outcn_mars_6_000100000000.dat'))

ida_7 = orbitalelements(np.loadtxt('./data/Outida_7_000100000000.dat'))
cn_7 = orbitalelements(np.loadtxt('./data/Outcn_7_000100000000.dat'))

ida_8 = orbitalelements(np.loadtxt('./data/Outida_8_000100000000.dat'))
cn_8 = orbitalelements(np.loadtxt('./data/Outcn_8_000100000000.dat'))

ida_9 = orbitalelements(np.loadtxt('./data/Outida_9_000100000000.dat'))
cn_9 = orbitalelements(np.loadtxt('./data/Outcn_9_000100000000.dat'))



# ida_4 = orbitalelements(np.loadtxt('./data/Outida_mars_4_000000000000.dat'))
# cn_4 = orbitalelements(np.loadtxt('./data/Outcn_mars_4_000000000000.dat'))

# ida_5 = orbitalelements(np.loadtxt('./data/Outida_mars_5_000000000000.dat'))
# cn_5 = orbitalelements(np.loadtxt('./data/Outcn_mars_5_000000000000.dat'))

# # ida_6 = orbitalelements(np.loadtxt('./data/Outida_mars_6_000000000000.dat'))
# # cn_6 = orbitalelements(np.loadtxt('./data/Outcn_mars_6_000000000000.dat'))

# ida_7 = orbitalelements(np.loadtxt('./data/Outida_7_000000000000.dat'))
# cn_7 = orbitalelements(np.loadtxt('./data/Outcn_7_000000000000.dat'))

# ida_8 = orbitalelements(np.loadtxt('./data/Outida_8_000000000000.dat'))
# cn_8 = orbitalelements(np.loadtxt('./data/Outcn_8_000000000000.dat'))

# ida_9 = orbitalelements(np.loadtxt('./data/Outida_9_000000000000.dat'))
# cn_9 = orbitalelements(np.loadtxt('./data/Outcn_9_000000000000.dat'))

no_damping_1 = orbitalelements(np.loadtxt('./data/Outno_damping_2_000100000000.dat'))
no_damping_1 = np.log10(no_damping_1)

ida_data = [ida_4, ida_5, ida_7, ida_8, ida_9]
cn_data = [cn_4, cn_5, cn_7, cn_8, cn_9]

for i, item in enumerate(ida_data):
    ida_data[i] = np.log10(item)
    # ida_data[i] = np.hstack((np.ones((len(item), 1))*i, item))
    
for i, item in enumerate(cn_data):
    cn_data[i] = np.log10(item)
    # cn_data[i] = np.hstack((np.ones((len(item), 1))*i, item))
    
ida_data_x = np.concatenate(ida_data)
cn_data_x = np.concatenate(cn_data)

df1 = pd.DataFrame(ida_4)
df2 = pd.DataFrame(cn_4)

df1.to_csv('./data/ida_collision.csv', index=False)
df2.to_csv('./data/cn_collision.csv', index=False)
# %%
ida_data = ida_8
cn_data = cn_8

fig, ax = plt.subplots(3, figsize=(9,10))

ax[0].scatter(ida_data[:,0], ida_data[:,1], s=ida_data[:,3]*100, alpha=0.7, label='IDA')
ax[0].scatter(cn_data[:,0], cn_data[:,1], s=cn_data[:,3]*100, alpha=0.7, label='CN')
ax[0].scatter(no_damping_1[:,0], no_damping_1[:,1], s=no_damping_1[:,3]*200, alpha=0.7, label='no damping')
ax[0].set_ylabel('eccentricity')
ax[0].set_ylim(0)
# ax[0].set_title(r'${10^6}$ years')
ax[0].legend()

ax[1].scatter(ida_data[:,0], ida_data[:,2], s=ida_data[:,3]*100, alpha=0.7, label='IDA')
ax[1].scatter(cn_data[:,0], cn_data[:,2], s=cn_data[:,3]*100, alpha=0.7, label='CN')
ax[1].scatter(no_damping_1[:,0], no_damping_1[:,2], s=no_damping_1[:,3]*200, alpha=0.7, label='no damping')
ax[1].set_ylabel('inclination')
ax[1].set_ylim(0)
ax[1].legend()

ax[2].scatter(ida_data[:,0], ida_data[:,3], s=ida_data[:,3]*100, alpha=0.7, label='IDA')
ax[2].scatter(cn_data[:,0], cn_data[:,3], s=cn_data[:,3]*100, alpha=0.7, label='CN')
ax[2].scatter(no_damping_1[:,0], no_damping_1[:,3], s=no_damping_1[:,3]*200, alpha=0.7, label='no damping')
ax[2].set_xlabel('semi-major axis (AU)')
ax[2].set_ylabel(r'mass (${M_e}$)')
ax[2].set_ylim(0)
ax[2].legend()

# fig.savefig('/home/john/summerproject/img/genga_scatterL.pdf', bbox_inches='tight')
# %%
ida_histograms = []
cn_histograms = []

for i, item in enumerate(ida_data):
    ida_histograms.append(np.histogram(item[:,1], bins=20, range=[0,3])[0])
    ranges = np.array(np.histogram(item[:,1], bins=20, range=[0,3])[1])
for i, item in enumerate(cn_data):
    cn_histograms.append(np.histogram(item[:,1], bins=20, range=[0,3])[0])
    ranges = np.array(np.histogram(item[:,1], bins=20, range=[0,3])[1])

ida_histograms = np.array(ida_histograms)
cn_histograms = np.array(cn_histograms)
# std1 = np.std(histograms, axis=0)
# means = np.mean(histograms, axis=0)

fig, ax = plt.subplots(1, figsize=(10,8))
for i, item in enumerate(ida_histograms):
    plt.plot(ranges[:-1], item)
# plt.scatter(ranges[:-1], ida_histograms[4])
# %%
fig, ax = plt.subplots(3, figsize=(9,10))
fig.subplots_adjust(hspace=0.3)

bins1 = np.linspace(0,3,20)
bins2 = np.linspace(0,0.1,20)
bins3 = np.linspace(0,0.02,20)

kde = True
hist = False

# for i, item in enumerate(ida_data):
#     sns.distplot(item[:,0], hist=hist, kde=kde, color='tab:blue', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[0])
#     sns.distplot(item[:,1], hist=hist, kde=kde, color='tab:blue', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[1])
#     sns.distplot(item[:,2], hist=hist, kde=kde, color='tab:blue', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[2])
# for i, item in enumerate(cn_data):
#     sns.distplot(item[:,0], hist=hist, kde=kde, color='tab:orange', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[0])
#     sns.distplot(item[:,1], hist=hist, kde=kde, color='tab:orange', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[1])
#     sns.distplot(item[:,2], hist=hist, kde=kde, color='tab:orange', hist_kws={"histtype": "step", "linewidth": 3, "alpha": 1}, ax=ax[2])
    # ax[0].hist([ida_data[i][:,0], cn_data[i][:,0]], 20, color=['blue','orange'], rwidth=0.3, alpha=0.8)    
    # ax[0].hist(ida_data[i][:,1], 20, color='blue', rwidth=0.3, alpha=0.8)
    # ax[0].hist(cn_data[i][:,1], 20, color='orange', rwidth=0.3, alpha=0.8)
    # ax[1].hist([ida_data[i][:,1], cn_data[i][:,1]], 20, color=['blue','orange'], rwidth=0.4, alpha=0.8)
    # ax[2].hist([ida_data[i][:,2], cn_data[i][:,2]], 20, color=['blue','orange'], rwidth=0.4, alpha=0.8)
    

 

sns.distplot(ida_data_x[:,0], color='tab:blue', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[0], label=['IDA'])
sns.distplot(ida_data_x[:,1], color='tab:blue', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[1])
sns.distplot(ida_data_x[:,2], color='tab:blue', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[2])
sns.distplot(cn_data_x[:,0], color='tab:orange', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[0], label=['CN'])
sns.distplot(cn_data_x[:,1], color='tab:orange', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[1])
sns.distplot(cn_data_x[:,2], color='tab:orange', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[2])
sns.distplot(no_damping_1[:,0], color='tab:green', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[0], label='no damping')
sns.distplot(no_damping_1[:,1], color='tab:green', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[1])
sns.distplot(no_damping_1[:,2], color='tab:green', hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1}, ax=ax[2])

# ax[0].hist([ida_data_x[:,0], cn_data_x[:,0]], 20, color=['tab:blue','tab:orange'])
# ax[1].hist([ida_data_x[:,1], cn_data_x[:,1]], 20, color=['tab:blue','tab:orange'])
# ax[2].hist([ida_data_x[:,2], cn_data_x[:,2]], 20, color=['tab:blue','tab:orange'])


ax[0].set_xlabel('log(a/AU)')
# ax[0].set_xlim(0)
ax[0].legend(loc='upper left')

ax[1].set_xlabel('log(e)')
ax[1].axvline(np.log10(0.05), linestyle='--', c='black')
# ax[1].axvline(0.05, linestyle='--', c='black')

# ax[2].hist([ida_4[:,2], cn_4[:,2]], 20, label=['IDA', 'CN'])
# ax[2].hist([ida_4[:,2], cn_4[:,2]], 20, label=['IDA', 'CN'])
ax[2].set_xlabel('log(i)')
ax[2].axvline(np.log10(0.05), linestyle='--', c='black')
# ax[2].axvline(0.05, linestyle='--', c='black')
# ax[2].set_ylim(0,3000)

fig.savefig('/home/john/summerproject/img/dist_sum.pdf', bbox_inches='tight')
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