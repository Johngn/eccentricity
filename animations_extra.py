#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:01:12 2020

@author: johngillan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
# %%
data_ida = np.array(pd.read_csv('./data/dataida_1body.csv'))
data_cn = np.array(pd.read_csv('./data/datacn_1body.csv'))

time_steps = 50001
n_planets = 1

data_ida_new = np.ones((time_steps*n_planets, 6))*-1
data_cn_new = np.ones((time_steps*n_planets, 6))*-1

ida_mask = np.zeros(time_steps*n_planets, dtype=bool)
cn_mask = np.zeros(time_steps*n_planets, dtype=bool)

for i in range(n_planets):
    planet0_mask = data_ida[:,0] == i
    mask_len = np.count_nonzero(planet0_mask)
    
    data_mask = np.zeros(time_steps*n_planets, dtype=bool)
    data_mask[i:n_planets*mask_len:n_planets] = True
    data_ida_new[data_mask] = data_ida[planet0_mask]
    
    planet0_mask = data_cn[:,0] == i
    mask_len = np.count_nonzero(planet0_mask)
    
    data_mask = np.zeros(time_steps*n_planets, dtype=bool)
    data_mask[i:n_planets*mask_len:n_planets] = True
    data_cn_new[data_mask] = data_cn[planet0_mask]
    
data_ida = np.reshape(data_ida_new, (-1, n_planets, 6))[::100]
data_cn = np.reshape(data_cn_new, (-1, n_planets, 6))[::100]
# %%
# data_nd = np.array(pd.read_csv('./data/datano_damping_2.csv'))

# time_steps = 50006
# n_planets = 100

# data_nd_new = np.ones((time_steps*n_planets, 6))*-1
# nd_mask = np.zeros(time_steps*n_planets, dtype=bool)

# for i in range(n_planets):
    
#     planet0_mask = data_nd[:,0] == i
#     mask_len = np.count_nonzero(planet0_mask)
    
#     data_mask = np.zeros(time_steps*n_planets, dtype=bool)
#     data_mask[i:n_planets*mask_len:n_planets] = True
#     data_nd_new[data_mask] = data_nd[planet0_mask]
    
# data_nd = np.reshape(data_nd_new, (-1, n_planets, 6))[::100]
# %%
solar_to_earth_mass = 1.99e30/5.97e24
times = data_cn[:,0,1]

fig, axes = plt.subplots(1, figsize=(8, 4))
# fig.subplots_adjust(hspace=0)
axes.set_xlim(-0,8)
# axes[0].set_ylim(0,0.35)
axes.set_ylabel("eccentricity")

# axes[1].set_xlim(-0,8)
# # axes[1].set_ylim(0,0.2)
# axes[1].set_ylabel("inclination")

# axes[2].set_xlim(-0,8)
# # axes[2].set_ylim(0,2.5)
# axes[2].set_xlabel("semi-major axis")
# axes[2].set_ylabel("mass")

ms = 10
alpha = 0.9

planets_ida_e, = axes.plot([], [], "o", ms=ms, alpha=alpha, label="Ida")
planets_cn_e, = axes.plot([], [], "o", ms=ms, alpha=alpha, label="CN")
# planets_nd_e, = axes[0].plot([], [], "o", ms=ms, alpha=alpha, label="no damping")

# planets_ida_i, = axes[1].plot([], [], "o", ms=ms, alpha=alpha, label="Ida")
# planets_cn_i, = axes[1].plot([], [], "o", ms=ms, alpha=alpha, label="CN")
# # planets_nd_i, = axes[1].plot([], [], "o", ms=ms, alpha=alpha, label="no damping")

# planets_ida_m, = axes[2].plot([], [], "o", ms=ms, alpha=alpha, label="Ida")
# planets_cn_m, = axes[2].plot([], [], "o", ms=ms, alpha=alpha, label="CN")
# planets_nd_m, = axes[2].plot([], [], "o", ms=ms, alpha=alpha, label="no damping")

text = axes.text(2.2, .8, '', fontsize=15)

axes.legend()
# axes[1].legend()
# axes[2].legend()

def animate(i):
    planets_ida_e.set_data(data_ida[i,:,3], data_ida[i,:,4])
    planets_cn_e.set_data(data_cn[i,:,3], data_cn[i,:,4])
    # planets_nd_e.set_data(data_nd[i,:,3], data_nd[i,:,4])
    
    planets_ida_i.set_data(data_ida[i,:,3], data_ida[i,:,5])
    planets_cn_i.set_data(data_cn[i,:,3], data_cn[i,:,5])
    # planets_nd_i.set_data(data_nd[i,:,3], data_nd[i,:,5])
    
    planets_ida_m.set_data(data_ida[i,:,3], data_ida[i,:,2]*solar_to_earth_mass)
    planets_cn_m.set_data(data_cn[i,:,3], data_cn[i,:,2]*solar_to_earth_mass)
    # planets_nd_m.set_data(data_nd[i,:,3], data_nd[i,:,2]*solar_to_earth_mass)
    
    text.set_text(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    return planets_cn_e, planets_cn_i
    
im_ani = FuncAnimation(fig, animate, frames=len(data_ida), interval=1)
# %%
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
f = f'vid/genga_scatter_a1.mp4' 
writervideo = FFMpegWriter(fps=30) # ffmpeg must be installed
im_ani.save(f, writer=writervideo)
# %%
data_all = [data_ida, data_cn, data_nd]

fig, axes = plt.subplots(3, figsize=(8,12))
fig.subplots_adjust(hspace=0.4)

def animate(i, data_all):
    axes[0].cla()
    axes[0].set_xlabel("semi-major axis")
    axes[0].set_ylim(0,30)
    axes[0].hist([data_all[0][i,:,3], data_all[1][i,:,3], data_all[2][i,:,3]], bins=40, range=(0,7), label=['IDA', 'CN', 'no damping'])
    axes[0].set_title(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    axes[0].legend()
    
    axes[1].cla()
    axes[1].set_xlabel("eccentricity")
    axes[1].set_ylim(0,30)
    axes[1].hist([data_all[0][i,:,4], data_all[1][i,:,4], data_all[2][i,:,4]], bins=40, range=(0,0.06))
    
    axes[2].cla()
    axes[2].set_xlabel("inclination")
    axes[2].set_ylim(0,30)
    axes[2].hist([data_all[0][i,:,5], data_all[1][i,:,5], data_all[2][i,:,5]], bins=40, range=(0,0.06))

animation = FuncAnimation(fig, animate, frames=len(data_ida), fargs=(data_all,), interval=1)
# %%
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
writervideo = FFMpegWriter(fps=30) # ffmpeg must be installed
animation.save('vid/genga_hist_3.mp4', writer=writervideo)

# animation.save('vid/genga_hist_4.gif', writer='imagemagick', fps=60)
