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
data_ida = np.array(pd.read_csv('./data/dataida_mars_3.csv'))
data_cn = np.array(pd.read_csv('./data/datacn_mars_3.csv'))

data_ida_new = np.ones((5050101, 6))*-1
data_cn_new = np.ones((5050101, 6))*-1

ida_mask = np.zeros(5050101, dtype=bool)
cn_mask = np.zeros(5050101, dtype=bool)

for i in range(101):
    planet0_mask = data_ida[:,0] == i
    mask_len = np.count_nonzero(planet0_mask)
    
    data_mask = np.zeros(5050101, dtype=bool)
    data_mask[i:101*mask_len:101] = True
    data_ida_new[data_mask] = data_ida[planet0_mask]
    
    planet0_mask = data_cn[:,0] == i
    mask_len = np.count_nonzero(planet0_mask)
    
    data_mask = np.zeros(5050101, dtype=bool)
    data_mask[i:101*mask_len:101] = True
    data_cn_new[data_mask] = data_cn[planet0_mask]
    
data_ida = np.reshape(data_ida_new, (-1, 101, 6))[::100]
data_cn = np.reshape(data_cn_new, (-1, 101, 6))[::100]
# %%
solar_to_earth_mass = 1.99e30/5.97e24
times = data_ida[:,0,1]

fig, axes = plt.subplots(3, figsize=(8, 10))
# fig.subplots_adjust(hspace=0)
axes[0].set_xlim(-0,7)
axes[0].set_ylim(0,0.15)
axes[0].set_ylabel("eccentricity")

axes[1].set_xlim(-0,7)
axes[1].set_ylim(0,0.15)
axes[1].set_ylabel("inclination")

axes[2].set_xlim(-0,7)
axes[2].set_ylim(0,2.5)
axes[2].set_xlabel("semi-major axis")
axes[2].set_ylabel("mass")

planets_ida_e, = axes[0].plot([], [], "o", markersize=8, alpha=0.8, label="Ida")
planets_cn_e, = axes[0].plot([], [], "o", markersize=8, alpha=0.8, label="CN")

planets_ida_i, = axes[1].plot([], [], "o", markersize=8, alpha=0.8, label="Ida")
planets_cn_i, = axes[1].plot([], [], "o", markersize=8, alpha=0.8, label="CN")

planets_ida_m, = axes[2].plot([], [], "o", markersize=8, alpha=0.8, label="Ida")
planets_cn_m, = axes[2].plot([], [], "o", markersize=8, alpha=0.8, label="CN")
text = axes[0].text(2.2, 0.12, '', fontsize=15)

axes[0].legend()
axes[1].legend()
axes[2].legend()

def animate(i):
    planets_ida_e.set_data(data_ida[i,:,3], data_ida[i,:,4])
    planets_cn_e.set_data(data_cn[i,:,3], data_cn[i,:,4])
    
    planets_ida_i.set_data(data_ida[i,:,3], data_ida[i,:,5])
    planets_cn_i.set_data(data_cn[i,:,3], data_cn[i,:,5])
    
    planets_ida_m.set_data(data_ida[i,:,3], data_ida[i,:,2]*solar_to_earth_mass)
    planets_cn_m.set_data(data_cn[i,:,3], data_cn[i,:,2]*solar_to_earth_mass)
    
    text.set_text(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    return planets_ida_e, planets_cn_e, planets_ida_i, planets_cn_i
    
im_ani = FuncAnimation(fig, animate, frames=len(data_ida), interval=1)
# %%
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
f = f'vid/genga_scatter_3.mp4' 
writervideo = FFMpegWriter(fps=30) # ffmpeg must be installed
im_ani.save(f, writer=writervideo)
# %%
data_all = [data_ida, data_cn]

fig, axes = plt.subplots(3, figsize=(8,12))
fig.subplots_adjust(hspace=0.4)

def animate(i, data_all):
    axes[0].cla()
    axes[0].set_xlabel("semi-major axis")
    axes[0].set_ylim(0,20)
    axes[0].hist([data_all[0][i,:,3], data_all[1][i,:,3]], bins=40, range=(0,7), label=['IDA', 'CN'])
    axes[0].set_title(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    axes[0].legend()
    
    axes[1].cla()
    axes[1].set_xlabel("eccentricity")
    axes[1].set_ylim(0,20)
    axes[1].hist([data_all[0][i,:,4], data_all[1][i,:,4]], bins=40, range=(0,0.06), label=['IDA', 'CN'])
    axes[1].legend()
    
    axes[2].cla()
    axes[2].set_xlabel("inclination")
    axes[2].set_ylim(0,20)
    axes[2].hist([data_all[0][i,:,5], data_all[1][i,:,5]], bins=40, range=(0,0.06), label=['IDA', 'CN'])
    axes[2].legend()

animation = FuncAnimation(fig, animate, frames=len(data_ida), fargs=(data_all,), interval=1)
# %%
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
writervideo = FFMpegWriter(fps=30) # ffmpeg must be installed
animation.save('vid/genga_hist_3.mp4', writer=writervideo)

# animation.save('vid/genga_hist_4.gif', writer='imagemagick', fps=60)
