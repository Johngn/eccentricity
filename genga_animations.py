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

data_ida = np.reshape(np.array(pd.read_csv('./data/dataida_mars_3.csv')), (-1, 101, 6))[::40]
data_cn = np.reshape(np.array(pd.read_csv('./data/datacn_mars_3.csv')), (-1, 101, 6))[::40]

# %%
times = data_ida[:,0,1]

fig, axes = plt.subplots(1, figsize=(8, 5))
axes.set_xlim(-0,6)
axes.set_ylim(0,0.15)
axes.set_xlabel("semi-major axis")
axes.set_ylabel("inclination")

planets_ida, = axes.plot([], [], "o", markersize=8, alpha=0.8, label="Ida")
planets_cn, = axes.plot([], [], "o", markersize=8, alpha=0.8, label="CN")
text = axes.text(2.2, 0.12, '', fontsize=15)

axes.legend()

def animate(i):
    planets_ida.set_data(data_ida[i,:,3], data_ida[i,:,5])
    planets_cn.set_data(data_cn[i,:,3], data_cn[i,:,5])
    text.set_text(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    return planets_ida, planets_cn
    
im_ani = FuncAnimation(fig, animate, frames=len(data_ida), interval=1)
# %%
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

f = f'vid/genga_mars_i_1.mp4' 
writervideo = FFMpegWriter(fps=30) # ffmpeg must be installed
im_ani.save(f, writer=writervideo)

# %%
data_ida = np.reshape(np.array(pd.read_csv('./data/dataida_mars_3.csv')), (-1, 101, 6))[::20]
data_cn = np.reshape(np.array(pd.read_csv('./data/datacn_mars_3.csv')), (-1, 101, 6))[::20]
times = data_ida[:,0,1]
data_all = [data_ida, data_cn]
# %%
fig, axes = plt.subplots(3, figsize=(8,12))
fig.subplots_adjust(hspace=0.4)

def animate(i, data_all):
    axes[0].cla()
    axes[0].set_xlabel("semi-major axis")
    axes[0].hist([data_all[0][i,:,3], data_all[1][i,:,3]], 20, label=['IDA', 'CN'])
    axes[0].set_title(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    axes[0].legend()
    
    axes[1].cla()
    axes[1].set_xlabel("eccentricity")
    axes[1].hist([data_all[0][i,:,4], data_all[1][i,:,4]], 20, label=['IDA', 'CN'])
    axes[1].legend()
    
    axes[2].cla()
    axes[2].set_xlabel("inclination")
    axes[2].hist([data_all[0][i,:,5], data_all[1][i,:,5]], 20, label=['IDA', 'CN'])
    axes[2].legend()

animation = FuncAnimation(fig, animate, frames=len(data_ida), fargs=(data_all,))
# %%
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

f = f'vid/genga_mars_hist_1.mp4' 
writervideo = FFMpegWriter(fps=30) # ffmpeg must be installed
animation.save(f, writer=writervideo)