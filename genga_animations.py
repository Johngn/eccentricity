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
data_ida = np.reshape(np.array(pd.read_csv('./data/dataida_mars_3.csv')), (-1, 101, 6))[::3000]
data_ida = np.reshape(np.array(pd.read_csv('./data/datacn_mars_3.csv')), (-1, 101, 6))[::3000]
# data_ida = np.array(pd.read_csv('./data/datacn_1body.csv'))
# data_cn = np.array(pd.read_csv('./data/datacn_1body.csv'))
times = data_ida[:,0,1]
data_all = [data_ida, data_cn]
# %%
fig, axes = plt.subplots(3, figsize=(8,12))
fig.subplots_adjust(hspace=0.4)

def animate(i, data_all):
    axes[0].cla()
    axes[0].set_xlabel("semi-major axis")
    axes[0].set_ylim(0,20)
    axes[0].hist([data_all[0][i,:,3], data_all[1][i,:,3]], bins=40, range=(0,7), label=['IDA', 'CN'])
    # axes[0].set_title(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
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

animation = FuncAnimation(fig, animate, frames=len(data_ida), fargs=(data_all,))
# %%
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

f = f'vid/genga_mars_hist_1.mp4' 
writervideo = FFMpegWriter(fps=30) # ffmpeg must be installed
animation.save(f, writer=writervideo)