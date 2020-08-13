import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

data_ida = np.array(pd.read_csv('./data/dataida_mars_4.csv'))
data_cn = np.array(pd.read_csv('./data/datacn_mars_4.csv'))
# %%
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
    
data_ida_new = np.reshape(data_ida_new, (-1, 101, 6))[::100]
data_cn_new = np.reshape(data_cn_new, (-1, 101, 6))[::100]
# %%
times = data_ida_new[:,0,1]

fig, axes = plt.subplots(2, figsize=(8, 8))
axes[0].set_xlim(-0,6)
axes[0].set_ylim(0,0.15)
axes[0].set_ylabel("eccentricity")

axes[1].set_xlim(-0,6)
axes[1].set_ylim(0,0.15)
axes[1].set_xlabel("semi-major axis")
axes[1].set_ylabel("inclination")

planets_ida_e, = axes[0].plot([], [], "o", markersize=8, alpha=0.8, label="Ida")
planets_cn_e, = axes[0].plot([], [], "o", markersize=8, alpha=0.8, label="CN")
planets_ida_i, = axes[1].plot([], [], "o", markersize=8, alpha=0.8, label="Ida")
planets_cn_i, = axes[1].plot([], [], "o", markersize=8, alpha=0.8, label="CN")
text = axes[0].text(2.2, 0.12, '', fontsize=15)

axes[0].legend()
axes[1].legend()

def animate(i):
    planets_ida_e.set_data(data_ida_new[i,:,3], data_ida_new[i,:,4])
    planets_cn_e.set_data(data_cn_new[i,:,3], data_cn_new[i,:,4])
    planets_ida_i.set_data(data_ida_new[i,:,3], data_ida_new[i,:,5])
    planets_cn_i.set_data(data_cn_new[i,:,3], data_cn_new[i,:,5])
    text.set_text(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    return planets_ida_e, planets_cn_e, planets_ida_i, planets_cn_i
    
im_ani = FuncAnimation(fig, animate, frames=len(data_ida_new), interval=1)
# %%
solar_to_earth_mass = 1.99e30/5.97e24

fig, axes = plt.subplots(2, figsize=(8,8))
# fig.subplots_adjust(hspace=0.4)
# i = 10
# markers = data_ida_new[-1,:,2]*solar_to_earth_mass*200
# plt.scatter(data_ida_new[0,:,3], data_ida_new[0,:,4], s=markers)

def animate(i, data_ida_new):
    axes[0].cla()
    axes[0].set_title(r'{} $\times 10^3$ years'.format(int(times[i]/1e3)))
    axes[0].set_xlim(-0,6)
    axes[0].set_ylim(0,0.15)
    markers = data_ida_new[0,:,2]*solar_to_earth_mass*200
    # print(markers[0])
    
    axes[0].scatter(data_ida_new[i,:,3], data_ida_new[i,:,4], s=markers, alpha=0.8, label="Ida")
    
    # axes[1].cla()
    # axes[1].plot(data_cn_new[i,:,3], data_cn_new[i,:,4], "o", markersize=8, alpha=0.8, label="CN")

animation = FuncAnimation(fig, animate, frames=len(data_ida_new), fargs=(data_ida_new,))

# %%
x = [1,2,3,4,5]
y = [a**2 for a in x]
s = np.array([10*4**n for n in range(len(x))])
plt.scatter(x,y,s=s)
plt.title('Doubling width of marker in scatter plot')
plt.xlabel('x')
plt.ylabel('x**2')
plt.xlim(0,6)
plt.ylim(0,30)