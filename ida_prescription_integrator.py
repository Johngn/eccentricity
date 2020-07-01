# %%
import numpy as np
import matplotlib.pyplot as plt

G = 6.674e-11                                           # gravitational constanst
au = 1.496e11                                           # astronomical unit
a = au*1                                                # radial position of planet
mstar = 1.989e30                                         # mass of star
year = 365.25*24.*60.*60.                               # year
mplanet = 5.972e25                                        # mass of orbiting planet
omegak = np.sqrt(G*(mstar)/a**3)                # Keplerian frequency
sigma = 17000*(a/au)**(-3/2)                              # surface density of disc
h = 0.05                                                # disc scale height H/r
t_wave = mstar/mplanet*mstar/sigma/a**2*h**4/omegak   # equation 7 from Ida 2020

noutputs = 100000
dt = year*0.01        # time step

e_results1 = np.zeros(noutputs)
e_results2 = np.zeros(noutputs)
e_results3 = np.zeros(noutputs)
e_results4 = np.zeros(noutputs)

e = 0.15

for i in range(noutputs):    
    ehat = e/h
    # t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
    de = -e/t_ecc*dt
    e = e+de
    e_results1[i] = e
    
# e = 0.25
# for i in range(noutputs):    
#     ehat = e/h    
#     # t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
#     t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
#     de = -e/t_ecc*dt
#     e = e+de
#     e_results2[i] = e
    
# e = 0.15
# for i in range(noutputs):    
#     ehat = e/h    
#     # t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
#     t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
#     de = -e/t_ecc*dt
#     e = e+de
#     e_results3[i] = e

# e = 0.05
# for i in range(noutputs):    
#     ehat = e/h    
#     # t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
#     t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
#     de = -e/t_ecc*dt
#     e = e+de
#     e_results4[i] = e

fig, ax = plt.subplots(1, figsize=(7,5))

ax.plot(np.arange(0,noutputs,1)*dt/year, e_results1, label='analytical')
ax.plot(np.arange(0,noutputs,1)*dt/year, e_total, label='numerical')
# ax.plot(np.arange(0,noutputs,1)*dt/year, e_results2)
# ax.plot(np.arange(0,noutputs,1)*dt/year, e_results3)
# ax.plot(np.arange(0,noutputs,1)*dt/year, e_results4)
ax.axhline(h, c='black', label='e = H/r')
ax.set_xlabel('time (years)')
ax.set_ylabel('eccentricity')
ax.set_xlim(0, noutputs*0.01)
# ax.set_ylim(0, 0.1)
ax.tick_params(which='both', direction="in", top=True, right=True)
ax.grid()
ax.legend()

# fig.savefig('/home/john/Desktop/summerproject/img/higheccentricity_num+analyticCN.png', bbox_inches='tight')
