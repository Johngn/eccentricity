# %%
import numpy as np
import matplotlib.pyplot as plt

G = 6.674e-11                                           # gravitational constanst
au = 1.496e11                                           # astronomical unit
r = 1*au                                                # radial position of planet
mstar = 2.0e30                                          # mass of star
year = 365.25*24.*60.*60.                               # year
mplanet = 6.0e25                                        # mass of orbiting planet
omegak = np.sqrt(G*(mstar+mplanet)/r**3)                # Keplerian frequency
sigma = 17000*(r/au)**(-1)                              # surface density of disc
h = 0.05                                                # disc scale height H/r
t_wave = mstar/mplanet*mstar/sigma/r**2*h**(4)/omegak   # equation 7 from Ida 2020

noutputs = 3000
e = 0.3             # initial eccentricity
dt = year*1         # time step

e_results = np.zeros(noutputs)

for i in range(noutputs):    
    ehat = e/h
    t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    # t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
    de = -e/t_ecc*dt
    e = e+de
    e_results[i] = e
    
fig, ax = plt.subplots(1, figsize=(8,6))

ax.plot(np.arange(0,noutputs,1)*dt/year, e_results)
ax.set_xlabel('time')
ax.set_ylabel('eccentricity')
ax.set_xlim(0)
ax.set_ylim(0)
ax.tick_params(which='both', direction="in", top=True, right=True)

