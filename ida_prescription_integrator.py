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
dt = year*1         # time step

e_results1 = np.zeros(noutputs)
e_results2 = np.zeros(noutputs)
e_results3 = np.zeros(noutputs)
e_results4 = np.zeros(noutputs)

e = 0.3  
for i in range(noutputs):    
    ehat = e/h
    t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    de = -e/t_ecc*dt
    e = e+de
    e_results1[i] = e
    
e = 0.25
for i in range(noutputs):    
    ehat = e/h    
    # t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
    t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    de = -e/t_ecc*dt
    e = e+de
    e_results2[i] = e
    
e = 0.15
for i in range(noutputs):    
    ehat = e/h    
    # t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
    t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    de = -e/t_ecc*dt
    e = e+de
    e_results3[i] = e

e = 0.05
for i in range(noutputs):    
    ehat = e/h    
    # t_ecc = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3)   # equation 11 from Creswell+Nelson 08
    t_ecc = t_wave/0.78*(1+(1/15)*ehat**3)              # equation 34 from Ida 20
    de = -e/t_ecc*dt
    e = e+de
    e_results4[i] = e
    
fig, ax = plt.subplots(1, figsize=(7,5))

ax.plot(np.arange(0,noutputs,1)*dt/year, e_results1)
ax.plot(np.arange(0,noutputs,1)*dt/year, e_results2)
ax.plot(np.arange(0,noutputs,1)*dt/year, e_results3)
ax.plot(np.arange(0,noutputs,1)*dt/year, e_results4)
ax.set_xlabel('time (years)')
ax.set_ylabel('eccentricity')
ax.set_xlim(0, 3000)
ax.set_ylim(0)
ax.tick_params(which='both', direction="in", top=True, right=True)

# fig.savefig('/home/john/Desktop/summerproject/img/eccentricity_over_time.png', bbox_inches='tight')
