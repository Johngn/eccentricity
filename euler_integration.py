#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:48:24 2020

@author: johngillan
"""
import numpy as np

def euler(noutputs, e, i, h, dt, t_wave, ida):   
    e_euler = np.zeros(noutputs)
    i_euler = np.zeros(noutputs)
    
    for j in range(noutputs):    
        ehat = e/h    
        ihat = i/h
        if ida:
            t_e = t_wave/0.78*(1+1/15*(ehat**2+ihat**2)**(3/2))                   # equation 34 from Ida 20
            t_i = t_wave/0.544*(1+1/21.5*(ehat**2+ihat**2)**(3/2))                # Ida
        else:
            t_e = t_wave/0.78*(1-0.14*ehat**2+0.06*ehat**3+0.18*ehat*ihat**2)     # equation 11 from Creswell+Nelson 08
            t_i = t_wave/0.544*(1-0.30*ihat**2+0.24*ihat**3+0.14*ihat*ehat**2)    # CN    
        de = -e/t_e*dt
        e = e+de
        e_euler[j] = e    
        di = -i/t_i*dt
        i = i+di
        i_euler[j] = i
        
    return [e_euler, i_euler]