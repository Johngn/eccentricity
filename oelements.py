#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:30:36 2020

@author: john
"""

import sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

files = glob(f'./data/Outcn_1body*.dat')
files.sort()
test = np.loadtxt(files[0])
data = [np.loadtxt(i) for i in files]
data = np.concatenate([np.reshape(item, (1, len(item))) for i, item in enumerate(data)])


k = 0.0172020989
G = k**2
mu = G

t = [item[0] for i, item in enumerate(data)]