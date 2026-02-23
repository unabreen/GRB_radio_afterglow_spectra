#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 21:42:31 2026

@author: unabreen

observations from VLA(GHz?), JCMT(sub mm?), OVRO(MHz?)
- must convert all frequencies into Hz
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import linregress
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import luminosity_func


col_names = ['name',
             'telescope',
             'Month',
             '3',
             'year',
             'Days',    #maybe day
             'Frequency(GHz)',  # GHz for VLA at least
             'Flux',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             '9'
             ]





filename = '../old_grb_sample/980329.dat'


redshift= 5   # approximation
lum_dist = 47647.9 * 3.08568*10e30  # Mpc -> pc -> cm

df = pd.read_csv(filename, sep='\s+', header=None, names = col_names)



# Drop rows where any object (string) column contains 'some_value'
#obj_cols = df.select_dtypes(include=['object']).columns
#df = df[df[obj_cols].eq('VLA').any(axis=1)]

'''
freq_85 = (df[df.iloc[:, 6] == 8.46])
#print(freq_85)
plt.scatter(freq_85.iloc[:, 5], freq_85.iloc[:, 7])
plt.title("GRB 980329 possible 8.5 light curve")
plt.xlabel('Days')
plt.ylabel('Flux-unknown units')
plt.xscale('log')
plt.yscale('log')

'''
radio_bands(df, 
            redshift, 
            upper_13=1.5, 
            lower_13=1.2, 
            upper_49= 4.5, 
            lower_49= 5.0, 
            upper_85 = 8.6, 
            lower_85 = 8.4, 
            upper_150 = 15, 
            lower_150 = 15, 
            lum_dist = lum_dist)






