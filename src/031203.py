#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 22:16:00 2026

@author: unabreen


All VLA observations


not many points in each band

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import linregress
from grb_functions import radio_bands, spectral_index, light_curves, spectral_index_plots

col_names = ['name',
             'telescope',
             'year',
             'Month',
             '4',
             'Day',    #maybe day
             'Frequency(GHz)',  # GHz for VLA at least
             'Flux',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             '9'
             ]

flux_col = 'Flux'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'

filename = '../old_grb_sample/031203.dat'


redshift= 0.1055   

df_031203 = pd.read_csv(filename, 
                        sep='\s+', 
                        header=None, 
                        names = col_names
                        )


freq_85 = (df_031203[df_031203.iloc[:, 6] == 8.46 ])
#print(freq_85)
plt.scatter(freq_85.iloc[:, 5], freq_85.iloc[:, 7])
plt.title("GRB 031203 possible 8.5 light curve")
plt.xlabel("days")
plt.ylabel('Flux(unknown units')
plt.xscale('log')
plt.yscale('log')







