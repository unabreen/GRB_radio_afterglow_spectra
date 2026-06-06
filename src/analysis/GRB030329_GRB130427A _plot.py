#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 12:06:49 2026

@author: unabreen
"""
# working with radio data from GRB 130427A and GRB 030329
# make spectral index plot from 4.9 and 8.5 GHz
# flux scatter plot for 4.9 and 8.5 GHz

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import linregress

#first columnm observed freq in Hz
#second column is
col_names = ['Frequency(Hz)','Days','Flux(mJy)','Flux uncertainty']

GRB030329 = 'radioall.dat'
GRB130427A = 'broadbandalldata.dat'



dist_030 = 816.4  * 3.08568e22   # Mpc to m
dist_130 = 1811.7 * 3.08568e22   # Mpc to m


df_030 = pd.read_csv(GRB030329, sep='\s+', header = None, names = col_names)
df_130 = pd.read_csv(GRB130427A, sep='\s+', header = None, names = col_names)





freq_49_130 = df_130[(df_130['Frequency(Hz)'] <= 4.86e9) & (df_130['Frequency(Hz)'] >= 4.75e9)]  
#dataframe for GRB 130427A for frequency = 4.9GHz
freq_85_130 = df_130[(df_130['Frequency(Hz)'] >= 8.45e9) & (df_130['Frequency(Hz)'] <= 8.55e9)]
#dataframe for GRB 130427A for frequency = 8.5GHz
freq_14_130 = df_130[(df_130['Frequency(Hz)'] >= 1.4e9) & (df_130['Frequency(Hz)'] <= 1.45e9)]


freq_49_030 = df_030[(df_030['Frequency(Hz)'] <= 4.86e9) & (df_030['Frequency(Hz)'] >= 4.8e9)]  
#dataframe for GRB 030329 for frequency = 4.9GHz
freq_85_030 = df_030[df_030['Frequency(Hz)'] == 8.46e9]
#dataframe for GRB 030329 for frequency = 8.5GHz


#spectral index calc
radio_data_030 =  df_030[(df_030['Frequency(Hz)'] <= 8.55e9) & (df_030['Frequency(Hz)'] >= 4.8e9)] 
#radio_data_030 = radio_data_030.sort_values(by = 'Days')


#luminosity calculations
flux_total_030 = np.trapz(df_030['Flux(mJy)'], df_030['Frequency(Hz)'] )
luminosity_030 = flux_total_030 * 4 * np.pi * (dist_030)**2

flux_total_130 = np.trapz(df_130['Flux(mJy)'], df_130['Frequency(Hz)'] )
luminosity_130 = flux_total_130 * 4 * np.pi * (dist_130)**2
print(f"The total flux for 030 is {flux_total_030:.3e}, and for 130 {flux_total_130:.3e}")
print(f"The luminosity of 030 is {luminosity_030:.3e}, and for 130 {luminosity_130:.3e}")

'''
bins_len = [ ]
# first bin
for i in range(len(radio_data_030['Days'])):
   if radio_data_030['Frequency(Hz)'][i] != radio_data_030['Frequency(Hz)'][0]:
       bins_len.append([0, radio_data_030['Days'][i]])
    
    
for i in range(len(radio_data_030['Days'])):
  
    if radio_data_030['Frequency(Hz)'][i+1] != radio_data_030['Frequency(Hz)'][i]:
        bins_len.append([i+1, bins_len[-1]])
'''

# Flux plots
fig , ((ax1, ax2)) = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 6))

ax1.errorbar(freq_49_030['Days'], freq_49_030['Flux(mJy)'], 
             yerr= freq_49_030['Flux uncertainty'], fmt = 'o', c = 'r', label = '4.9 GHz')

ax1.errorbar(freq_85_030['Days'], freq_85_030['Flux(mJy)'], 
             freq_85_030['Flux uncertainty'], fmt = 'o',
          c ='b', label = '8.5 GHz')
ax1.set_title('GRB 030329 Flux')
ax1.set_xlabel('Days')
ax1.set_ylabel('Flux(mJy)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0, 10e3)
ax1.set_ylim(0, 10e2)



ax2.errorbar(freq_49_130['Days'], freq_49_130['Flux(mJy)'], 
             yerr= freq_49_130['Flux uncertainty'], fmt = 'o', c = 'r', label = '4.9 GHz')

ax2.errorbar(freq_85_130['Days'], freq_85_130['Flux(mJy)'], 
             freq_85_130['Flux uncertainty'], fmt = 'o',
          c ='b', label = '8.5 GHz')
ax2.errorbar(freq_14_130['Days'], freq_14_130['Flux(mJy)'], 
             freq_14_130['Flux uncertainty'], fmt = 'o',
          c ='g', label = '1.4 GHz')
ax2.set_title('GRB 130427A Flux')
ax2.set_xlabel('Days')
ax2.set_ylabel('Flux(mJy)')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0, 10e3)
ax2.set_ylim(0, 10e2)

ax1.legend(loc=2, fancybox=True, shadow=False, prop={"size": 12})
ax2.legend(loc=2, fancybox=True, shadow=False, prop={"size": 12})




#correcting times for redshift- 1 + z 
#

#luminosity - need luminosity distance - 






 
