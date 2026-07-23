#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 15:43:19 2026

@author: unabreen

Examining the evolution of the spectrum over time. Data broken into early 
and late times to more clearly see. 

"""
import pandas as pd
import numpy as np
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import luminosity_func
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean
import matplotlib.pyplot as plt


col_names = ['name',
             'telescope',
             'Month',
             'Column 3',  # unknown, named as column with index 3
             'year',
             'Days',    
             'Frequency(GHz)',  # GHz 
             'Flux(mJy)',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             'Column 9'   # unknown, named as column with index 9
             ]


flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'


filename = '../old_grb_sample/980425.dat'


redshift= 0.0087   # approximation
lum_dist =  37.7 * 3.08568*10e26  # Mpc-> cm


df = pd.read_csv(filename, sep='\s+', header=None, names = col_names)

#keep copy of non cleaned dataframe for reference
df_old = df

# convert microJy to milliJy- 1e-6 to 1e-3
df[flux_col] = df[flux_col] * 1e-3
df[flux_err_col] = df[flux_err_col] * 1e-3

uncertainty_clean(df, flux_col, flux_err_col)
clean_GRB_data(df, flux_col, flux_err_col, time_col)
             
# adjusting for redshift      
upper_13=1.6*(1+redshift)
lower_13=1.2*(1+redshift) 
upper_49= 5.0*(1+redshift) 
lower_49= 4.6*(1+redshift) 
upper_85 = 8.9*(1+redshift) 
lower_85 = 8.4*(1+redshift) 
upper_150 = 15*(1+redshift) 
lower_150 = 15*(1+redshift) 
          

freq_13, freq_49, freq_85, freq_150, df= radio_bands(
                                df,  
                                upper_13, 
                                lower_13, 
                                upper_49, 
                                lower_49, 
                                upper_85, 
                                lower_85, 
                                upper_150, 
                                lower_150, 
                                lum_dist= lum_dist,
                                redshift= redshift,
                                time_col= time_col,
                                freq_col= freq_col,
                                flux_col= flux_col,
                                flux_err_col = flux_err_col
                                )
  

    
light_curves('980425', 
             Freq_85 = freq_85,
             Freq_49 = freq_49,
             Freq_13 = freq_13,
             time_col = time_col, 
             freq_col = freq_col, 
             flux_col= flux_col,
             flux_err_col= flux_err_col)



#early time spectrum
fig, ((ax1, ax2)) = plt.subplots(nrows = 1, ncols = 2, figsize= (10, 5))


early_df = df[df[time_col] < 16]
late_df = df[df[time_col] > 70]

ax1.scatter(early_df[freq_col], early_df[flux_col])
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title('Spectrum early times')

ax2.scatter(late_df[freq_col], late_df[flux_col])
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title('Spectrum late times')