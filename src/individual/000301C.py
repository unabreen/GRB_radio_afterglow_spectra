#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:00:46 2026

@author: unabreen


    
inferring column names from dataset, then set generic names for easier calculations
Cleaning:
    first clean for values with error above 3 standard deviations
    after separating into frequency bands, clean for nonnegative values and
    get rid of non numeric columns
    
Plotting: 
    plot light curves for general sense of brightness/timescale
    plot spectral index to observe for quality and get rid of outliers

Calculating spectral index;
    calculate and store in seperate dataframes, then download to csv files for
    future comparison
    
    
specific dataset information:
    originally in microJy
    
not enough data for use
    
    
"""

import pandas as pd

from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean



filename = '../old_grb_sample/000301C.dat'

redshift = 2.0335
lum_dist= 16144.5 *3.08568e24  # Mpc to pc


#inferred names
col_names = ['name',
             'telescope',
             'Month',
             '3',  #unknown, named as column with index 3
             'year',
             'Days',    
             'Frequency(GHz)',  # GHz 
             'Flux(mJy)',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             '9'  #unknown, named as column with index 9
             ]
# generic names of cols for operations
flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'

df = pd.read_csv(filename, sep='\s+', header=None, names = col_names)

upper_13=1.6*(1+redshift) 
lower_13=1.2*(1+redshift) 
upper_49= 5.0*(1+redshift) 
lower_49= 4.6*(1+redshift) 
upper_85 = 8.9*(1+redshift) 
lower_85 = 8.4*(1+redshift) 
upper_150 = 15*(1+redshift) 
lower_150 = 15*(1+redshift) 

uncertainty_clean(df, flux_col, flux_err_col)
clean_GRB_data(df, flux_col, flux_err_col, time_col)

# convert to mJy from microJy
df[flux_col] = df[flux_col] * 1e-3
df[flux_err_col] = df[flux_err_col] * 1e-3

freq_13, freq_49, freq_85, freq_150, df= radio_bands(
                                df,  
                                upper_13=1.6*(1+redshift), 
                                lower_13=1.2*(1+redshift), 
                                upper_49= 5.0*(1+redshift), 
                                lower_49= 4.6*(1+redshift), 
                                upper_85 = 8.9*(1+redshift), 
                                lower_85 = 8.4*(1+redshift), 
                                upper_150 = 15*(1+redshift), 
                                lower_150 = 15*(1+redshift), 
                                lum_dist= lum_dist,
                                redshift= redshift,
                                time_col= time_col,
                                freq_col= freq_col,
                                flux_col= flux_col,
                                flux_err_col = flux_err_col)

light_curves('0000301C', 
             Freq_85 = freq_85,
             Freq_49 = freq_49,
             Freq_13 = freq_13,
             time_col = time_col, 
             freq_col = freq_col, 
             flux_col= flux_col,
             flux_err_col= flux_err_col)

spix_49_85 = spectral_index(
    df, 
    lower = lower_49, 
    upper = upper_85, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col
    )

spix_13_49 = spectral_index(
    df, 
    lower = lower_13, 
    upper = upper_49, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col
)
spix_85_150 = spectral_index(df, 
    lower = lower_85, 
    upper = upper_150, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col)    

spectral_index_plots('000301C', 
                     freq_13_49 = spix_13_49,
                     freq_49_85 = spix_49_85,
                     freq_85_150= spix_85_150
                     )










