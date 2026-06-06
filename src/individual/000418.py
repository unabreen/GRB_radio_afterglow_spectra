#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 12:41:25 2026

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
    no 1.3 GHz data
    originally in microJy, GHz
"""


import pandas as pd

from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean


col_names = ['name',
             'telescope',
             'Month',
             '3', # unknown, named as column with index 3
             'year',
             'Days',    
             'Frequency(GHz)',  # GHz 
             'Flux(mJy)',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             '9' # unknown, named as column with index 3
             ]


flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'


filename = '../old_grb_sample/000418.dat'


redshift= 1.119
lum_dist =  7705.9 * 3.08568*10e26  # Mpc  -> cm



df = pd.read_csv(filename, sep='\s+', header=None, names = col_names)

# original dataframe for reference
df_old = df

# band limits
upper_13 = 1.41*(1 + redshift)
lower_13 = 1.25*(1 + redshift)
upper_49 = 5.1*(1 + redshift)
lower_49 = 4.75*(1 + redshift)
upper_85 = 8.7*(1 + redshift)
lower_85 = 8.4*(1 + redshift) 
upper_150 = 16*(1 + redshift)
lower_150 = 14.5*(1 + redshift)


# convert microJy to milliJy- 1e-6 to 1e-3
df[flux_col] = df[flux_col] * 1e-3
df[flux_err_col] = df[flux_err_col] * 1e-3

uncertainty_clean(df, flux_col, flux_err_col)
clean_GRB_data(df, flux_col, flux_err_col, time_col)

df_post_clean = df

freq_13, freq_49, freq_85, freq_150, df= radio_bands(df, 
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
                                                     freq_col = freq_col, 
                                                     flux_col= flux_col,
                                                     flux_err_col= flux_err_col)


light_curves('000418', 
            
             Freq_85 = freq_85,
             Freq_49 = freq_49,
             Freq_150= freq_150,
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

spix_85_150 = spectral_index(
    df, 
    lower = lower_85, 
    upper = upper_150, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col
)



spectral_index_plots('000418', 
                     freq_85_150 = spix_85_150,
                     freq_49_85 = spix_49_85
                     )

# make sure folder exists before downloading
# lines commented out when running and redownload unwanted
#spix_13_49.to_csv('spectral data files/980425 spectra/980425_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('spectral data files/000418 spectra/000418_spectral_index(4.9-8.5).csv', index=True)










