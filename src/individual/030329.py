#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 15:46:07 2026

@author: unabreen
"""
# working with radio data from GRB 030329
# make spectral index plot from 4.9 and 8.5 GHz
# flux scatter plot for 4.9 and 8.5 GHz

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import linregress
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import luminosity_func
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean

"""
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
    
"""
#first columnm observed freq in Hz- changed to GHz

col_names = ['Frequency(GHz)','Days','Flux(mJy)','Flux uncertainty']

flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'

filename = '../old_grb_sample/radio030329.dat'   

redshift = 0.1685


# radio specific band limits
upper_13 = 1.41*(1 + redshift)
lower_13 = 1.25*(1 + redshift)
upper_49 = 5.1*(1 + redshift)
lower_49 = 4.75*(1 + redshift)
upper_85 = 8.7*(1 + redshift)
lower_85 = 8.4*(1 + redshift) 
upper_150 = 16*(1 + redshift)
lower_150 = 14.5*(1 + redshift)

lum_dist = 816.4  * 3.08568e24  # Mpc to cm

df = pd.read_csv(filename, sep='\s+', header = None, names = col_names)


#change Hz to GHz
df[freq_col] = df[freq_col]* 1e-9

df = uncertainty_clean(df, flux_col, flux_err_col)



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


light_curves('030329', 
            
             Freq_85 = freq_85,
             Freq_49 = freq_49,
             Freq_13 = freq_13,
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


spix_13_49 = spix_13_49[(spix_13_49['alpha'] > -3) & (spix_13_49['alpha'] < 4)]


spectral_index_plots('030329', 
                     freq_13_49 = spix_13_49,
                     freq_49_85 = spix_49_85,
                     freq_85_150= spix_85_150
                     )


# make sure there is a folder for {GRB} spectra before saving
#spix_13_49.to_csv('spectral data files/030329 spectra/030329_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('spectral data files/030329 spectra/030329_spectral_index(4.9-8.5).csv', index=True)
#spix_85_150.to_csv('spectral data files/030329 spectra/030329_spectral_index(8.5-15).csv', index=True)






