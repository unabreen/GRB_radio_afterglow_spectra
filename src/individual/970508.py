

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 15:46:07 2026

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
    
    
Specific notes about dataset:
    flux originally in microJanskys- must convert to milliJanskys
    frequency in Hz
    no data in 15 GHz
"""

import pandas as pd
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean

col_names = ['Frequency(GHz)', 'Days', 'Flux(mJy)', 'Flux uncertainty']

flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'

GRB970508 = '../old_grb_sample/radio970508.dat'

redshift = 0.835
lum_dist = 5358.8 * 3.08568e24   # Mpc to m

# dataset specific limits on radio bands
upper_13 = 1.41*(1 + redshift)
lower_13 = 1.25*(1 + redshift)
upper_49 = 5.1*(1 + redshift)
lower_49 = 4.75*(1 + redshift) 
upper_85 = 8.6*(1 + redshift)
lower_85 = 8.4*(1 + redshift)
upper_150 = 16*(1 + redshift)
lower_150 = 14.5*(1 + redshift)

df = pd.read_csv(GRB970508, sep='\s+', header=None, names=col_names)

df[freq_col] = df[freq_col] * 1e-9
df[flux_col] = df[flux_col] * 1e-3
df[flux_err_col] = df[flux_err_col] * 1e-3

uncertainty_clean(df, flux_col, flux_err_col)
clean_GRB_data(df, flux_col, flux_err_col, time_col)

freq_13, freq_49, freq_85, freq_150, df = radio_bands(df, 
                                                     upper_13, 
                                                     lower_13,
                                                     upper_49,
                                                     lower_49, 
                                                     upper_85, 
                                                     lower_85, 
                                                     upper_150, 
                                                     lower_150, 
                                                     lum_dist=lum_dist, 
                                                     redshift=redshift, 
                                                     time_col=time_col,
                                                     freq_col=freq_col, 
                                                     flux_col=flux_col,
                                                     flux_err_col=flux_err_col)

light_curves('970508', 
             Freq_85=freq_85,
             Freq_49=freq_49,
             Freq_13=freq_13,
             time_col=time_col, 
             freq_col=freq_col, 
             flux_col=flux_col,
             flux_err_col=flux_err_col)

spix_49_85 = spectral_index(
    df, 
    lower=lower_49, 
    upper=upper_85, 
    freq_col=freq_col,
    time_col=time_col,
    flux_col=flux_col
    )

spix_13_49 = spectral_index(
    df, 
    lower=lower_13, 
    upper=upper_49, 
    freq_col=freq_col,
    time_col=time_col,
    flux_col=flux_col
    )

spix_13_49 = spix_13_49[(spix_13_49['alpha'] < 4)]

spectral_index_plots('970508',
                     freq_49_85=spix_49_85,
                     freq_13_49=spix_13_49)

#spix_13_49.to_csv('spectral data files/970508 spectra/970508_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('spectral data files/970508 spectra/970508_spectral_index(4.9-8.5).csv', index=True)










