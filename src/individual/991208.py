#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:26:34 2026

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
    - originally in microJanskys
    - same column format as ~1990s GRBs- used 980425
    - large uncertanties for 1.3 GHz band
    - 4.9-8.5 and 8.5-15 bands used
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
             'col_3(placeholder)',
             'year',
             'Days',
             'Frequency(GHz)',
             'Flux(mJy)',
             'Flux uncertainty',
             'col_9(placeholder)'
             ]

flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'

filename = '../old_grb_sample/991208.dat'

redshift = 0.706
lum_dist = 4357.2 * 3.08568*10e26  # Mpc-> cm

df = pd.read_csv(filename, sep='\s+', header=None, names=col_names)

df[flux_col] = df[flux_col] * 1e-3
df[flux_err_col] = df[flux_err_col] * 1e-3

uncertainty_clean(df, flux_col, flux_err_col)
clean_GRB_data(df, flux_col, flux_err_col, time_col)

upper_13 = 1.5*(1 + redshift)
lower_13 = 1.25*(1 + redshift)
upper_49 = 5.1*(1 + redshift)
lower_49 = 4.75*(1 + redshift)
upper_85 = 8.7*(1 + redshift)
lower_85 = 8.35*(1 + redshift)
upper_150 = 16*(1 + redshift)
lower_150 = 14.5*(1 + redshift)

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

light_curves('991208', 
             Freq_150=freq_150,
             Freq_85=freq_85,
             Freq_49=freq_49,
             Freq_13=freq_13,
             time_col=time_col, 
             freq_col=freq_col, 
             flux_col=flux_col,
             flux_err_col=flux_err_col)

spix_85_150 = spectral_index(
    df, 
    lower=lower_85, 
    upper=upper_150, 
    freq_col=freq_col,
    time_col=time_col,
    flux_col=flux_col
    )
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

spix_85_150 = spix_85_150[(spix_85_150['alpha'] > -3) & (spix_85_150['alpha'] < 4)]

spectral_index_plots('991208', 
                     freq_85_150=spix_85_150,
                     freq_13_49=spix_13_49,
                     freq_49_85=spix_49_85
                     )

#spix_49_85.to_csv('spectral data files/991208 spectra/991208_spectral_index(4.9-8.5).csv', index=True)
#spix_85_150.to_csv('spectral data files/991208 spectra/991208_spectral_index(8.5-15).csv', index=True)



