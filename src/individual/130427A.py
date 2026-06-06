#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:51:51 2026

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

specific dataset information
- frequency originally in Hz
- 1.3-4.9 and 4.9-8.5 bands used
"""


import pandas as pd
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean


col_names = ['Frequency(GHz)', #originally in Hz
             'Days',
             'Flux(mJy)',
             'Flux uncertainty']



flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'


filename = '../old_grb_sample/broadbandalldata.dat'


redshift = 0.34
lum_dist = 1811.7  * 3.08568*10e26   # Mpc to cm


df = pd.read_csv(filename, sep='\s+', header=None, names = col_names)

#df[freq_col] = df[freq_col] / 1e9

#keep copy of non cleaned dataframe for reference
df_old = df



#uncertainty_clean(df, flux_col, flux_err_col)
#clean_GRB_data(df, flux_col, flux_err_col, time_col)
             
# adjusting for redshift  
"""    
upper_13 = 1.41*(1 + redshift)
lower_13 = 1.25*(1 + redshift)
upper_49 = 5.2*(1 + redshift)
lower_49 = 4.75*(1 + redshift)
upper_85 = 8.6*(1 + redshift)
lower_85 = 8.4*(1 + redshift)
upper_150 = 16.1*(1 + redshift)
lower_150 = 13.9*(1 + redshift)
"""

upper_13 = 1.41e9*(1 + redshift)
lower_13 = 1.25e9*(1 + redshift)
upper_49 = 5.2e9*(1 + redshift)
lower_49 = 4.75e9*(1 + redshift)
upper_85 = 8.6e9*(1 + redshift)
lower_85 = 8.4e9*(1 + redshift)
upper_150 = 1.61e10*(1 + redshift)
lower_150 = 1.39e10*(1 + redshift)
          

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
  

    
light_curves('130427A', 
             Freq_85 = freq_85,
             Freq_49 = freq_49,
             Freq_13 = freq_13,
             Freq_150= freq_150,
             time_col = time_col, 
             freq_col = freq_col, 
             flux_col= flux_col,
             flux_err_col= flux_err_col)


spix_85_150 = spectral_index(
    df, 
    lower = lower_85, 
    upper = upper_150, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col
    )

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

spix_13_49= spix_13_49[(spix_13_49['alpha'] <= 5) ]
spix_13_49 = spix_13_49[(spix_13_49['alpha'] >= -3)]  

spectral_index_plots('130427A', 
                     freq_13_49 = spix_13_49,
                     freq_49_85 = spix_49_85,
                     freq_85_150 = spix_85_150
                     )
# make sure folder exists before downloading
# lines commented out when running and redownload unwanted
#spix_13_49.to_csv('spectral data files/130427A spectra/130427A_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('spectral data files/130427A spectra/130427A_spectral_index(4.9-8.5).csv', index=True)





