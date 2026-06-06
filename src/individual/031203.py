#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 22:16:00 2026

@author: unabreen


All VLA observations


not many points in each band

Workflow:
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
             'year',
             'Month',
             '4', # unknown, named as column with index 4
             'Days',    #maybe day
             'Frequency(GHz)',  # GHz for VLA at least
             'Flux',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             '9' # unknown, named as column with index 4
             ]


flux_col = 'Flux'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'

filename = '../old_grb_sample/031203.dat'


redshift= 0.1055   
lum_dist= 490.7  * 3.08568e24  # Mpc to cm 

df = pd.read_csv(filename, 
                        sep='\s+', 
                        header=None, 
                        names = col_names
                        )



uncertainty_clean(df, flux_col, flux_err_col)
clean_GRB_data(df, flux_col, flux_err_col, time_col)
                   
                   

freq_13, freq_49, freq_85, freq_150, df_red = radio_bands(
                                df,  
                                upper_13=1.43*(1+redshift), 
                                lower_13=1.43*(1+redshift), 
                                upper_49= 5.0*(1+redshift), 
                                lower_49= 4.6*(1+redshift), 
                                upper_85 = 8.6*(1+redshift), 
                                lower_85 = 8.4*(1+redshift), 
                                upper_150 = 15*(1+redshift), 
                                lower_150 = 15*(1+redshift), 
                                lum_dist= lum_dist,
                                redshift= redshift,
                                time_col= time_col,
                                freq_col= freq_col,
                                flux_col= flux_col,
                                flux_err_col = flux_err_col
                                )
freq_bands = [freq_13, freq_49, freq_85, freq_150]

# must clean after making frequency band frames
freq_49 = clean_GRB_data(freq_49, flux_col, flux_err_col, time_col)
freq_13 = clean_GRB_data(freq_13, flux_col, flux_err_col, time_col)
freq_85 = clean_GRB_data(freq_85, flux_col, flux_err_col, time_col)
freq_150 = clean_GRB_data(freq_150, flux_col, flux_err_col, time_col)
  
    
light_curves('031203', 
             
             Freq_85 = freq_85,
             Freq_49 = freq_49,
             Freq_13 = freq_13,
             time_col = time_col, 
             freq_col = freq_col, 
             flux_col= flux_col,
             flux_err_col= flux_err_col)


upper_13=1.43*(1+redshift) 
lower_13=1.43*(1+redshift) 
upper_49= 5.0*(1+redshift) 
lower_49= 4.6*(1+redshift) 
upper_85 = 8.6*(1+redshift) 
lower_85 = 8.4*(1+redshift) 
upper_150 = 15*(1+redshift) 
lower_150 = 15*(1+redshift)


spix_49_85 = spectral_index(
    df_red, 
    lower = lower_49, 
    upper = upper_85, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col
    )

spix_13_49 = spectral_index(
    df_red, 
    lower = lower_13, 
    upper = upper_49, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col
    )

spix_49_85 = spix_49_85[spix_49_85['alpha'] > -3]

spectral_index_plots('031203',
                     freq_49_85 = spix_49_85, 
                     freq_13_49 = spix_13_49)

#spix_13_49.to_csv('031203_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('031203_spectral_index(4.9-8.5).csv', index=True)










