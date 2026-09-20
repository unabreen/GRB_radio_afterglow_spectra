#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 17:43:22 2026

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
    - no data in 15 GHz
    


"""
import sys
from pathlib import Path

# Points to .../src/analysis, relative to this script's own location
sys.path.append(str(Path(__file__).resolve().parent.parent / "analysis"))


import pandas as pd
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import luminosity_func
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean



col_names = ['name',
             'telescope',
             'Month',
             '3',  # unknown, named as column with index 3
             'year',
             'Days',    
             'Frequency(GHz)',  # GHz 
             'Flux(mJy)',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             '9'   # unknown, named as column with index 9
             ]


flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'



''
filename = '../../data/raw/980425.dat'


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



spectral_index_plots('980425', 
                     freq_13_49 = spix_13_49,
                     freq_49_85 = spix_49_85
                     )
# make sure folder exists before downloading
# lines commented out when running and redownload unwanted
#spix_13_49.to_csv('spectral data files/980425 spectra/980425_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('spectral data files/980425 spectra/980425_spectral_index(4.9-8.5).csv', index=True)

luminosity_func(df, lum_dist, redshift, flux_col)



















