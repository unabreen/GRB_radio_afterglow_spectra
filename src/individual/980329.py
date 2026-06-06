#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 21:42:31 2026

@author: unabreen

observations from VLA(GHz?), JCMT(sub mm?), OVRO(MHz?)
- must convert all frequencies into Hz
"""

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


col_names = ['name',
             'telescope',
             'Month',
             '3',
             'year',
             'Days',    #maybe day
             'Frequency(GHz)',  # GHz 
             'Flux(erg/s cm^2)',   # in microJy
             'Flux uncertainty',  #possibly
             '9'
             ]


flux_col = 'Flux(erg/s cm^2)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'


filename = '../old_grb_sample/980329.dat'


redshift= 5   # approximation
lum_dist = 47647.9 * 3.08568*10e26  # Mpc -> pc -> cm



df = pd.read_csv(filename, sep='\s+', header=None, names = col_names)

df[flux_col] = df[flux_col] * 10**(-29)
df[flux_col] = df[flux_col] * 10**(-29)
df[flux_err_col] = df[flux_err_col] * 10**(-29)
df[flux_err_col] = df[flux_err_col] * 10**(-29)


df= uncertainty_clean(df, flux_col, flux_err_col)
                   
#print(sigma_three)


freq_13, freq_49, freq_85, freq_150, df = radio_bands(
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
df =  clean_GRB_data(df, flux_col, flux_err_col, time_col)
  
  
light_curves('980329', 
             Freq_150 = freq_150,
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
spix_85_150 = spectral_index(
    df, 
    lower = lower_85, 
    upper = upper_150, 
    freq_col= freq_col,
    time_col= time_col,
    flux_col= flux_col
    )



spectral_index_plots('980329', spix_85_150, spix_49_85, spix_13_49)






