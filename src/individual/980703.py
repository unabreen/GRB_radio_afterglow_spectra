#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 12:32:02 2026
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
    - originally in milliJanskys
    - no data in 15 GHz
    - use 4.9-8.5 band
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
             'Flux(unknown units)',
             'Flux uncertainty',
             'col_9(placeholder)'
             ]

flux_col = 'Flux(unknown units)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'

filename = '../old_grb_sample/980703.dat'

redshift = 0.966   # approximation
lum_dist = 6419.4 * 3.08568*10e26  # Mpc -> pc -> cm

df = pd.read_csv(filename, sep='\s+', header=None, names=col_names)

uncertainty_clean(df, flux_col, flux_err_col)

freq_13, freq_49, freq_85, freq_150, df = radio_bands(
                                df,  
                                upper_13=1.6*(1+redshift), 
                                lower_13=1.2*(1+redshift), 
                                upper_49=5.0*(1+redshift), 
                                lower_49=4.6*(1+redshift), 
                                upper_85=8.9*(1+redshift), 
                                lower_85=8.4*(1+redshift), 
                                upper_150=15*(1+redshift), 
                                lower_150=15*(1+redshift), 
                                lum_dist=lum_dist,
                                redshift=redshift,
                                time_col=time_col,
                                freq_col=freq_col,
                                flux_col=flux_col,
                                flux_err_col=flux_err_col
                                )

clean_GRB_data(freq_49, flux_col, flux_err_col, time_col)
clean_GRB_data(freq_13, flux_col, flux_err_col, time_col)
clean_GRB_data(freq_85, flux_col, flux_err_col, time_col)

light_curves('980703', 
             Freq_85=freq_85,
             Freq_49=freq_49,
             Freq_13=freq_13,
             time_col=time_col, 
             freq_col=freq_col, 
             flux_col=flux_col,
             flux_err_col=flux_err_col)

upper_13 = 1.43*(1+redshift) 
lower_13 = 1.43*(1+redshift) 
upper_49 = 5.0*(1+redshift) 
lower_49 = 4.6*(1+redshift) 
upper_85 = 8.6*(1+redshift) 
lower_85 = 8.4*(1+redshift) 
upper_150 = 15*(1+redshift) 
lower_150 = 15*(1+redshift)

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

spectral_index_plots('980703', 
                     freq_49_85=spix_49_85,
                     freq_13_49=spix_13_49
                     )

#spix_13_49.to_csv('980703_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('980703spectral_index(4.9-8.5).csv', index=True)



