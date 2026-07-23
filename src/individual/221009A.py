#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 14:42:31 2026

@author: unabreen

working with radio data from GRB 221009A
make spectral index plot from 4.9 and 8.5 GHz
flux scatter plot for 4.9 and 8.5 GHz

column 1: ignore (time in Modified Julian Date)
- column 2: observing frequency in GHz
- column 3: flux in milliJy
- column 4: ignore (statistical uncertainty only)
- column 5: uncertainty in milli Jy
- column 6: ignore (telescope name)
- column 7: time in days
"""

import pandas as pd
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean

col_names = ['time(Modified Julian Date)', 
             'Frequency(Hz)', 
             'Flux(mJy)', 
             'statistical uncertainty flux', 
             'Flux uncertainty', 
             'col_5(placeholder)',
             'telescope name', 
             'Days']

filepath = '../old_grb_sample/all_data.csv'

redshift = 0.151
lum_dist = 723.6 * 3.08568e22   # Mpc to m

# specific radio band limits
upper_13 = 1.41e9*(1 + redshift)
lower_13 = 1.25e9*(1 + redshift)
upper_49 = 5.2e9*(1 + redshift)
lower_49 = 4.75e9*(1 + redshift)
upper_85 = 8.6e9*(1 + redshift)
lower_85 = 8.4e9*(1 + redshift)
upper_150 = 1.61e10*(1 + redshift)
lower_150 = 1.39e10*(1 + redshift)

df = pd.read_csv(filepath, header=None, names=col_names)

time_col = 'Days'
flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(Hz)'

# Drop rows where any object (string) column contains 'some_value'
# specific to 221009A data file- missing values have '-', 
# Days column not numeric
obj_cols = df.select_dtypes(include=['object']).columns
df = df[~df[obj_cols].eq('-').any(axis=1)]

df['Days'] = pd.to_numeric(df['Days'], errors='coerce')

df['Frequency(Hz)'] = df['Frequency(Hz)'] * 10e8

# Set display options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)

uncertainty_clean(df, flux_col, flux_err_col)
clean_GRB_data(df, flux_col, flux_err_col, time_col)

freq_13, freq_49, freq_85, freq_150, df = radio_bands(
                                df,  
                               upper_13=1.41e9*(1 + redshift),
                               lower_13=1.25e9*(1 + redshift),
                               upper_49=5.2e9*(1 + redshift),
                               lower_49=4.75e9*(1 + redshift),
                               upper_85=8.6e9*(1 + redshift),
                               lower_85=8.4e9*(1 + redshift),
                               upper_150=1.61e10*(1 + redshift),
                               lower_150=1.39e10*(1 + redshift),
                                lum_dist=lum_dist,
                                redshift=redshift,
                                time_col=time_col,
                                freq_col=freq_col,
                                flux_col=flux_col,
                                flux_err_col=flux_err_col
                                )

light_curves('221009A', 
             Freq_150=freq_150,
             Freq_85=freq_85,
             Freq_49=freq_49,
             Freq_13=freq_13,
             time_col=time_col, 
             freq_col=freq_col, 
             flux_col=flux_col,
             flux_err_col=flux_err_col)

spix_85_150 = spectral_index(df, 
                             lower_85, 
                             upper_150,
                             freq_col=freq_col,
                             time_col=time_col,
                             flux_col=flux_col
                             )    

spix_49_85 = spectral_index(df, lower_49, upper_85,
                            freq_col=freq_col,
                            time_col=time_col,
                            flux_col=flux_col
                            )
spix_13_49 = spectral_index(df, lower_13, upper_49,
                            freq_col=freq_col,
                            time_col=time_col,
                            flux_col=flux_col
                            )

spix_85_150 = spix_85_150[(spix_85_150['alpha'] <= 20)]
spix_85_150 = spix_85_150[(spix_85_150['alpha'] >= -3)]

spectral_index_plots('221009A', spix_85_150, spix_49_85, spix_13_49)

#spix_13_49.to_csv('spectral data files/221009A spectra/221009A_spectral_index(1.3-4.9).csv', index=True)
#spix_49_85.to_csv('spectral data files/221009A spectra/221009A_spectral_index(4.9-8.5).csv', index=True)
#spix_85_150.to_csv('spectral data files/221009A spectra/221009A_spectral_index(8.5-15).csv', index=True)










