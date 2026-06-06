#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:14:38 2026

@author: unabreen

This script is focused on modeling the light curves of GRB 980425 to 
investigate their double peaks. It first follows the established pipeline 
of this project to create the light curves for each frequency. The fitting uses
the linear least squares fitter in astropy.modeling, as well as the one 
dimensional broken power law function for the 1.3 GHz single peak.

The next step is to model the double peaks of 4.9 and 8.5 GHz, which I separate
into early and late times. The double peak fitter is currently has bugs which
I am investigating.
"""



import pandas as pd
import numpy as np
from grb_functions import light_curves
from grb_functions import spectral_index_plots
from grb_functions import radio_bands
from grb_functions import spectral_index
from grb_functions import luminosity_func
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean
import matplotlib.pyplot as plt


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


filename = '../old_grb_sample/980425.dat'


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

# curve fits- light curves
# 1.3 GHz
peak_flux = freq_13[flux_col].max()
peak_time = freq_13.loc[freq_13[flux_col].idxmax(), time_col]

# data to fit
x_data = freq_13[time_col]
y_data = freq_13[flux_col]


#astropy fit

from astropy.modeling.models import BrokenPowerLaw1D
from astropy.modeling.fitting import LevMarLSQFitter


initial_model_13 = BrokenPowerLaw1D(
    amplitude = peak_flux,
    x_break = peak_time,
    alpha_1 = 1.2,
    alpha_2 = -1.2 )

#no masks for 1.3 GHz- 1 peak

fitter = LevMarLSQFitter()
fitted = fitter(initial_model_13, x_data, y_data )

x_fit = np.linspace(freq_13[time_col].min(), freq_13[time_col].max())
y_fit = fitted(x_fit)

fig, ax = plt.subplots()

ax.errorbar(x= freq_13[time_col], y = freq_13[flux_col], 
            yerr= freq_13[flux_err_col], fmt = 'o')
ax.set_title('1.3 GHz band with fit')
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(x_fit, y_fit, "r-", label="Broken power law fit")



#4.9 light curve fit

peak_flux_49 = freq_49[flux_col].max()
peak_time_49 = freq_49.loc[freq_49[flux_col].idxmax(), time_col]

second_peak_table = freq_49[freq_49[time_col] > 21.44]
second_peak_flux_49 = second_peak_table[flux_col].max()
second_peak_time_49 = freq_49.loc[second_peak_table[flux_col].idxmax(), time_col]

min_peaks_49 = freq_49[(freq_49[time_col] > peak_time_49) & (freq_49[time_col] > peak_flux_49)].max()

x_data_49 = freq_49[freq_col]
y_data_49 = freq_49[flux_col]

break1 = peak_time_49
break2 = 21.443
break3= 32.36

initial_model_49_1 = BrokenPowerLaw1D(
    amplitude = peak_flux_49,
    x_break = break1,
    alpha_1 = 1.2,
    alpha_2 = -1.2 )

initial_model_49_2 = BrokenPowerLaw1D(
    amplitude = 24,
    x_break = break2,
    alpha_1 = -1,
    alpha_2 = 1 )

initial_model_49_3 = BrokenPowerLaw1D(
    amplitude = 30,
    x_break = break3,
    alpha_1 = 1,
    alpha_2 = -1.2 )

#masks to filter fit bands
mask1 = x_data < break1
mask2= (x_data >= break1) & x_data < break2
mask3 = x_data >= break2

fitter = LevMarLSQFitter()
fit1 = fitter(initial_model_49_1, x_data[mask1], y_data[mask1])
fit2 = fitter(initial_model_49_2, x_data[mask2], y_data[mask2])
fit3 = fitter(initial_model_49_3, x_data[mask3], y_data[mask3])

x_fit_1 = np.linspace(freq_49[time_col].min(), peak_time_49)
y_fit_1 = fitted(x_fit_1)

x_fit_2 = np.linspace(peak_time_49, second_peak_time_49)
y_fit_2 = fitted(x_fit_2)

x_fit_3 = np.linspace(second_peak_time_49, freq_49[time_col].max())
y_fit_3 = fitted(x_fit_3)


fig, ax = plt.subplots()

ax.errorbar(x= freq_49[time_col], y = freq_49[flux_col], 
            yerr= freq_13[flux_err_col], fmt = 'o')
ax.set_title('4.9 GHz band with fit')
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(x_fit_1, y_fit_1, "r-", label="Broken power law fit segment 1")
ax.plot(x_fit_2, y_fit_2, "r-", label="Broken power law fit segment 2")
ax.plot(x_fit_3, y_fit_3, "r-", label="Broken power law fit segment 3")







#early time spectrum
fig, ((ax1, ax2)) = plt.subplots(nrows = 1, ncols = 2, figsize= (10, 5))


early_df = df[df[time_col] < 16]
late_df = df[df[time_col] > 70]

ax1.scatter(early_df[freq_col], early_df[flux_col])
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title('Spectrum early times')

ax2.scatter(late_df[freq_col], late_df[flux_col])
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title('Spectrum late times')

















