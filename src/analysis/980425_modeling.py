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
from grb_functions import radio_bands
from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean
import matplotlib.pyplot as plt


col_names = ['name',
             'telescope',
             'Month',
             'Column 3',  # unknown, named as column with index 3
             'year',
             'Days',    
             'Frequency(GHz)',  # GHz 
             'Flux(mJy)',   # possibly, unknown units
             'Flux uncertainty',  #possibly
             'Column 9'   # unknown, named as column with index 9
             ]


flux_col = 'Flux(mJy)'
flux_err_col = 'Flux uncertainty'
freq_col = 'Frequency(GHz)'
time_col = 'Days'


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









#astropy fit

from astropy.modeling.models import BrokenPowerLaw1D
from astropy.modeling.fitting import LMLSQFitter

# data to fit
x_data = freq_13[time_col]
y_data = freq_13[flux_col]
y_err_data = freq_13[flux_err_col]

x_max = x_data.max()
x_min = x_data.min()
nu_norm = np.median(y_data)
break_time = freq_13.loc[y_data.idxmax(), time_col]


initial_model = BrokenPowerLaw1D(
    amplitude = x_min,  #look at initial guess here
    x_break = break_time,
    alpha_1 = 1,
    alpha_2 = -1 )

#no masks for 1.3 GHz- 1 peak

fitter = LMLSQFitter()
fitted = fitter(initial_model, x_data, y_data )

x_fit = np.linspace(x_min, x_max, 500)
y_fit = fitted(x_fit)

slope_1 = fitted.alpha_1.value * -1  # slope before the break
slope_2 = fitted.alpha_2.value  * -1 # slope after the break
amplitude = fitted.amplitude.value
x_break = fitted.x_break.value

# plot with fitted line
fig, ax = plt.subplots()

ax.errorbar(x= freq_13[time_col], y = freq_13[flux_col], 
            yerr= freq_13[flux_err_col], fmt = 'o')
ax.set_title(f'1.3 GHz band with fit\nslope 1 = {slope_1:.4f}, slope 2 = {slope_2:.4f}')
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(x_fit, y_fit, "r-", label="Broken power law fit")


"""
two peaks fitter function

Currently manually breaking sets- have to automate break points
break 4.9 GHz set at 21.4434 seconds
break 8.5 GHz set at 23.4361
freq_set is frame, set_break is float of time break, freq_name is string of 
frequency in GHz
"""
break_set_49 = 21.4434
break_set_85 = 23.4361
name_49 = 4.9
name_85 = 8.5
def light_curve_fit(freq_set, set_break, set_name): 
    
    
    first_peak_mask = freq_set[time_col] < set_break
    second_peak_mask = freq_set[time_col] >= set_break
    first_peak = freq_set[first_peak_mask]
    second_peak = freq_set[second_peak_mask]
    
    x_data_1 = first_peak[time_col]
    y_data_1 = first_peak[flux_col]
    
    x_data_2 = second_peak[time_col]
    y_data_2 = second_peak[flux_col]
    
    x_max_1 = x_data_1.max()
    x_min_1 = x_data_1.min()
    nu_norm_1 = y_data_1.max()  #using max instead of median due to data spread
    break_time_1 = first_peak.loc[y_data_1.idxmax(), time_col]
   
    x_min_2 = x_data_2.min()
    nu_norm_2 = y_data_2.max()
    break_time_2 = second_peak.loc[y_data_2.idxmax(), time_col]
    x_max_2 = ((break_time_2-x_min_2) * 2) + x_min_2
    
    initial_model_1 = BrokenPowerLaw1D(
        amplitude = nu_norm_1,
        x_break = break_time_1,
        alpha_1 = 1,
        alpha_2 = -1 )
    
    initial_model_2 = BrokenPowerLaw1D(
        amplitude = nu_norm_2,
        x_break = break_time_2,
        alpha_1 = 1,
        alpha_2 = -1 )
    
    
    fitter = LMLSQFitter()
    
    #first peak fit
    fitted_1 = fitter(initial_model_1, x_data_1, y_data_1, filter_non_finite= True)
    
    x_fit_1 = np.linspace(x_min_1, x_max_1, 500)
    y_fit_1 = fitted_1(x_fit_1)
    
    #second peak fit
    fitted_2 = fitter(initial_model_2, x_data_2, y_data_2, filter_non_finite= True )
    
    x_fit_2 = np.linspace(x_min_2, x_max_2, 500)
    y_fit_2 = fitted_2(x_fit_2)
    
    first_slope_1 = fitted_1.alpha_1.value *-1  # slope before the break in first peak
    first_slope_2 = fitted_1.alpha_2.value  *-1 # slope after the break in first peak
    first_amplitude = fitted.amplitude.value
    first_x_break = fitted.x_break.value
    
    second_slope_1 = fitted_2.alpha_1.value  *-1 # slope before the break in second peak
    second_slope_2 = fitted_2.alpha_2.value  *-1 # slope after the break in second peak
    second_amplitude = fitted.amplitude.value
    second_x_break = fitted.x_break.value
    
    # plot with fitted line
    fig, ax = plt.subplots()
    
    ax.errorbar(x= first_peak[time_col], y = first_peak[flux_col], 
                yerr= first_peak[flux_err_col], fmt = 'o', zorder =2)
    ax.errorbar(x= second_peak[time_col], y = second_peak[flux_col], 
                yerr= second_peak[flux_err_col], fmt = 'o', c = 'r', zorder=2)
    ax.set_title(f' {set_name} GHz band with fit\nfirst slope set 1 = {first_slope_1:.3f}, first slope set 2 = {first_slope_2:.3f},\nsecond slope set 1 = {second_slope_1:.3f} , second slope set 2 = {second_slope_2:.3f}')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(x_fit_1, y_fit_1, "b-", label="Broken power law fit 1", zorder=5)
    ax.plot(x_fit_2, y_fit_2, "k-", label="Broken power law fit 2", zorder=5)
    ax.legend()
    plt.show()
    
    from collections import namedtuple

    PeakFit = namedtuple('PeakFit', ['a1', 'a2', 'amp', 'brk'])
    
    first_peak_fits = PeakFit(a1=first_slope_1, a2=first_slope_2, amp=first_amplitude, brk=first_x_break)
    second_peak_fits = PeakFit(a1=second_slope_1, a2=second_slope_2, amp=second_amplitude, brk=second_x_break)


    return fitted_1, fitted_2
    
first_peak_fits_49, second_peak_fits_49 = light_curve_fit(freq_49, break_set_49, name_49)
first_peak_fits_85, second_peak_fits_85 = light_curve_fit(freq_85, break_set_85, name_85)

print(first_peak_fits_49, second_peak_fits_49)

# superposition test
break_time_1 = first_peak_fits_49.x_break
break_time_2 = second_peak_fits_49.x_break

x_fit_1 = np.linspace(freq_49[time_col].min(), break_time_1, 500)
x_fit_2 = np.linspace(break_time_1, break_time_2, 500)
x_fit_3 = np.linspace(break_time_2, freq_49[time_col].max(), 500)

y_fit_1_1 = first_peak_fits_49(x_fit_1)
y_fit_1_2 = second_peak_fits_49(x_fit_1)
y_fit_1 = y_fit_1_1 + y_fit_1_2

y_fit_2_1 = first_peak_fits_49(x_fit_2)
y_fit_2_2 = second_peak_fits_49(x_fit_2)
y_fit_2 = y_fit_2_1 + y_fit_2_2

y_fit_3_1 = first_peak_fits_49(x_fit_3)
y_fit_3_2 = second_peak_fits_49(x_fit_3)
y_fit_3 = y_fit_3_1 + y_fit_3_2

fig, ax = plt.subplots()
ax.errorbar(x= freq_49[time_col], y = freq_49[flux_col], 
            yerr= freq_49[flux_err_col], fmt = 'o')
ax.set_title('4.9 GHz with superposition of models between peaks')
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(x_fit_1, y_fit_1, "k-", label="Broken power law fit 1", zorder=5)
ax.plot(x_fit_2, y_fit_2, "k-", label="Broken power law fit 2- superposition", zorder=5)
ax.plot(x_fit_3, y_fit_3, "k-", label="Broken power law fit 3", zorder=5)
plt.show()

