#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 14:30:30 2026

@author: unabreen
"""

import sys
from pathlib import Path

# Points to .../src/analysis, relative to this script's own location
sys.path.append(str(Path(__file__).resolve().parent.parent / "analysis"))


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from grb_functions import clean_GRB_data
from grb_functions import uncertainty_clean
from scipy.optimize import curve_fit



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





"""
Free parameters to fit:
    Fnu_a(calculated), nu_a(calcualted), p, alpha_nu, alpha_F
"""

# assumes freq_col and flux_col are established
def spectra_model(df, #GRB dataframe
                  alpha_F, #guess value
                  alpha_nu, #guess value
                  p,
                  x = freq_col,
                  y = flux_col,
                  t = time_col
                  ):  # p guess
    
    
    nu = df[x]
    F= df[y]
    time = df[t]
    Fnu_a_init = F[0]
    nu_a_init = nu[0]
    t_init = time[0]
    
    
    Fnu_a = Fnu_a_init * (time/t_init)** -alpha_F
    nu_a = nu_a_init * ((time/t_init)** -alpha_nu)
    
    alpha1 = 2.5
    alpha2 = -(p-1)/2
    
    return np.where(nu < nu_a, 
             Fnu_a * (nu/nu_a)** alpha1,
             Fnu_a * (nu/nu_a)** alpha2
             )

z = spectra_model(df = df, alpha_F= 1, alpha_nu = 1, p = 2)  

   
#set x data for fit reference later
# y data goes in fitter function
x_data = df[freq_col].to_numpy()
y_data = df[flux_col].to_numpy()
yerr = df[flux_err_col].to_numpy()

#initial guesses: alpha_F= 1, alpha_nu = 1, p = 2
# x data is entire dataframe for model to work
popt, pcov = curve_fit(spectra_model, 
                       df, 
                       y_data,
                       p0 = [1, 1, 2],
                       sigma= yerr)




x_fit = np.linspace(x_data.min(), x_data.max(), 100)
alpha_F, alpha_nu, p = popt

y_fit = spectra_model(df, alpha_F, alpha_nu, p)

#calculating chi squared
y = df[[flux_col]]
y_error = df[[flux_err_col]]
residuals = (y_data - y_fit)/yerr
chi2 = np.sum(residuals**2)
dof = y_data - len(popt)
reduced_chi2 = chi2/dof
    

mask = df[time_col].between(11, 12, inclusive="neither")
row_indices = df.index[mask]
test_epoch = df[mask]

y_test = y_fit[mask]

fig, ax = plt.subplots()

ax.scatter(test_epoch[freq_col], test_epoch[flux_col])
ax.plot(test_epoch[freq_col], y_test)
ax.set_xlabel('Frequency(GHz)')
ax.set_ylabel('Flux(mJy)')
ax.set_xscale('log')
ax.set_yscale('log')


        
    