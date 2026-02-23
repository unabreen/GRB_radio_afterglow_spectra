

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 15:46:07 2026

@author: unabreen
"""
# working with radio data from GRB 030329
# make spectral index plot from 4.9 and 8.5 GHz
# flux scatter plot for 4.9 and 8.5 GHz
#flux in micro janskys

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import linregress
from grb_functions import light_curves, spectral_index_plots
#first columnm observed freq in Hz
#second column is
col_names = ['Frequency(Hz)','Days','Flux(mJy)','Flux uncertainty']

GRB970508 = 'radio970508.dat'   


redshift = 0.835

lum_dist = 5358.8  * 3.08568e22   # Mpc to m    FIX FOR LUMINOSITY CALC

# dataset specific limits on radio bands
upper_13 = 1.41e9*(1 + redshift)
lower_13 = 1.25e9*(1 + redshift)
upper_49 = 5.1e9*(1 + redshift)
lower_49 = 4.75e9*(1 + redshift) 
upper_85 = (8.6e9*(1+ redshift))
lower_85 =  8.4e9 * (1 + redshift)
upper_150 = 1.6e10*(1 + redshift)
lower_150 = 1.45e10*(1 + redshift)

#data frame of all data of GRB 970508
df_970 = pd.read_csv(GRB970508, sep='\s+', header = None, names = col_names)

df_970 = df_970[df_970['Flux(mJy)'] >= 0]
df_970 = df_970[df_970['Flux uncertainty'] >= 0]


def radio_bands(df_all):

    #redshift corrections
    df_all['Days'] = df_all['Days'] / (1 + redshift)
    df_all['Flux(mJy)'] = df_all['Flux(mJy)'] / (1 + redshift)
    df_all['Frequency(Hz)'] = df_all['Frequency(Hz)'] *(1 + redshift)
    
    
    
    # dataframe catches adjusted for redshift
    #dataframe for GRB 030329 for 1.3 GHz (1.25-1.41)
    freq_13 = (df_all[(df_all['Frequency(Hz)'] <= upper_13) &
                         (df_all['Frequency(Hz)'] >= lower_13)]).sort_values(by='Days')
                         
    #dataframe for GRB 030329 for 4.9 GHz (4.75-5.0)                   
    freq_49 = (df_all[(df_all['Frequency(Hz)'] <= upper_49) &
                         (df_all['Frequency(Hz)'] >= lower_49)]).sort_values(by='Days')
    
    #dataframe for GRB 030329 for frequency = 8.5GHz(8.4-8.6)
    freq_85 = (df_all[(df_all['Frequency(Hz)'] <= upper_85 ) &
                         (df_all['Frequency(Hz)'] >= lower_85)]).sort_values(by='Days')
    
    #dataframe for GRB 030329 for frequency = 15Hz
    freq_150 = (df_all[(df_all['Frequency(Hz)'] <= upper_150) &
                         (df_all['Frequency(Hz)'] >= lower_150)]).sort_values(by='Days')
    
    #luminosity calculations
    flux_total = np.trapezoid(df_all['Flux(mJy)'], df_all['Frequency(Hz)'] )
    luminosity = flux_total * 4 * np.pi * (lum_dist)**2
    
    return freq_13, freq_49, freq_85, freq_150, flux_total, luminosity 



def spectral_index(
        df_full,
        lower,
        upper,
        bin_width = 1.0
        
        ):
    
    df = df_full[(df_full['Frequency(Hz)'] >= lower) & 
                 (df_full['Frequency(Hz)'] <= upper)]
    df = df.sort_values(by= 'Days')
    
    t_min, t_max = df['Days'].min(), df['Days'].max()
    bins = np.arange(t_min, t_max, bin_width)
    
    
    results = []
    
    for i in range(len(bins)- 1):
        
        t_low, t_high = bins[i], bins[i+1]
        t_center = (t_low + t_high)/2

        in_bin = df[(df['Days']>= t_low) & (df['Days'] <= t_high)]
        
        if len(np.unique(in_bin['Frequency(Hz)'])) < 2:
            continue
        
        x = np.log(in_bin['Frequency(Hz)'].values)
        y = np.log(in_bin['Flux(mJy)'].values)

        # Error propagation: sigma_logF = sigma_F / F
        #sigma_y = in_bin['Flux uncertainty'].values / in_bin['Flux(mJy)'].values

        # Weighted linear fit: y = alpha * x + C
        #w = 1.0 / sigma_y**2
        
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        
        results.append({'alpha': slope, 
                       'bin center': t_center,
                       'alpha err' : std_err
                       
                       })
    result_frame = pd.DataFrame(results)
    result_frame['Days'] = df['Days']
    return result_frame



freq_13, freq_49, freq_85, freq_150, flux_total, luminosity = radio_bands(df_970)

spix_49_85 = spectral_index(df_970, lower_49, upper_85)
spix_13_49 = spectral_index(df_970, lower_13, upper_49)


light_curves('970508', 
             Freq_85=freq_85, 
             Freq_49=freq_49, 
             Freq_13=freq_13)

spectral_index_plots('970508',
                     freq_49_85=spix_49_85,
                     freq_13_49= spix_13_49)






"""
# Flux plots
fig , ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    nrows = 2, 
    ncols = 2, 
    figsize = (10, 10))
    

#15     no data



#8.5
ax2.errorbar(freq_85['Days'], freq_85['Flux(mJy)'], 
             freq_85['Flux uncertainty'], 
             fmt = 'o', 
             c ='c', 
             label = '8.5 GHz')

ax2.set_title('GRB 970508 Flux- 8.5 GHz')
ax2.set_xlabel('Days')
ax2.set_ylabel('Flux(mJy)')
ax2.set_xscale('log')
ax2.set_yscale('log')





#4.9
ax3.errorbar(freq_49['Days'], freq_49['Flux(mJy)'], 
             yerr= freq_49['Flux uncertainty'], 
             fmt = 'o', 
             c = 'r', 
             label = '4.9 GHz')


ax3.set_title('GRB 970508 Flux- 4.9 GHz')
ax3.set_xlabel('Days')
ax3.set_ylabel('Flux(mJy)')
ax3.set_xscale('log')
ax3.set_yscale('log')


# 1.3

ax4.errorbar(freq_13['Days'], freq_13['Flux(mJy)'], 
             freq_13['Flux uncertainty'], 
             fmt = 'o', 
             c ='g', 
             label = '1.3 GHz')

ax4.set_title('GRB 970508 Flux- 1.3 GHz')
ax4.set_xlabel('Days')
ax4.set_ylabel('Flux(mJy)')
ax4.set_xscale('log')
ax4.set_yscale('log')



#spectral index plots

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    nrows = 2, 
    ncols = 2, 
    figsize = (10, 10))


#8.5- 15   # no data



# 4.9-8.5
ax2.errorbar(x= spix_49_85['Days'], 
             y=spix_49_85['alpha'],
             yerr = spix_49_85['alpha err'],
             fmt = 'o',
             c = 'c'
             
             )
ax2.set_title("GRB 970508- spectral index 4.9-8.5 GHz")
ax2.set_xscale('log')
ax2.set_xlabel('Days')
ax2.set_ylabel('Spectral index')
ax2.axhline(y=2, linestyle='--')
ax2.axhline(y=(1/3), linestyle='--')
ax2.axhline(y=-0.6, linestyle='--')


#1.3-4.9
ax3.errorbar(x= spix_13_49['Days'], 
             y=spix_13_49['alpha'],
             yerr = spix_13_49['alpha err'],
             fmt = 'o',
             c = 'r'
             
             )
ax3.set_title("GRB 970508- spectral index 1.3-4.9 GHz")
ax3.set_xscale('log')
ax3.set_xlabel('Days')
ax3.set_ylabel('Spectral index')
ax3.axhline(y=2, linestyle='--')
ax3.axhline(y=(1/3), linestyle='--')
ax3.axhline(y=-0.6, linestyle='--')

"""










