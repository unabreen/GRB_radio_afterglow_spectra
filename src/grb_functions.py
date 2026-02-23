#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 20:59:08 2026

@author: unabreen

Storing all the functions I use to view GRB data in this script
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress


def radio_bands(
        df_all,
        redshift,
        upper_13,
        lower_13,
        upper_49,
        lower_49,
        upper_85,
        lower_85,
        upper_150,
        lower_150,
        lum_dist

                ):
    
    # Keep only rows where all numeric columns have valid numbers
    numeric_cols = df_all.select_dtypes(include=['number']).columns
    df_all = df_all[df_all[numeric_cols].notna().all(axis=1)]

    #redshift corrections
    df_all['Days'] = df_all['Days'] / (1 + redshift)
    df_all['Flux(erg cm^-2'] = df_all['Flux(erg cm^-2'] / (1 + redshift)
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
    flux_total = np.trapezoid(df_all['Flux(erg cm^-2'], df_all['Frequency(Hz)'] )
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
        y = np.log(in_bin['Flux(erg cm^-2'].values)

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




def light_curves(
        grb_name,
        Freq_150 = None,
        Freq_85 = None,
        Freq_49 = None,
        Freq_13 = None,
        ):
    
    fig , ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        nrows = 2, 
        ncols = 2, 
        figsize = (10, 10))
        
    
    #15
    if Freq_150 is not None:
        ax1.errorbar(Freq_150['Days'], Freq_150['Flux(erg cm^-2'], 
                     Freq_150['Flux uncertainty'], 
                     fmt = 'o', 
                     c ='b', 
                     label = '15 GHz')
        
        ax1.set_title(f'GRB {grb_name} Flux- 15 GHz')
        ax1.set_xlabel('Days')
        ax1.set_ylabel('Flux(erg cm^-2')
        ax1.set_xscale('log')
        ax1.set_yscale('log')
    
    
    
    #8.5
    if Freq_85 is not None:
        ax2.errorbar(Freq_85['Days'], Freq_85['Flux(erg cm^-2'], 
                     Freq_85['Flux uncertainty'], 
                     fmt = 'o', 
                     c ='c', 
                     label = '8.5 GHz')
        
        ax2.set_title(f'GRB {grb_name} Flux- 8.5 GHz')
        ax2.set_xlabel('Days')
        ax2.set_ylabel('Flux(erg cm^-2')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
    
    #4.9
    if Freq_49 is not None:
        ax3.errorbar(Freq_49['Days'], Freq_49['Flux(erg cm^-2'], 
                     yerr= Freq_49['Flux uncertainty'], 
                     fmt = 'o', 
                     c = 'r', 
                     label = '4.9 GHz')
        
        
        ax3.set_title(f'GRB {grb_name} Flux- 4.9 GHz')
        ax3.set_xlabel('Days')
        ax3.set_ylabel('Flux(erg cm^-2')
        ax3.set_xscale('log')
        ax3.set_yscale('log')
    
    
    # 1.3
    if Freq_13 is not None:
        ax4.errorbar(Freq_13['Days'], Freq_13['Flux(erg cm^-2'], 
                     Freq_13['Flux uncertainty'], 
                     fmt = 'o', 
                     c ='g', 
                     label = '1.3 GHz')
        
        ax4.set_title(f'GRB {grb_name} Flux- 1.3 GHz')
        ax4.set_xlabel('Days')
        ax4.set_ylabel('Flux(erg cm^-2')
        ax4.set_xscale('log')
        ax4.set_yscale('log')



def spectral_index_plots(
        grb_name,
        freq_85_150= None,
        freq_49_85= None,
        freq_13_49= None
        
        ):
    #spectral index plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        nrows = 2, 
        ncols = 2, 
        figsize = (10, 10))
    
    
    #8.5- 15
    if freq_85_150 is not None:
        ax1.errorbar(x= freq_85_150['bin center'], 
                     y=freq_85_150['alpha'],
                     yerr = freq_85_150['alpha err'],
                     fmt = 'o',
                     c = 'b'
                     
                     )
        ax1.set_title(f"GRB {grb_name}- spectral index 8.5-15 GHz")
        ax1.set_xscale('log')
        ax1.set_xlabel('Days')
        ax1.set_ylabel('Spectral index')
        ax1.axhline(y=2, linestyle='--')
        ax1.axhline(y=(1/3), linestyle='--')
        ax1.axhline(y=-0.6, linestyle='--')
    
    # 4.9-8.5
    if freq_49_85 is not None:
        ax2.errorbar(x= freq_49_85['bin center'], 
                     y=freq_49_85['alpha'],
                     yerr = freq_49_85['alpha err'],
                     fmt = 'o',
                     c = 'c'
                     
                     )
        ax2.set_title(f"GRB {grb_name}- spectral index 4.9-8.5 GHz")
        ax2.set_xscale('log')
        ax2.set_xlabel('Days')
        ax2.set_ylabel('Spectral index')
        ax2.axhline(y=2, linestyle='--')
        ax2.axhline(y=(1/3), linestyle='--')
        ax2.axhline(y=-0.6, linestyle='--')
    
    #1.3-4.9
    if freq_13_49 is not None:
        ax3.errorbar(x= freq_13_49['bin center'], 
                     y=freq_13_49['alpha'],
                     yerr = freq_13_49['alpha err'],
                     fmt = 'o',
                     c = 'r'
                     
                     )
        ax3.set_title(f"GRB {grb_name}- spectral index 1.3-4.9 GHz")
        ax3.set_xscale('log')
        ax3.set_xlabel('Days')
        ax3.set_ylabel('Spectral index')
        ax3.axhline(y=2, linestyle='--')
        ax3.axhline(y=(1/3), linestyle='--')
        ax3.axhline(y=-0.6, linestyle='--')







def luminosity_func(df_all,
               lum_dist,
               redshift):
    df_all["Luminosity"] = df_all['Flux(erg cm^-2'] * 4 * np.pi * (lum_dist/(1 + redshift))**2













