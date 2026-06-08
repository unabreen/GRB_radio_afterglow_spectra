#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:13:58 2026

@author: unabreen

This script takes the spectral data calculated from the GRB datasets
which I have decided have sufficient data to analyze. I define a function
that plots the spectral data in sets corresponding to the defined frequency
bands (1.3-4.9, 4.9-8.5, 8.5-15). The spectral data files are sorted into 
sets to plot depending on the frequency band. the 4.9-8.5 band is split into 
two sets in order for the plots to be legible. Colors known to contrast well 
on plot are assigned to each GRB for consistent comparison across plots. 
In addition to the plot combining temporal evolution of spectral index 
directly, I plot the spectral index datasets with their centers shifted so the 
patterns can be observed visually. This is conducted by locating the data 
point closest to zero(in the middle of the spectral index evolution) and 
scaling the time column(bin center) so that the zero index maps to t= 1 day.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from grb_functions import spectral_index_plots


# function to plot spectral data files
# takes the set of GRB names and the radio band to locate the file and plot 
# each dataset
# takes radio band in form of '1.3-4.9')
def plot_spectral_data(GRB_set, band, GRB_color, split = False):

    #importing spectral index files for the 1.3-4.9 GHz band
    spectra_dict = {}
    for GRB in GRB_set:
        filename = f'spectral data files/{GRB} spectra/{GRB}_spectral_index({band}).csv'
        spectra_dict[GRB] = pd.read_csv(filename, sep=',')
    
   
    
    #spectral index plot- generate figure and then plot each spectral dataset 
    # within  the loop
    
    # creates a new column in each frame with times centered around t=1 day
    for GRB in spectra_dict:
        df = spectra_dict[GRB]
        n = len(df)
        lower = int(n * 0.1)
        upper = int(n * 0.75)
        
        # Find the index where alpha is closest to 0, restricted to the middle 50%
        middle_df = df.iloc[lower:upper]
        zero_idx = middle_df['alpha'].abs().idxmin()
        zero_bin_center = df.loc[zero_idx, 'bin center']
        
        
        # Calculate the shift needed so that zero_bin_center maps to 10e0(1) = 10.0
        target_value = 10
        scale = target_value/zero_bin_center
        
        df['shifted bin center'] = df['bin center'] *scale 
    
    
    fig, (ax1) = plt.subplots(figsize = (10, 10))
    for GRB in spectra_dict:
        frame = spectra_dict[GRB]
        ax1.errorbar(x= frame['bin center'], 
                     y= frame['alpha'],
                     yerr = frame['alpha err'],
                     fmt = 'o',
                     c = GRB_color[GRB],
                     label = f"GRB {GRB}",
                     markersize = 8
                     )
    ax1.set_title(f"spectral indices {band} GHz- All GRBs", fontsize = 20)
    ax1.set_xscale('log')
    ax1.set_xlabel('Days')
    ax1.set_ylabel('Spectral index')
    ax1.axhline(y=2, linestyle='--')
    ax1.axhline(y=(1/3), linestyle='--')
    ax1.axhline(y=-0.6, linestyle='--')
    ax1.legend(
        fontsize=14,     
        markerscale=2) 
    ax1.tick_params(axis='both', labelsize=14)
    
    plt.tight_layout()
    
 
    
    
    
    fig, (ax1) = plt.subplots(figsize = (10, 10))
    for GRB in spectra_dict:
        frame = spectra_dict[GRB]
        ax1.errorbar(x= frame['shifted bin center'], 
                     y= frame['alpha'],
                     yerr = frame['alpha err'],
                     fmt = 'o',
                     c = GRB_color[GRB],
                     label = f"GRB {GRB}",
                     markersize = 8
                     )
    ax1.set_title(f"spectral indices {band} GHz- All GRBs- centered around alpha = 0", 
                  fontsize = 20)
    ax1.set_xscale('log')
    ax1.set_xlabel('Days')
    ax1.set_ylabel('Spectral index')
    ax1.axhline(y=2, linestyle='--')
    ax1.axhline(y=(1/3), linestyle='--')
    ax1.axhline(y=-0.6, linestyle='--')
    ax1.axvline(x=10, linestyle='--')
    ax1.legend(
        fontsize=14,     
        markerscale=2) 
    ax1.tick_params(axis='both', labelsize=14)
    
    plt.tight_layout()
 
    plt.show()
    
    return(spectra_dict)


GRB_color = {'030329': 'red', 
             '221009A': 'turquoise',      
             '980703': 'orange', 
             '970508': 'green',      
             '130427A': 'goldenrod',
             '980425': 'blue',
             '000418': 'purple',
             "021004" : 'pink', 
             "030723" : 'black',
             "991208" : 'dodgerblue',
             '980329' : 'crimson'
             }



GRB_set_13_49 = ['030329', '221009A', '980703', '970508', '980425']
GRB_set_49_85_one = ['030329', '221009A', '980703', '970508', 
                 '980425']
GRB_set_49_85_two = ['130427A', '000418', '021004', '030723', '991208',
                    '980329' ] 
GRB_set_85_150 = ['030329', '221009A', '991208' ]





    
    
frames_13_49 = plot_spectral_data(GRB_set_13_49, '1.3-4.9', GRB_color)
frames_49_85 = plot_spectral_data(GRB_set_49_85_one, '4.9-8.5', GRB_color, 
                                  split = 'one')
plot_spectral_data(GRB_set_49_85_two, '4.9-8.5', GRB_color,
                   split = 'two')    



# no spectral data for 8.5- 15 for  
plot_spectral_data(GRB_set_85_150, '8.5-15', GRB_color)

    
    
    
    
