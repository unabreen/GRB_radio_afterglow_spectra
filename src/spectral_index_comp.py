#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:13:58 2026

@author: unabreen
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from grb_functions import spectral_index_plots


# collection of GRB spectral index files(extracted from raw data earlier)
# will end up plotting this collection on the same figure
GRB_set = ['030329', '221009A', '031203', '980703', '970508']


# function to plot spectral data files
# takes the set of GRB names and the radio band to locate the file and plot 
# each dataset
# takes radio band in form of '1.3-4.9'
def plot_spectral_data(GRB_set, band, GRB_color):

    #importing spectral index files for the 1.3-4.9 GHz band
    spectra_dict = {}
    for GRB in GRB_set:
        filename = f'spectral data files/{GRB} spectra/{GRB}_spectral_index({band}).csv'
        spectra_dict[GRB] = pd.read_csv(filename, sep=',')
    
    # color generation
    #colors = cm.hsv(np.linspace(0, 0.9, len(GRB_set)))  # stop at 0.9 to avoid red repeating
    
    #spectral index plot- generate figure and then plot each spectral dataset 
    # within  the loop
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
    
    plt.tight_layout()
    plt.show()
    


GRB_color = {'030329': 'red', 
             '221009A': 'turquoise',      
             '980703': 'orange', 
             '970508': 'goldenrod',      
             '130427A': 'green',
             '980425': 'blue',
             '000418': 'purple'
             }



GRB_set_13_49 = ['030329', '221009A', '980703', '970508', '980425']
GRB_set_49_85 = ['030329', '221009A', '130427A','980703', '970508', 
                 '980425', '000418']
GRB_set_85_150 = ['030329', '221009A', ]

plot_spectral_data(GRB_set_13_49, '1.3-4.9', GRB_color)
plot_spectral_data(GRB_set_49_85, '4.9-8.5', GRB_color)   


# no spectral data for 8.5- 15 for  
plot_spectral_data(GRB_set_85_150, '8.5-15', GRB_color)

    
    
