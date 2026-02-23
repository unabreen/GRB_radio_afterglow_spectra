#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 00:17:29 2026

@author: unabreen
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


    
    #luminosity = flux_total * 4 * np.pi * (lum_dist)**2


#030329
col_names_030 = ['Frequency(Hz)','Days','Flux(erg cm^-2','Flux uncertainty']

GRB030329 = 'radio030329.dat'   

redshift_030 = 0.1685


# radio specific band limits


lum_dist_030 = 816.4  * 3.08568e20   # Mpc to cm

df_030 = pd.read_csv(GRB030329, sep='\s+', header = None, names = col_names_030)
#milli Janskys to ergs
df_030['Flux(erg cm^-2'] = df_030['Flux(erg cm^-2'] * 10e-26
df_030['Flux uncertainty'] = df_030['Flux uncertainty'] * 10e-26

freq_13_030, freq_49_030, freq_85_030, freq_150_030, flux_total, luminosity= radio_bands(df_030,
            redshift=redshift_030,
            upper_13 = 1.41e9*(1 + redshift_030),
            lower_13 = 1.25e9*(1 + redshift_030),
            upper_49 = 5.1e9*(1 + redshift_030),
            lower_49 = 4.75e9*(1 + redshift_030),
            upper_85 = 8.7e9*(1 + redshift_030),
            lower_85 = 8.4e9*(1 + redshift_030),
            upper_150 = 1.6e10*(1 + redshift_030),
            lower_150 = 1.45e10*(1 + redshift_030),
            lum_dist=lum_dist_030
            )


#221009A

col_names_221 = ['time(Modified Julian Date)', 
             'Frequency(Hz)', 
             'Flux(erg cm^-2', 
             'statistical uncertainty flux', 
             'Flux uncertainty', 
             'idk',
             'telescope name', 
             'Days']

GRB221009A = 'all_data.csv'

redshift_221 = 0.151
lum_dist_221 = 723.6 * 3.08568e20   # Mpc to m
df_221 = pd.read_csv(GRB221009A, header = None, names = col_names_221)
# specific to 221009A data file- missing values have '-', 
#Days column not numeric
obj_cols = df_221.select_dtypes(include=['object']).columns
df_221 = df_221[~df_221[obj_cols].eq('-').any(axis=1)]

df_221['Days'] = pd.to_numeric(df_221['Days'], errors='coerce')
df_221['Frequency(Hz)'] = df_221['Frequency(Hz)'] * (10**9)

df_221['Flux(erg cm^-2'] = df_221['Flux(erg cm^-2'] * 10e-26
df_221['Flux uncertainty'] = df_221['Flux uncertainty'] * 10e-26


freq_13_221, freq_49_221, freq_85_221, freq_150_221, flux_total, luminosity =radio_bands(df_221, 
            upper_13 = 1.41e9*(1 + redshift_221),
            lower_13 = 1e9*(1 + redshift_221),
            upper_49 = 5.2e9*(1 + redshift_221),
            lower_49 = 4.75e9*(1 + redshift_221),
            upper_85 = 8.6e9*(1 + redshift_221),
            lower_85 = 8.4e9*(1 + redshift_221),
            upper_150 = 1.61e10*(1 + redshift_221),
            lower_150 = 1.39e10*(1 + redshift_221),
            lum_dist=lum_dist_221,
            redshift= redshift_221

            )

#970508
col_names_970 = ['Frequency(Hz)','Days','Flux(erg cm^-2','Flux uncertainty']

GRB970508 = 'radio970508.dat'   


redshift_970 = 0.835

lum_dist_970 = 5358.8  * 3.08568e20  # Mpc to cm
df_970 = pd.read_csv(GRB970508, sep='\s+', header = None, names = col_names_970)

# negative value filtering
df_970 = df_970[df_970['Flux(erg cm^-2'] >= 0]
df_970 = df_970[df_970['Flux uncertainty'] >= 0]

# micro Janskys to ergs
df_970['Flux(erg cm^-2'] = df_970['Flux(erg cm^-2'] * 10e-29
df_970['Flux uncertainty'] = df_970['Flux uncertainty'] * 10e-29


freq_13_970, freq_49_970, freq_85_970, freq_150_970, flux_total, luminosity = radio_bands(df_970,
                            upper_13 = 1.41e9*(1 + redshift_970),
                            lower_13 = 1.25e9*(1 + redshift_970),
                            upper_49 = 5.1e9*(1 + redshift_970),
                            lower_49 = 4.75e9*(1 + redshift_970), 
                            upper_85 = 8.6e9*(1+ redshift_970),
                            lower_85 =  8.4e9 * (1 + redshift_970),
                            upper_150 = 1.6e10*(1 + redshift_970),
                            lower_150 = 1.45e10*(1 + redshift_970),
                            redshift= redshift_970,
                            lum_dist= lum_dist_970             
                            )                                                                      
                                                                                                  

#130427A
col_names_130 = ['Frequency(Hz)','Days','Flux(erg cm^-2','Flux uncertainty']

GRB130427A = 'broadbandalldata.dat'

redshift_130 = 0.34


df_130 = pd.read_csv(GRB130427A, sep='\s+', header = None, names = col_names_130)

lum_dist_130 = 1811.7  * 3.08568e20  # Mpc to cm   #FIX

# mJy to erg
df_130['Flux(erg cm^-2'] = df_130['Flux(erg cm^-2'] * 10e-26
df_130['Flux uncertainty'] = df_130['Flux uncertainty'] * 10e-26


freq_13_130, freq_49_130, freq_85_130, freq_150_130, flux_total, luminosity = radio_bands(df_130,
                       upper_13 = 1.41e9*(1 + redshift_130),
                       lower_13 = 1.25e9*(1 + redshift_130),
                       upper_49 = 5.2e9*(1 + redshift_130),
                       lower_49 = 4.75e9*(1 + redshift_130),
                       upper_85 = 8.6e9*(1 + redshift_130),
                       lower_85 = 8.4e9*(1 + redshift_130),
                       upper_150 = 1.61e10*(1 + redshift_130),
                       lower_150 = 1.39e10*(1 + redshift_130),
                        redshift=redshift_130,
                        lum_dist=lum_dist_130       )
                            
                               





#light_curves('221009A', freq_150, freq_85, freq_49, freq_13)
#spectral_index_plots('221009A', spix_85_150, spix_49_85, spix_13_49)

fig , ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    nrows = 2, 
    ncols = 2, 
    figsize = (10, 10))

ax1.errorbar(freq_150_030['Days'], freq_150_030['Flux(erg cm^-2'], 
             freq_150_030['Flux uncertainty'], 
             fmt = 'o', 
             c ='c', 
             label = '030329')

ax1.set_title('comparing Flux- 8.5 GHz')
ax1.set_xlabel('Days')
ax1.set_ylabel('Flux(erg cm^-2')
ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.errorbar(freq_150_221['Days'], freq_150_221['Flux(erg cm^-2'], 
             freq_150_221['Flux uncertainty'], 
             fmt = 'o', 
             c ='b', 
             label = '221009A')
ax1.set_xlabel('Days')
ax1.set_ylabel('Flux(erg cm^-2')
ax1.set_xscale('log')
ax1.set_yscale('log')


ax1.errorbar(freq_150_130['Days'], freq_150_130['Flux(erg cm^-2'], 
             freq_150_130['Flux uncertainty'], 
             fmt = 'o', 
             c ='g', 
             label = '130427A')
ax1.set_xlabel('Days')
ax1.set_ylabel('Flux(erg cm^-2')
ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.legend() 

ax2.errorbar(freq_85_030['Days'], freq_85_030['Flux(erg cm^-2'], 
             freq_85_030['Flux uncertainty'], 
             fmt = 'o', 
             c ='c', 
             label = '030329')

ax2.set_title('comparing Flux- 8.5 GHz')
ax2.set_xlabel('Days')
ax2.set_ylabel('Flux(erg cm^-2')
ax2.set_xscale('log')
ax2.set_yscale('log')

          

ax2.errorbar(freq_85_221['Days'], freq_85_221['Flux(erg cm^-2'], 
             freq_85_221['Flux uncertainty'], 
             fmt = 'o', 
             c ='b', 
             label = '221009A')


ax2.set_xlabel('Days')
ax2.set_ylabel('Flux(erg cm^-2')
ax2.set_xscale('log')
ax2.set_yscale('log')


ax2.errorbar(freq_85_970['Days'], freq_85_970['Flux(erg cm^-2'], 
             freq_85_970['Flux uncertainty'], 
             fmt = 'o', 
             c ='r', 
             label = '970508')


ax2.set_xlabel('Days')
ax2.set_ylabel('Flux(erg cm^-2')
ax2.set_xscale('log')
ax2.set_yscale('log')


ax2.errorbar(freq_85_130['Days'], freq_85_130['Flux(erg cm^-2'], 
             freq_85_130['Flux uncertainty'], 
             fmt = 'o', 
             c ='g', 
             label = '130427A')


ax2.set_xlabel('Days')
ax2.set_ylabel('Flux(erg cm^-2')
ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.legend() 

ax3.errorbar(freq_49_030['Days'], freq_49_030['Flux(erg cm^-2'], 
             freq_49_030['Flux uncertainty'], 
             fmt = 'o', 
             c ='c', 
             label = '030329')
ax3.set_title('comparing Flux- 4.9 GHz')
ax3.set_xlabel('Days')
ax3.set_ylabel('Flux(erg cm^-2')
ax3.set_xscale('log')
ax3.set_yscale('log')

          

ax3.errorbar(freq_49_221['Days'], freq_49_221['Flux(erg cm^-2'], 
             freq_49_221['Flux uncertainty'], 
             fmt = 'o', 
             c ='b', 
             label = '221009A')
ax3.set_xlabel('Days')
ax3.set_ylabel('Flux(erg cm^-2')
ax3.set_xscale('log')
ax3.set_yscale('log')

ax3.errorbar(freq_49_970['Days'], freq_49_970['Flux(erg cm^-2'], 
             freq_49_970['Flux uncertainty'], 
             fmt = 'o', 
             c ='r', 
             label = '970508')
ax3.set_xlabel('Days')
ax3.set_ylabel('Flux(erg cm^-2')
ax3.set_xscale('log')
ax3.set_yscale('log')

ax3.errorbar(freq_49_130['Days'], freq_49_130['Flux(erg cm^-2'], 
             freq_49_130['Flux uncertainty'], 
             fmt = 'o', 
             c ='g', 
             label = '130427A')
ax3.set_xlabel('Days')
ax3.set_ylabel('Flux(erg cm^-2')
ax3.set_xscale('log')
ax3.set_yscale('log')

ax3.legend() 


ax4.errorbar(freq_13_030['Days'], freq_13_030['Flux(erg cm^-2'], 
             freq_13_030['Flux uncertainty'], 
             fmt = 'o', 
             c ='c', 
             label = '030329')
ax4.set_title('comparing Flux- 1.3 GHz')
ax4.set_xlabel('Days')
ax4.set_ylabel('Flux(erg cm^-2')
ax4.set_xscale('log')
ax4.set_yscale('log')

          

ax4.errorbar(freq_13_221['Days'], freq_13_221['Flux(erg cm^-2'], 
             freq_13_221['Flux uncertainty'], 
             fmt = 'o', 
             c ='b', 
             label = '221009A')


ax4.errorbar(freq_13_970['Days'], freq_13_970['Flux(erg cm^-2'], 
             freq_13_970['Flux uncertainty'], 
             fmt = 'o', 
             c ='r', 
             label = '970508')


ax4.errorbar(freq_13_130['Days'], freq_13_130['Flux(erg cm^-2'], 
             freq_13_130['Flux uncertainty'], 
             fmt = 'o', 
             c ='g', 
             label = '130427A')


ax4.legend() 


#spectral index comp
'''
upper_13 = 1.41e9*(1 + redshift_030),
lower_13 = 1.25e9*(1 + redshift_030),
upper_49 = 5.1e9*(1 + redshift_030),
lower_49 = 4.75e9*(1 + redshift_030),
upper_85 = 8.7e9*(1 + redshift_030),
lower_85 = 8.4e9*(1 + redshift_030),
upper_150 = 1.6e10*(1 + redshift_030),
lower_150 = 1.45e10*(1 + redshift_030),
'''
spix_85_150_030 = spectral_index(df_030, 
                                 lower = 8.4e9*(1 + redshift_030), 
                                 upper = 1.6e10*(1 + redshift_030))    
spix_49_85_030 = spectral_index(df_030, 
                                lower = 4.75e9*(1 + redshift_030), 
                                upper = 8.7e9*(1 + redshift_030))
spix_13_49_030 = spectral_index(df_030, 
                                lower = 1.25e9*(1 + redshift_030), 
                                upper = 5.1e9*(1 + redshift_030))

'''
upper_13 = 1.41e9*(1 + redshift_130),
lower_13 = 1.25e9*(1 + redshift_130),
upper_49 = 5.2e9*(1 + redshift_130),
lower_49 = 4.75e9*(1 + redshift_130),
upper_85 = 8.6e9*(1 + redshift_130),
lower_85 = 8.4e9*(1 + redshift_130),
upper_150 = 1.61e10*(1 + redshift_130),
lower_150 = 1.39e10*(1 + redshift_130),
'''
spix_85_150_130 = spectral_index(df_130, 
                                 lower = 8.4e9*(1 + redshift_130),
                                 upper = 1.61e10*(1 + redshift_130))    
spix_49_85_130 = spectral_index(df_030, 
                                lower = 4.75e9*(1 + redshift_130),
                                upper = 8.6e9*(1 + redshift_130))
spix_13_49_130 = spectral_index(df_030, 
                                lower = 1.25e9*(1 + redshift_130),
                                upper = 5.2e9*(1 + redshift_130))

'''
upper_13 = 1.41e9*(1 + redshift_221),
lower_13 = 1.25e9*(1 + redshift_221),
upper_49 = 5.2e9*(1 + redshift_221),
lower_49 = 4.75e9*(1 + redshift_221),
upper_85 = 8.6e9*(1 + redshift_221),
lower_85 = 8.4e9*(1 + redshift_221),
upper_150 = 1.61e10*(1 + redshift_221),
lower_150 = 1.39e10*(1 + redshift_221),
'''
spix_85_150_221 = spectral_index(df_130, 
                                 lower = 8.4e9*(1 + redshift_221),
                                 upper = 1.61e10*(1 + redshift_221))  
spix_49_85_221 = spectral_index(df_030, 
                                lower = 4.75e9*(1 + redshift_221),
                                upper = 8.6e9*(1 + redshift_221))
spix_13_49_221 = spectral_index(df_030, 
                                lower = 1.25e9*(1 + redshift_221),
                                upper = 5.2e9*(1 + redshift_221))

fig , ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    nrows = 2, 
    ncols = 2, 
    figsize = (10, 10))

ax1.errorbar(x= spix_85_150_030['bin center'], 
             y=spix_85_150_030['alpha'],
             yerr = spix_85_150_030['alpha err'],
             fmt = 'o',
             c = 'c',
             label = '030329'
             
             )
ax1.set_title("spectral index 8.5-15 GHz")
ax1.set_xscale('log')
ax1.set_xlabel('Days')
ax1.set_ylabel('Spectral index')
ax1.axhline(y=2, linestyle='--')
ax1.axhline(y=(1/3), linestyle='--')
ax1.axhline(y=-0.6, linestyle='--')

ax1.errorbar(x= spix_85_150_130['bin center'], 
             y=spix_85_150_130['alpha'],
             yerr = spix_85_150_130['alpha err'],
             fmt = 'o',
             c = 'g',
             label = '030329'
             )


ax1.errorbar(x= spix_85_150_221['bin center'], 
             y=spix_85_150_221['alpha'],
             yerr = spix_85_150_221['alpha err'],
             fmt = 'o',
             c = 'b',
             label = '221009A'
             )
ax1.legend()








luminosity_func(df_all=df_030, lum_dist=lum_dist_030, redshift=redshift_030)
luminosity_func(df_all=df_130, lum_dist=lum_dist_130, redshift=redshift_130)
luminosity_func(df_all=df_221, lum_dist=lum_dist_221, redshift=redshift_221)
luminosity_func(df_all=df_970, lum_dist=lum_dist_970, redshift=redshift_970)


fig , (ax1) = plt.subplots(
    nrows = 1, 
    ncols = 1, 
    figsize = (10, 10))
ax1.scatter(df_130['Days'], df_130['Luminosity'],
            c='g',
            label= '130427A')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title('GRB Spectral Luminosity over time')

ax1.scatter(df_030['Days'], df_030['Luminosity'],
            c='c',
            label= '030329')



ax1.scatter(df_221['Days'], df_221['Luminosity'],
            c='b',
            label= '221009A')


ax1.scatter(df_970['Days'], df_970['Luminosity'],
            c='r',
            label= '970508')

ax1.legend() 



