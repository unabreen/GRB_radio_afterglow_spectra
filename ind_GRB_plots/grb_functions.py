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



"""
This function filters out all flux values with an uncertainty value over
three standard deviations. It takes the names of the
main dataframe and the flux, flux err, and time columns.
"""
def uncertainty_clean(
        df_all,
        flux_col,
        flux_err_col
        
        ):
    
    #upper limit uncertainty
    sigma_three = df_all[flux_col].std()
    
    
    df_all = df_all[abs(df_all[flux_err_col]) <= sigma_three]
    
    
    
"""
This function removes negative and non numeric values. It takes the names of the
main dataframe and the flux, flux err, and time columns.
""" 
def clean_GRB_data(
        df_all,
        flux_col,
        flux_err_col,
        time_col
        ):
    # Keep only rows where all numeric columns have valid numbers
    numeric_cols = df_all.select_dtypes(include=['number']).columns
    df_all = df_all[df_all[numeric_cols].notna().all(axis=1)]
    df_all = df_all[df_all[flux_col] >= 0]
    df_all = df_all[df_all[flux_err_col] >= 0]
   
   

"""
This function takes specific dataset limits and filters measurements from
 the main dataset into four sub dataframes using boolean indexing. It takes the
 name of the main dataframe and its columns, and variables with the redshift 
 adjusted band limits as values.
 
 The values are corrected for redshift before filtering to match the redshift 
 corrected band limits. It returns the four sub dataframes.
"""
def radio_bands(
        df_all,
        upper_13,
        lower_13,
        upper_49,
        lower_49,
        upper_85,
        lower_85,
        upper_150,
        lower_150,
        lum_dist,
        redshift,
        time_col,
        freq_col,
        flux_col,
        flux_err_col
                ):

    

    #redshift corrections
    df_all[time_col] = df_all[time_col] / (1 + redshift)
    df_all[flux_col] = df_all[flux_col] / (1 + redshift)
    df_all[freq_col] = df_all[freq_col] *(1 + redshift)
    
    
    
    # dataframe catches adjusted for redshift
    #dataframe for GRB 030329 for 1.3 GHz (1.25-1.41)
    freq_13 = (df_all[(df_all[freq_col] <= upper_13) &
                         (df_all[freq_col] >= lower_13)])
    freq_13 = freq_13.sort_values(by=time_col)
                         
    #dataframe for GRB 030329 for 4.9 GHz (4.75-5.0)                   
    freq_49 = (df_all[(df_all[freq_col] <= upper_49) &
                         (df_all[freq_col] >= lower_49)])
    freq_49 = freq_49.sort_values(by=time_col)
    
    #dataframe for GRB 030329 for frequency = 8.5GHz(8.4-8.6)
    freq_85 = (df_all[(df_all[freq_col] <= upper_85 ) &
                         (df_all[freq_col] >= lower_85)])
    freq_85 = freq_85.sort_values(by=time_col)
    
    #dataframe for GRB 030329 for frequency = 15Hz
    freq_150 = (df_all[(df_all[freq_col] <= upper_150) &
                         (df_all[freq_col] >= lower_150)])
    freq_150 = freq_150.sort_values(by=time_col)
    
    #luminosity calculations
    #flux_total = np.trapezoid(df_all[freq_col], df_all[freq_col] )
    #luminosity = flux_total * 4 * np.pi * (lum_dist)**2
    
    return freq_13, freq_49, freq_85, freq_150, df_all


"""
This function takes the full dataframe along with upper and lower frequency
limits and bin width to create a spectral index dataframe between the two 
frequency limits. The spectral index is calculated as F~ nu^alpha. To calculate
the instantaneous value, the slope of the Flux vs frequency curve between 
the stated frequencies is measured with data points in the given bin. The 
function loops through the time bins, each time calculating the slope of the
points assigned to the bin using linear regression. If there are fewer than 
two points in each bin, the loop continues to the next bin. The slope values
 are stored in result_frame with each slope value having a corresponding
time bin center. The function returns the result frame.
"""
def spectral_index(
        df_full,
        lower,
        upper,
        freq_col,
        time_col,
        flux_col,
        bin_width = 1.0
        ):
    
    df = df_full[(df_full[freq_col] >= lower) & 
                 (df_full[freq_col] <= upper)]
    df = df.sort_values(by= time_col)
    
    t_min, t_max = df[time_col].min(), df[time_col].max()
    bins = np.arange(t_min, t_max, bin_width)
    
    
    results = []
    
    for i in range(len(bins)- 1):
        
        t_low, t_high = bins[i], bins[i+1]
        t_center = (t_low + t_high)/2

        in_bin = df[(df[time_col]>= t_low) & (df[time_col] <= t_high)]
        
        if len(np.unique(in_bin[freq_col])) < 2:
            continue
        
        x = np.log(in_bin[freq_col].values)
        y = np.log(in_bin[flux_col].values)

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
    
    return result_frame


"""
This function takes the sub dataframe returned from radio_bands() and plots
them on a 2x2 figure with seperate plots for each frequency band. If a sub
dataframe for a certain frequency does not exist, then the plot is not made and
comments that there is no data to plot in the title.
"""

def light_curves(
        grb_name,
        time_col,
        freq_col,
        flux_col,
        flux_err_col,
        Freq_150 = None,
        Freq_85 = None,
        Freq_49 = None,
        Freq_13 = None
        ):
    
    fig , ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        nrows = 2, 
        ncols = 2, 
        figsize = (10, 10))
        
    
    #15
    if Freq_150 is not None:
        ax1.errorbar(Freq_150[time_col], Freq_150[flux_col], 
                     yerr=Freq_150[flux_err_col], 
                     fmt = 'o', 
                     c ='b', 
                     label = '15 GHz')
        
        ax1.set_title(f'GRB {grb_name} Flux- 15 GHz')
        ax1.set_xlabel('Days')
        ax1.set_ylabel('Flux(mJy)')
        ax1.set_xscale('log')
        ax1.set_yscale('log')
    if Freq_150 is None:
        ax1.set_title("No data in 15 GHz")
    
    
    #8.5
    if Freq_85 is not None:
        ax2.errorbar(Freq_85[time_col], Freq_85[flux_col], 
                     yerr= Freq_85[flux_err_col], 
                     fmt = 'o', 
                     c ='c', 
                     label = '8.5 GHz')
        
        ax2.set_title(f'GRB {grb_name} Flux- 8.5 GHz')
        ax2.set_xlabel('Days')
        ax2.set_ylabel('Flux(mJy)')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
    if Freq_85 is None:
        ax2.set_title("No data in 8.5 GHz")
    
    #4.9
    if Freq_49 is not None:
        ax3.errorbar(Freq_49[time_col], Freq_49[flux_col], 
                     yerr= Freq_49[flux_err_col], 
                     fmt = 'o', 
                     c = 'r', 
                     label = '4.9 GHz')
        
        
        ax3.set_title(f'GRB {grb_name} Flux- 4.9 GHz')
        ax3.set_xlabel('Days')
        ax3.set_ylabel('Flux(mJy)')
        ax3.set_xscale('log')
        ax3.set_yscale('log')
    if Freq_49 is None:
        ax3.set_title("No data in 4.9 GHz")
    
    # 1.3
    if Freq_13 is not None:
        ax4.errorbar(Freq_13[time_col], Freq_13[flux_col], 
                     yerr = Freq_13[flux_err_col], 
                     fmt = 'o', 
                     c ='g', 
                     label = '1.3 GHz')
        
        ax4.set_title(f'GRB {grb_name} Flux- 1.3 GHz')
        ax4.set_xlabel('Days')
        ax4.set_ylabel('Flux(mJy)')
        ax4.set_xscale('log')
        ax4.set_yscale('log')
    if Freq_13 is None:
        ax4.set_title("No data in 1.3 GHz")


"""
This function takes the GRB name(number) as a string, and the spectral index
dataframes if they exist. The default value is set to None and will result in
an empty plot if not assigned. The spectral index values are plotted over time,
and values corresponding to the synchrotron spectrum are plotted as horizontal
asymptotes to more clearly visualize the spectral evolution.
"""

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
    elif freq_85_150 is None:
        ax1.set_title('No data in 8.5-15 GHz')
    
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
    elif freq_49_85 is None:
        ax2.set_title('No data in 4.9-8.5 GHz')
    
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
    elif freq_13_49 is None:
        ax3.set_title('No data in 1.3-4.9 GHz')





"""
This function calculates the luminosity corresponding to each flux measurement
and adds a column to the main dataframe with the new values. It uses luminosity
distance and redshift defined in a variable and the names of the dataframe and
flux column.
"""

def luminosity_func(df_all,
               lum_dist,
               redshift,
               flux_col):
    df_all["Luminosity"] = df_all[flux_col] * 4 * np.pi * (lum_dist/(1 + redshift))**2













