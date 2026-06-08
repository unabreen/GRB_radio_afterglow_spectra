#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 09:44:38 2026

@author: unabreen
"""
import pandas as pd
import numpy as np
from scipy.stats import linregress



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
    
    t_range = df[time_col].max() - df[time_col].min()
    
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