# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes the actual trends necessary to plot Fig. 12 of 
Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python analysis_figure12(2).py
'''

# %% IMPORT MODULES

import numpy as np
import scipy
from scipy.stats import t
import rioxarray
from analysis_functions import SSA

root_folder = ...

# %% INITIAL SETTINGS

# Parameters
central_lon = -150
alpha = .05 # significance level
res_aggr_grid = 2 # aggregated resolution
W = 4 # spatial aggregation window
lat_thresh = 23 # threshold for the spatila separation into NH, SH, and tropics
# CC-scaling (estimates obtained from NASA Goddard Institute for Space Studies (GISS)):
Tlow, Thigh = 0.9, 1.0 # overall warming estimates (lower, upper)
Tclow, Tchigh = 1.2, 1.5 # continental warming estimates (lower, upper)
Tolow, Tohigh = 0.5, 0.8 # oceanic warming estimates (lower, upper)

# Folders and files
continents_folder = ...
input_folder = root_folder + 'data/figure12/'

# %% DATA
l_y0, l_y1 = [1940,1960,1980,2000], [1959,1979,1999,2023]
L = len(l_y0)

## NEW DATA
l_counts, l_dur, l_ARivt, l_ivt = [],[],[],[]
for n in range(L):
    l_counts.append(
        np.load(input_folder + '%s-%s_ar_events.npy' %(l_y0[n], l_y1[n])))
    l_dur.append(
        np.load(input_folder + '%s-%s_mean_dur.npy' %(l_y0[n], l_y1[n])))
    l_ARivt.append(
        np.load(input_folder + '%s-%s_meanivt_AR.npy' %(l_y0[n], l_y1[n])))
    l_ivt.append(
        np.load(input_folder + '%s-%s_meanivt_NOAR.npy' %(l_y0[n], l_y1[n])))
# Convert to arrays
a_counts, a_dur = np.vstack(l_counts), np.vstack(l_dur)
a_ARivt, a_ivt = np.vstack(l_ARivt), np.vstack(l_ivt)
l_allarrs = [a_counts, a_dur, a_ARivt, a_ivt]

## AGGREGATION
Tyears, Lla, Llo = a_counts.shape
a_years = np.arange(l_y0[0], l_y1[-1]+1)
# Grid coordinates
longitude = np.load(input_folder + 'lon_pikart.npy')
latitude = np.load(input_folder + 'lat_pikart.npy')
# Aggregate
lon_starts = np.arange(longitude[0], longitude[-1], res_aggr_grid)
lat_starts = np.arange(latitude[0], latitude[-1], res_aggr_grid)
Llo, Lla = lon_starts.size, lat_starts.size
aggr_lons_cen, aggr_lats_cen = lon_starts + res_aggr_grid/2, lat_starts + res_aggr_grid/2

# %% COMPUTE trends and significance

## LINEAR TRENDS
# 4 in order of the figure's panels
a_slopes, a_signif = np.zeros((4, Lla, Llo)), np.zeros((4, Lla, Llo)) 

for n in range(4):
    
    a_arr = l_allarrs[n]
    
    for i in range(Llo):
        for j in range(Lla):  
            
            tmp_cell = a_arr[:,j,i]
            
            # Find NAN or zero values: Perform fit over non-NAN/non-zero years 
            # only. Discard if more than half of the years are NAN/zero.
            tmp_nonan = np.where(~np.isnan(tmp_cell) * (tmp_cell > 0))[0]
            Nval = tmp_nonan.size
            
            if Nval > Tyears/2:
                slope, intercept, r, p, se = scipy.stats.linregress(
                    a_years[tmp_nonan], tmp_cell[tmp_nonan])
            else:
                slope, intercept, r, p, se = 0,0,0,0,0
            
            # Save slope
            a_slopes[n, j,i] = slope
            
            # Significance testing
            tinv = lambda p, df: abs(t.ppf(p/2, df)); ts = tinv(
                alpha, len(a_years)-2)
            if (slope + ts*se > 0) & (slope - ts*se > 0) & (slope != 0):
                a_signif[n, j,i] = 1
            elif (slope + ts*se < 0) & (slope - ts*se < 0) & (slope != 0):
                a_signif[n, j,i] = 1
            elif slope == 0:         
                a_signif[n, j,i] = np.nan
            else:         
                a_signif[n, j,i] = np.nan
    
## Comparison: AR vs non-AR related trends
# Show slopes only if AR-related IVT is higher than non-AR-related IVT
a_slopes[3,] = np.where(((a_slopes[2,] > 0) & (a_slopes[2,] > a_slopes[3,]) & 
                         (~np.isnan(a_slopes[3,]))),  a_slopes[3,], np.nan)
a_signif[3,] = np.where(((a_slopes[2,] > 0) & (a_slopes[2,] > a_slopes[3,]) & 
                         (~np.isnan(a_slopes[3,]))), a_signif[3,], np.nan)
    
# CLAUSIUS-CLAPEYRON ESTIMATION
a_slopes_CC = np.nan*np.ones((2, 3, 2, Lla, Llo))

for n in range(2):
    
    a_arr = l_allarrs[n+2]
    
    for i in range(Llo):
        for j in range(Lla):  
            
            tmp_cell = a_arr[:,j,i]
            
            # Find NAN-values: perform fit over non-NAN years only. 
            # Discard if more than half of the years are NAN.
            tmp_nonan = np.where(~np.isnan(tmp_cell))[0]
            
            if tmp_nonan.size > Tyears/2:
                
                # Clausius-Clapeyron
                
                ## CONTINENTAL
                tmp_x = np.linspace(0, Tclow, Tyears)
                a_slopes_CC[n, 0, 0, j,i] = 100*scipy.stats.linregress(
                    tmp_x[tmp_nonan], tmp_cell[tmp_nonan]).slope/tmp_cell[tmp_nonan][0]
                tmp_x = np.linspace(0, Tchigh, Tyears)
                a_slopes_CC[n, 0, 1, j,i] = 100*scipy.stats.linregress(
                    tmp_x[tmp_nonan], tmp_cell[tmp_nonan]).slope/tmp_cell[tmp_nonan][0]
                
                ## OCEANIC
                tmp_x = np.linspace(0, Tolow, Tyears)
                a_slopes_CC[n, 1, 0, j,i] = 100*scipy.stats.linregress(
                    tmp_x[tmp_nonan], tmp_cell[tmp_nonan]).slope/tmp_cell[tmp_nonan][0]
                tmp_x = np.linspace(0, Tohigh, Tyears)
                a_slopes_CC[n, 1, 1, j,i] = 100*scipy.stats.linregress(
                    tmp_x[tmp_nonan], tmp_cell[tmp_nonan]).slope/tmp_cell[tmp_nonan][0]
                
                ## OVERALL
                tmp_x = np.linspace(0, Tlow, Tyears)
                a_slopes_CC[n, 2, 0, j,i] = 100*scipy.stats.linregress(
                    tmp_x[tmp_nonan], tmp_cell[tmp_nonan]).slope/tmp_cell[tmp_nonan][0]
                tmp_x = np.linspace(0, Thigh, Tyears)
                a_slopes_CC[n, 2, 1, j,i] = 100*scipy.stats.linregress(
                    tmp_x[tmp_nonan], tmp_cell[tmp_nonan]).slope/tmp_cell[tmp_nonan][0]

# %% Clausius-Clapeyron scaling

# Here, I compute the CC estimates for continental and oceanic warming.
# Please insert them into the text.

# Load continents mask
a_conti = rioxarray.open_rasterio(continents_folder + 'world_continents_3000sqkm.tif')
a_conti = np.flip(np.squeeze(a_conti.values), axis=0)
a_mask = np.where(a_conti>0, 1, 0)
nlat, nlon = a_mask.shape
a_aggmask = np.add.reduceat(
    np.add.reduceat(a_mask, np.arange(0, nlat, W), axis=0, dtype=int), 
    np.arange(0, nlon, W), axis=1, dtype=int)
a_contimask = np.where(a_aggmask>0, 1, np.nan)
a_seamask = np.where(a_aggmask==0, 1, np.nan)

# MASK AND AVERAGE
## Dimensions: type x conti_ocean_overall x temperature estimate
print('## AR')
print('### OVERALL')
print(np.round(np.nanmean(a_slopes_CC[0, 2, 1, :, :]),2))
print(np.round(np.nanmean(a_slopes_CC[0, 2, 0, :, :]),2))
print('\n### CONTINENTAL')
print(np.round(np.nanmean(a_slopes_CC[0, 0, 1, :, :]*a_contimask),2))
print(np.round(np.nanmean(a_slopes_CC[0, 0, 0, :, :]*a_contimask),2))
print('\n### OCEANIC')
print(np.round(np.nanmean(a_slopes_CC[0, 1, 1, :, :]*a_seamask),2))
print(np.round(np.nanmean(a_slopes_CC[0, 1, 0, :, :]*a_seamask),2))

print('\n## NON-AR')
print('### OVERALL')
print(np.round(np.nanmean(a_slopes_CC[1, 2, 1, :, :]),2))
print(np.round(np.nanmean(a_slopes_CC[1, 2, 0, :, :]),2))
print('\n### CONTINENTAL')
print(np.round(np.nanmean(a_slopes_CC[1, 0, 1, :, :]*a_contimask),2))
print(np.round(np.nanmean(a_slopes_CC[1, 0, 0, :, :]*a_contimask),2))
print('\n### OCEANIC')
print(np.round(np.nanmean(a_slopes_CC[1, 1, 1, :, :]*a_seamask),2))
print(np.round(np.nanmean(a_slopes_CC[1, 1, 0, :, :]*a_seamask),2))

# %% NONLINEAR TRENDS

### Separation of three geographical zones by lat_thresh
a_nh, a_trop, a_sh, = np.zeros((4,Tyears)), np.zeros((4,Tyears)), np.zeros((4,Tyears))
a_nhtrend = np.zeros((4,Tyears))
a_troptrend = np.zeros((4,Tyears))
a_shtrend = np.zeros((4,Tyears))

for n in range(4):
    a_arr = l_allarrs[n]
    if n == 0:
        a_nh[n,] = np.nansum(np.nansum(a_arr[:,lat_starts>lat_thresh,:], 1), 1)
        a_trop[n,] = np.nansum(
            np.nansum(
                a_arr[:,np.where((lat_starts<=lat_thresh) & (lat_starts>=-lat_thresh)
                                 )[0],:], 1), 1)
        a_sh[n,] = np.nansum(np.nansum(a_arr[:,lat_starts<-lat_thresh,:], 1), 1)
    else:
        a_nh[n,] = np.nanmean(np.nanmean(a_arr[:,lat_starts>lat_thresh,:], 1), 1)
        a_trop[n,] = np.nanmean(
            np.nanmean(
                a_arr[:,np.where((lat_starts<=lat_thresh) & (lat_starts>=-lat_thresh)
                                 )[0],:], 1), 1)
        a_sh[n,] = np.nanmean(np.nanmean(a_arr[:,lat_starts<-lat_thresh,:], 1), 1)
    
    ### SSA
    a_nhtrend[n,] = SSA(a_nh[n,], 10).reconstruct(np.arange(1)).values
    a_troptrend[n,] = SSA(a_trop[n,], 10).reconstruct(np.arange(1)).values
    a_shtrend[n,] = SSA(a_sh[n,], 10).reconstruct(np.arange(1)).values
    
# %% SAVE OUTPUTS

np.save(input_folder + 'slopes.npy', a_slopes) 
np.save(input_folder + 'signif.npy', a_signif)
np.save(input_folder + 'aggr_lons.npy', aggr_lons_cen) 
np.save(input_folder + 'aggr_lats.npy', aggr_lats_cen)

np.save(input_folder + 'nh.npy', a_nh)
np.save(input_folder + 'trop.npy', a_trop)
np.save(input_folder + 'sh.npy', a_sh) 
np.save(input_folder + 'nh_trend.npy', a_nhtrend)
np.save(input_folder + 'trop_trend.npy', a_troptrend)
np.save(input_folder + 'sh_trend.npy', a_shtrend)
np.save(input_folder + 'years.npy', a_years)
