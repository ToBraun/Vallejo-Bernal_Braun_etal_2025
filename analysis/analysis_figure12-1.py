# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes time series necessary for the trend analysis in 
Vallejo-Bernal & Braun et al., (2025).

We do the calculations for four time chunks by defining year_start and year_end. 

Run this script from the command line with:

    python analysis_figure12(1).py
'''

# %% IMPORT MODULES

import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr
from analysis_functions import count_ARevents

# %% PARAMETERS

# Years 
year_start = 1940
year_end = 1959
Nyears = year_end - year_start + 1

# Print log
print('ANALYSIS FIGURE 12: %s-%s\n' %(year_start, year_end))

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Parameters
lon_var = config['lon_var']
lat_var = config['lat_var']
var = config['var']
data_base_name = config['data_base_name']
output_base_name = config['output_base_name']
W = 4 # spatial window that results in 2deg resolution for final maps

# Folders and files
output_base_name = config['output_base_name']
ivt_file = os.path.join(config['ivt_folder'], data_base_name %(var, '%s'))
eulerian_file = os.path.join(
    config['output_folder'], 'PIKARTV1_eulerian', 
    output_base_name %('PIKARTV1_eulerian', '%s', 'nc'))

output_folder = './output/figure12/'
os.makedirs(output_folder, exist_ok=True)

# %% LOAD DATA

l_ivt, l_mask, l_time_ivt, l_time_mask = [], [], [], []

# Loop through years
for year in range(year_start, year_end + 1):
    
    # Open IVT dataset
    with xr.open_dataset(ivt_file %year) as ds:
        
        # Load IVT
        ivt_year = ds[var].load()
        
        # Load time
        time_ivt_year = ds['time'].load()
        
        # Append the data of the year
        l_ivt.append(ivt_year)
        l_time_ivt.append(time_ivt_year)
         
    # Open eulerian dataset
    with xr.open_dataset(eulerian_file %year) as ds:
        
        # Load AR mask
        mask_year = ds['ar_mask'].load()
        
        # Load time
        time_mask_year = ds['time'].load()
        
        # Append the data of the year
        l_mask.append(mask_year)
        l_time_mask.append(time_mask_year)

# Concatenate the data along the time axis
a_ivt = np.concatenate(l_ivt, axis=0)
a_mask = np.concatenate(l_mask, axis=0)
time_ivt = pd.to_datetime(np.concatenate(l_time_ivt, axis=0))
time_mask = pd.to_datetime(np.concatenate(l_time_mask, axis=0))

# Load coordinates
eulerian_dataset = xr.open_dataset(eulerian_file %year_start)  
a_lats = eulerian_dataset[lat_var].values  # extract/copy the data
a_lons = eulerian_dataset[lon_var].values
eulerian_dataset.close()
nlat, nlon = a_lats.size, a_lons.size

# Make sure that the IVT and the AR mask share the same time axis
a_ivt = a_ivt[np.isin(time_ivt, time_mask),:,:]

# %% CALCULATIONS

##### Spatial aggregation to 2deeg x 2deeg grid

# Coordinates
aggr_lons = np.add.reduceat(a_lons, np.arange(0, nlon, W), axis=0)/W
aggr_lats = np.add.reduceat(a_lats, np.arange(0, nlat, W), axis=0)/W

# AR mask
aggr_mask = np.add.reduceat(
    np.add.reduceat(a_mask, np.arange(0, nlat, W), axis=1, dtype=int),
    np.arange(0, nlon, W), axis=2, dtype=int)
aggr_mask = np.where(aggr_mask > 0, 1, 0)
naggrlat = np.size(aggr_mask,1)
naggrlon = np.size(aggr_mask,2)

# Initialize variables
a_ann_counts, a_annavg_dur = np.zeros(
    (Nyears, naggrlat, naggrlon)), np.zeros((Nyears, naggrlat, naggrlon))
a_annivt_AR, a_annivt_NOAR = np.zeros(
    (Nyears, nlat, nlon)), np.zeros((Nyears, nlat, nlon))

# Loop through years
for i in range(Nyears):
    
    ##### Count and duration: AR events
    
    # Select a year 
    tmp_idx = (time_mask.year == (year_start + i))
    tmp_mask = aggr_mask[tmp_idx,:,:]
    
    # Compute the annual number of AR events and their average duration
    a_ann_counts[i,:,:], a_annavg_dur[i,:,:] = count_ARevents(tmp_mask)
    
    ##### IVT: AR and non-AR conditions
    
    # Split into AR and non-AR related IVT
    tmp_ivt_AR = a_ivt[tmp_idx,:,:] * a_mask[tmp_idx,:,:]
    tmp_ivt_NONAR = a_ivt[tmp_idx,:,:] * (1-a_mask[tmp_idx,:,:])
    
    # We only wanna average over the respective non-zero pixels, so let's mask 
    # the others with NANs
    tmp_ivt_AR = np.where(tmp_ivt_AR > 0, tmp_ivt_AR, np.nan)
    tmp_ivt_NONAR = np.where(tmp_ivt_NONAR > 0, tmp_ivt_NONAR, np.nan)
    
    # Temporal averaging
    a_annivt_AR[i,:,:] = np.nanmean(tmp_ivt_AR, 0)
    a_annivt_NOAR[i,:,:] = np.nanmean(tmp_ivt_NONAR, 0)
    
    print(year_start + i)
    
##### Spatial aggregation to 2deeg x 2deeg grid

# AR-related IVT
# Sum first
a_totalivt_AR = np.add.reduceat(
    np.add.reduceat(a_annivt_AR, np.arange(0, nlat, W), axis=1), 
    np.arange(0, nlon, W), axis=2)

# Only count AR-pixels to exclude effects due to AR vs non-AR frequency 
tmp_nonzero = np.where(a_annivt_AR > 0, 1, 0)
a_nonzero_AR = np.add.reduceat(
    np.add.reduceat(tmp_nonzero, np.arange(0, nlat, W), axis=1),
    np.arange(0, nlon, W), axis=2)
a_meanivt_AR = a_totalivt_AR/a_nonzero_AR

# NON-AR-related IVT
# Sum first
a_totalivt_NOAR = np.add.reduceat(
    np.add.reduceat(a_annivt_NOAR, np.arange(0, nlat, W), axis=1),
    np.arange(0, nlon, W), axis=2)

# Only count NON-AR-pixels to exclude effects due to AR vs non-AR frequency 
tmp_nonzero = np.where(a_annivt_NOAR > 0, 1, 0)
a_nonzero_NOAR = np.add.reduceat(
    np.add.reduceat(tmp_nonzero, np.arange(0, nlat, W), axis=1),
    np.arange(0, nlon, W), axis=2)
a_meanivt_NOAR = a_totalivt_NOAR/a_nonzero_NOAR

# %% SAVE RESULTS

# Write those big boys out
np.save(output_folder + 
        '%s-%s_ar_events.npy' %(year_start, year_end), a_ann_counts)
np.save(output_folder + 
        '%s-%s_mean_dur.npy' %(year_start, year_end), a_annavg_dur)
np.save(output_folder + 
        '%s-%s_meanivt_AR.npy' %(year_start, year_end), a_meanivt_AR)
np.save(output_folder + 
        '%s-%s_meanivt_NOAR.npy' %(year_start, year_end), a_meanivt_NOAR)
np.save(output_folder + 'lon_pikart.npy', a_lons)
np.save(output_folder + 'lat_pikart.npy', a_lats)