# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes the statistics necessary to plot Fig. S7 of 
Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python analysis_figureS7.py
'''

# %% IMPORT MODULES

import os
import yaml
import numpy as np
import xarray as xr

# %% PARAMETERS

# Define years to process
years = np.arange(1980, 2021) 

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Folders and files
output_base_name = config['output_base_name']
eulerian_file = os.path.join(
    config['output_folder'], 'PIKARTV1_eulerian', 
    output_base_name %('PIKARTV1_eulerian', '%s', 'nc'))

path = ...
target_file = path + '/target_v4_ERA5_1940-2024.nc'

output_folder = './output/figureS7/'
os.makedirs(output_folder, exist_ok=True)

# %% FREQUENCY OF AR CONDITIONS IN PIKART

# Load the years
partial_sum = []
size_mask = 0

for year in years:
    
    # Open dataset
    with xr.open_dataset(eulerian_file %year) as ds:
    
        # Load AR mask
        mask_year = ds['ar_mask'].load()
        
        # Calculate the total number of time steps with AR conditions
        partial_sum.append((mask_year > 0).sum(dim='time').values)
        
        # Calculate the total number of time steps in the year
        size_mask += mask_year.sizes['time']
    
    print(year)

# Concatenate the years along the time axis
sum_years = np.stack(partial_sum, axis=0)

# Calculate AR frequency
ar_freq_pikart = 100*np.sum(sum_years, axis=0)/size_mask
print('Total number of time steps PIKART: %s' %size_mask)

# Load latitude and longitude
eulerian_dataset = xr.open_dataset(eulerian_file %years[0])  
lat_pikart = eulerian_dataset['latitude'].values  # extract/copy the data
lon_pikart = eulerian_dataset['longitude'].values
eulerian_dataset.close()

np.save(output_folder + '%s-%s_ar_freq_pikart.npy' %(years[0], years[-1]), 
        ar_freq_pikart)
np.save(output_folder + 'lon_pikart.npy', lon_pikart)
np.save(output_folder + 'lat_pikart.npy', lat_pikart)

# %% FREQUENCY OF AR CONDITIONS IN TARGET

# Load dataset
target = xr.open_dataset(target_file)

# Extract coordinate variables
lons = target['lon'].values
lats = target['lat'].values

# Initialize array to accumulate AR days per grid cell
ar_mask_target = np.zeros((len(lats), len(lons)))

# Loop over years
size_mask = 0
for year in years:
    print(f"Processing year {year}...")
    
    # Select data for the given year
    target_year = target.sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
    ar_mask = (np.squeeze(target_year["shapemap"].values) > 0).astype(int)
    ar_mask_target += np.sum(ar_mask, axis=0)
    size_mask += np.size(ar_mask,0)
    target_year = None
    
# Frequency of AR conditions
ar_freq_target = 100*ar_mask_target/size_mask
print('Total number of time steps tARget: %s' %size_mask)

# Save results
np.save(output_folder + '%s-%s_ar_freq_target.npy' %(years[0], years[-1]), 
        ar_freq_target)
np.save(output_folder + 'lon_target.npy', lons)
np.save(output_folder + 'lat_target.npy', lats)
