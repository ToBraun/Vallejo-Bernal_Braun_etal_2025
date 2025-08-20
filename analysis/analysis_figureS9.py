# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes the statistics necessary to plot Fig. S9 of 
Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python analysis_figureS9.py
'''

# %% IMPORT MODULES

# Standard library imports
import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr

# %% PARAMETERS

# Parameters
years = range(1940, 2024)

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Folders and files
output_base_name = config['output_base_name']
eulerian_file = os.path.join(
    config['output_folder'], 'PIKARTV1_eulerian', 
    output_base_name %('PIKARTV1_eulerian', '%s', 'nc'))

lagrangian_file = '1940-2023_PIKARTV1_Vallejo-Bernal_Braun_etal_2025.csv'

output_folder = './output/figureS9/'
os.makedirs(output_folder, exist_ok=True)

# %% FREQUENCY OF AR CONDITIONS

# Load the years
partial_sum = []
size_mask = 0

for year in years:
    
    # Open dataset
    eulerian_pikart = xr.open_dataset(eulerian_file %year)
    
    # Load AR mask
    mask_year = eulerian_pikart['ar_mask'].isel().load()
    
    # Calculate the total number of time steps with AR conditions
    partial_sum.append((mask_year > 0).sum(dim='time').values)
    
    # Calculate the total number of time steps in the year
    size_mask += mask_year.sizes['time']
    
    # Close the dataset
    eulerian_pikart.close()
    print(year)

# Concatenate the years along the time axis
sum_years = np.stack(partial_sum, axis=0)

# Calculate AR frequency
ar_freq = 100*np.sum(sum_years, axis=0)/size_mask
print('Total number of time steps: %s' %size_mask)

# Load latitude and longitude
eulerian_dataset = xr.open_dataset(eulerian_file %years[0])  
lat = eulerian_dataset['latitude'].values  # extract/copy the data
lon = eulerian_dataset['longitude'].values

# Save results
np.save(output_folder + '%s-%s_ar_freq.npy' %(years[0], years[-1]), ar_freq)
np.save(output_folder + 'lon_pikart.npy', lon)
np.save(output_folder + 'lat_pikart.npy', lat)

# %% LAGRANGIAN LANDFALLS WILL BE BUBBLES

# Load the Lagrangian catalog:
d_ars = pd.read_csv(lagrangian_file)
d_ars['time'] = pd.to_datetime(d_ars['time'])

# Subset to 1979-2023
d_ars = d_ars[(d_ars['time'] >= '%s-01-01 00:00:00' %years[0]) & (d_ars['time'] < '%s-01-01 00:00:00' %(years[-1]+1))]

# We normalise the bubbles by the total number of unique time instances to get a percentage:
Tnorm = np.unique(d_ars['time'].values).size

# Step 1: Filter out masked values
df_valid = d_ars[(~np.isnan(d_ars['lf_lat'])) & (~np.isnan(d_ars['lf_lon']))]

# Step 2: Get unique latitudes and longitudes
unique_lats = np.sort(df_valid['lf_lat'].unique())
unique_lons = np.sort(df_valid['lf_lon'].unique())

# Create mappings from lat/lon to matrix indices
lat_to_index = {lat: i for i, lat in enumerate(unique_lats)}
lon_to_index = {lon: i for i, lon in enumerate(unique_lons)}

# Step 3: Define landfall histogram matrix
n_lat = len(unique_lats)
n_lon = len(unique_lons)
a_lf_matrix = np.zeros((n_lat, n_lon), dtype=int)

# Step 4: Count landfalls: how frequently is a grid cell subjected to landfall?
for _, row in df_valid.iterrows():
    lat_idx = lat_to_index[row['lf_lat']]
    lon_idx = lon_to_index[row['lf_lon']]
    a_lf_matrix[lat_idx, lon_idx] += 1
    
# %% SAVE OUTPUTS

np.save(output_folder + '%s-%s_unique_lons.npy' %(years[0], years[-1]), unique_lons)
np.save(output_folder + '%s-%s_unique_lats.npy' %(years[0], years[-1]), unique_lats)
np.save(output_folder + '%s-%s_lf_map.npy' %(years[0], years[-1]), a_lf_matrix)