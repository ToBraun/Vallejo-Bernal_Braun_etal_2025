# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script merges the longitudinal chunks and generates the input data
for Figs. 9 and 10 of Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python merge_figures9-10.py
'''

# %% IMPORT MODULES

import os
import yaml
import numpy as np
import xarray as xr
import time as tm
import datetime as dt

# %% PARAMETERS

# Starting time
start = tm.time()
print('Starting time: %s' %dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print('FIGURE 9 MERGER\n')

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Parameters
lon_var = config['lon_var']
lat_var = config['lat_var']
output_base_name = config['output_base_name']

# Folders and files
input_folder = './output/'
eulerian_file = os.path.join(
    config['output_folder'], 'PIKARTV1_eulerian', 
    output_base_name %('PIKARTV1_eulerian', '1979', 'nc'))

# %% MERGE RESULTS FOR FIGURE 9

# Load the first Eulerian file
nc_ecat = xr.open_dataset(eulerian_file)

# Load coordinates
lat = nc_ecat[lat_var].values
lon = nc_ecat[lon_var].values

# Load the first chunk
cont_rank_count0 = np.load(input_folder + 'figure09/cont_rank_count_0.npy')

# Initialize variables
cont_rank_count = np.zeros((np.size(cont_rank_count0,0), np.size(cont_rank_count0,1), len(lat), len(lon)), dtype=int)

# Load the chunks
for i in range(36):
    chunk = [20*i, 20*(i+1)]
    cont_rank_count[:,:,:,chunk[0]:chunk[1]] = np.load(input_folder + 'figure09/cont_rank_count_%s.npy' %i)
    
assert chunk[1] == len(lon), 'Incomplete longitudinal chunks'

# Save outputs
np.save(input_folder + 'figure09/cont_rank_count.npy', cont_rank_count)

# %% MERGE RESULTS FOR FIGURE 10

# Load the first chunk
rank_count0 = np.load(input_folder + 'figure10/rank_count_0.npy')

# Initialize variables
rank_count = np.zeros((np.size(rank_count0,0), len(lat), len(lon)), dtype=int)

# Load the chunks
for i in range(36):
    chunk = [20*i, 20*(i+1)]
    rank_count[:,:,chunk[0]:chunk[1]] = np.load(input_folder + 'figure10/rank_count_%s.npy' %i)
    
assert chunk[1] == len(lon), 'Incomplete longitudinal chunks'

# Save outputs
np.save(input_folder + 'figure10/rank_count.npy', rank_count)
np.save(input_folder + 'figure10/lon_pikart.npy', lon)
np.save(input_folder + 'figure10/lat_pikart.npy', lat)

# Ending time
end = tm.time()

# Print processing time
print('Ending time: %s' %dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print('Processing time: %s hours' %round((end-start)/3600, 2))
