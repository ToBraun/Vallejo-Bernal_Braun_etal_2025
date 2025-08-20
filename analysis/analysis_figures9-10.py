# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes the statistics necessary to plot Figs. 9 and 10 of 
Vallejo-Bernal & Braun et al., (2025) for a longitudinal chunk of the data.

Run this script from the command line with:

    python analysis_figures9-10.py
'''

# %% IMPORT MODULES

import os
import sys
import yaml
import numpy as np
import xarray as xr
import rioxarray
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))
from ardt_functions import rle

# %% PARALLELIZATION

# Get task
task = int(os.environ['SLURM_ARRAY_TASK_ID'])

# Calculate chunk
chunk = [20*task, 20*(task+1)]
print('ANALYSIS FIGURES 9 AND 10: lon=%s\n' %chunk)

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Folders and files
output_base_name = config['output_base_name']
eulerian_file = os.path.join(
    config['output_folder'], 'PIKARTV1_eulerian', 
    output_base_name %('PIKARTV1_eulerian', '%s', 'nc'))

continents_folder = ...
output_folder = './output/'

# %% FREQUENCY OF AR RANKS

# Load the years
ar_ranks = []

for year in range(1979,2024):
    
    # Open dataset
    with xr.open_dataset(eulerian_file %year) as ds:
        
        # Load AR rank
        rank_year = ds['rank'].load()
        
        # Append the data of the year
        ar_ranks.append(rank_year[:,:,chunk[0]:chunk[1]])    

# Concatenate the years along the time axis
ar_ranks = np.concatenate(ar_ranks, axis=0)

# Load latitude and longitude
eulerian_dataset = xr.open_dataset(eulerian_file %1979)  
lat = eulerian_dataset['latitude'].values  # extract/copy the data
lon = eulerian_dataset['longitude'].values
eulerian_dataset.close()

# Load continents mask
cont_mask = rioxarray.open_rasterio(continents_folder + 'world_continents_3000sqkm.tif')

# Check grids and extract data
if np.array_equal(cont_mask['x'].values, lon) and np.array_equal(
        np.flip(cont_mask['y'].values, axis=0), lat):
    cont_mask = np.flip(np.squeeze(cont_mask.values), axis=0)
else:
   print('Continents mask is on an incorrect grid')
   
# Extract the chunk 
cont_mask = cont_mask[:,chunk[0]:chunk[1]]
lon = lon[chunk[0]:chunk[1]]
nlat, nlon = lat.size, lon.size

# %% COMPUTATION OF HISTOGRAMS 

## CONTINENTAL LOOP
conti_rank_count = np.zeros((8, 5, nlat, nlon), dtype=int)

for k in range(5):
    
    print('Continental loop --- Rank %s' %(k+1))
    
    # Constrain to continent
    for nconti in range(1,9):
        
        a_mask = np.where(cont_mask==nconti, 1, 0)
        conti_rank = ar_ranks*a_mask
        
        # Constrain to rank
        spec_rank = np.where(conti_rank == (k+1), conti_rank, 0)
        
        # Calculate number of AR events
        for i in range(nlat):
            for j in range(nlon):
                
                # Grab rank time series
                tmp_rank = spec_rank[:, i, j]
                tmp_bool, _ = rle(tmp_rank)
                conti_rank_count[nconti-1, k, i,j] = np.sum(tmp_bool>0)
               
np.save(output_folder + 'figure09/cont_rank_count_%s.npy' %task, conti_rank_count)

# RANK LOOP
rank_count = np.zeros((5, nlat, nlon), dtype=int)

for k in range(5):
    
    print('Rank loop --- Rank %s' %(k+1))
    
    # Mask specific rank
    rank_k = np.where(ar_ranks==(k+1), 1, 0)
    
    # Calculate number of AR events    
    for i in range(nlat):
        for j in range(nlon):
            
            # Extract rank time series
            rank_loc = rank_k[:, i, j]
            
            # Compute number of AR events
            tmp_bool, _ = rle(rank_loc)
            rank_count[k, i, j] = np.sum(tmp_bool>0)
                   
np.save(output_folder + 'figure10/rank_count_%s.npy' %task, rank_count)
