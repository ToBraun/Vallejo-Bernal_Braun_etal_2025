# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes the statistics of Table S2 in 
Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python analysis_figureS8-1.py
'''

# %% MODULES
 
# Standard library imports
import numpy as np
import pandas as pd

# %% PARAMETERS

# Folders
data_folder = ...
output_folder = ...

# Files
pikart_file_p85 = '1940-2023_PIKARTV1_p85_Vallejo-Bernal_Braun_etal_2025.csv'
pikart_file_p875 = '1940-2023_PIKARTV1_Vallejo-Bernal_Braun_etal_2025.csv'
pikart_file_p90 = '1940-2023_PIKARTV1_p90_Vallejo-Bernal_Braun_etal_2025.csv'

# %% PIKART PERCENTILE 85TH

# Initialize variables
id_count = np.zeros(3)
median_ivt = np.zeros(3)
q_ivt = np.zeros(3)
median_dur = np.zeros(3)
q_dur = np.zeros(3)
lat_north = np.zeros(3)
lat_south = np.zeros(3)
lf_count = np.zeros(3)
lf_frac = np.zeros(3)
files = [pikart_file_p85, pikart_file_p875, pikart_file_p90]

for n in range(3):
    
    # Open csv file with our AR catalog
    pikart = pd.read_csv(data_folder + files[n])
    pikart['time'] = pd.to_datetime(pikart['time'])
    pikart = pikart[(pikart['time'] >= '1983-01-01 00:00:00') & 
                    (pikart['time'] < '2017-01-01 00:00:00')]
    
    # Track IDs
    all_ids = np.unique(pikart['trackid'])
    nids = all_ids.size
    
    ## LOOP
    nlf = 0
    ivt, lifetime, lat, level = np.zeros(nids), np.zeros(nids), np.zeros(nids), np.zeros(nids)
    for i in range(nids):
        trackid = all_ids[i]
        ar = pikart[pikart['trackid'] == trackid]
        
        # strength
        ar_ivt = ar['mean_ivt'].mean()
        ivt[i] = ar_ivt
    
        #duration
        dur = ar['lifetime'].iloc[0]
        lifetime[i] = dur
        
        #latitude
        tmp_centr = ar['centroid_lat']
        lat[i] = np.nanmean(tmp_centr)
        
        #landfall: comment-in with new catalog AFTER classification!
        if np.any(~np.isnan(ar.lf_lon)): 
            nlf += 1
            
    # Store results            
    id_count[n] = nids
    median_ivt[n] = np.nanmedian(ivt)
    q_ivt[n] = np.nanquantile(ivt,.95)
    median_dur[n] = np.nanmedian(lifetime/24)
    q_dur[n] = np.nanquantile(lifetime/24,.95)
    lat_north[n] = np.mean(lat[lat>23.5])
    lat_south[n] = np.mean(lat[lat<-23.5])
    lf_count[n] = nlf
    lf_frac[n] = nlf/nids

# %% STATS for table 1

d_stats = pd.DataFrame({'time period': ['1940-2023', '1940-2023', '1940-2023'],
                        'spatial extent': ['global', 'global', 'global'],
                        'temporal resolution': ['6h', '6h', '6h'],
                        'spatial resolution': ['0.5 x 0.5', '0.5 x 0.5', '0.5 x 0.5'],
                        '# AR tracks': id_count,
                        'median_IVT': median_ivt,
                        'q_IVT': q_ivt,
                        'median_dur': median_dur,
                        'q_dur': q_dur,
                        'lat_north': lat_north,
                        'lat_south': lat_south,
                        'nlf': lf_count,
                        'nlf_frac': lf_frac
                        },
                       index=['PIKART P85', 'PIKART P87.5', 'PIKART P90'])
d_stats.to_csv(output_folder + 'sensitivity_analsys.csv', index=True)
