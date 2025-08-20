# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes the statistics necessary to plot Fig. 7 of 
Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python analysis_figure7.py
'''

# %% MODULES:
 
# Standard library imports
import numpy as np
import pandas as pd

# %% PARAMETERS

# Folders
data_folder = ...
validation_folder = ...
output_folder = ...

# Files
target_file = 'tARget_globalARcatalog_ERA5_1940-2023_v4_converted_including6h.pkl'
arconn_file = 'AR-CONNECT.MERRA2.1983-2016.Object_characteristics.csv'
ipart_file = '1983-2016_0,75_IPART_catalog.csv'
pikart_file = '1940-2023_PIKARTV1_Vallejo-Bernal_Braun_etal_2025.csv'

# %% tARget V4:

# Import catalog
d_target = pd.read_pickle(validation_folder + target_file)
d_target['time'] = pd.to_datetime(d_target['time'])
d_target = d_target[(d_target['time'] >= '1983-01-01') & (d_target['time'] < '2017-01-01 00:00:00')]

# Track IDs
a_allids = np.unique(d_target['trackid'])
Nids = a_allids.size

## LOOP
nlf_target = 0
a_ivt_target, a_lifetime_target, a_lat_target = np.zeros(Nids), np.zeros(Nids), np.zeros(Nids)
for i in range(Nids):
    ID = a_allids[i]
    AR = d_target[d_target['trackid'] == ID]
    
    # strength
    IVT = AR['IVT_mean'].mean()
    a_ivt_target[i] = IVT

    #duration
    dur = (pd.to_datetime(AR['time']).iloc[-1] - pd.to_datetime(AR['time']).iloc[0])/pd.Timedelta('1 hour')+6
    a_lifetime_target[i] = int(dur)
    
    #latitude
    tmp_centr = AR.centroid_y
    a_lat_target[i] = np.nanmean(tmp_centr)
    
    # landfall
    if np.any(~np.isnan(AR.lf_lon)):
        nlf_target += 1

# %% AR-CONNECT

## Load for ARCONNECT
d_arconn = pd.read_csv(validation_folder + arconn_file, encoding='latin1')
a_ivt_arconn = d_arconn['Average IVT [kg/(ms)]'].values
a_lifetime_arconn = d_arconn['Duration [hr]'].values
a_lat_arconn = d_arconn['Lifecycle Central Latitude'].values
## AR-CONNECT ain't got no landfall...

# %% IPART

# Import catalog
d_ipart = pd.read_csv(validation_folder + ipart_file)
d_ipart['time'] = pd.to_datetime(d_ipart['time'])

# Track IDs
a_allids = np.unique(d_ipart['trackid'])
Nids = a_allids.size

## LOOP
nlf_ipart = 0
a_ivt_ipart, a_lifetime_ipart, a_lat_ipart = np.zeros(Nids), np.zeros(Nids), np.zeros(Nids)
for i in range(Nids):
    ID = a_allids[i]
    AR = d_ipart[d_ipart['trackid'] == ID]
    
    # strength
    IVT = AR['IVT_mean'].mean()
    a_ivt_ipart[i] = IVT

    #duration
    dur = (pd.to_datetime(AR['time']).iloc[-1] - pd.to_datetime(AR['time']).iloc[0])/pd.Timedelta('1 hour')+6
    a_lifetime_ipart[i] = int(dur)
    
    #latitude
    tmp_centr = AR['centroid_y']
    a_lat_ipart[i] = np.nanmean(tmp_centr)

# no landfall for IPART...

# %% PIKART

# Open csv file with our AR catalog
d_pikart = pd.read_csv(data_folder + pikart_file)
d_pikart['time'] = pd.to_datetime(d_pikart['time'])
d_pikart = d_pikart[(d_pikart['time'] >= '1983-01-01 00:00:00') & (d_pikart['time'] < '2017-01-01 00:00:00')]

# Track IDs
a_allids = np.unique(d_pikart['trackid'])
Nids = a_allids.size

## LOOP
nlf_pikart = 0
a_ivt_pikart, a_lifetime_pikart, a_lat_pikart, a_level = np.zeros(Nids), np.zeros(Nids), np.zeros(Nids), np.zeros(Nids)
for i in range(Nids):
    ID = a_allids[i]
    AR = d_pikart[d_pikart['trackid'] == ID]
    
    # strength
    IVT = AR['mean_ivt'].mean()
    a_ivt_pikart[i] = IVT

    #duration
    dur = AR['lifetime'].iloc[0]
    a_lifetime_pikart[i] = dur
    
    #latitude
    tmp_centr = AR['centroid_lat']
    a_lat_pikart[i] = np.nanmean(tmp_centr)
    
    #landfall: comment-in with new catalog AFTER classification!
    if np.any(~np.isnan(AR.lf_lon)): 
        nlf_pikart += 1

# %% STATS for table 1

"""
You can write the table out as a csv and then simply copy & paste the new values.
Or just use the Spyder variable explorer...
"""

d_stats = pd.DataFrame({'time period': ['1979-2020', '1983-2016', '1940-2024', '1940-2023'],
                        'spatial extent': ['extratropics', 'extratropics', 'global', 'global'],
                        'temporal resolution': ['6h', '3h', '6h', '6h'],
                        'spatial resolution': ['0.75 x 0.75', '0.50 x 0.625', '0.25 x 0.25', '0.5 x 0.5'],
                        '# AR tracks': [a_lat_ipart.size, a_lat_arconn.size, a_lat_target[~np.isnan(a_lat_target)].size, Nids],
                        'median_IVT': [np.nanmedian(a_ivt_ipart), np.nanmedian(a_ivt_arconn), np.nanmedian(a_ivt_target), np.nanmedian(a_ivt_pikart)],
                        'q_IVT': [np.nanquantile(a_ivt_ipart,.95), np.nanquantile(a_ivt_arconn,.95), np.nanquantile(a_ivt_target,.95), np.nanquantile(a_ivt_pikart,.95)],
                        'median_dur': np.hstack([np.nanmedian(a_lifetime_ipart), np.nanmedian(a_lifetime_arconn), np.nanmedian(a_lifetime_target[a_lifetime_target>6]), np.nanmedian(a_lifetime_pikart)])/24,
                        'q_dur': np.hstack([np.nanquantile(a_lifetime_ipart,.95), np.nanquantile(a_lifetime_arconn,.95), np.nanquantile(a_lifetime_target,.95), np.nanquantile(a_lifetime_pikart,.95)])/24,
                        'lat_north': [np.mean(a_lat_ipart[a_lat_ipart>23.5]), np.mean(a_lat_arconn[a_lat_arconn>23.5]), np.mean(a_lat_target[a_lat_target>23.5]), np.mean(a_lat_pikart[a_lat_pikart>23.5])],
                        'lat_south': [np.mean(a_lat_ipart[a_lat_ipart<-23.5]), np.mean(a_lat_arconn[a_lat_arconn<-23.5]), np.mean(a_lat_target[a_lat_target<-23.5]), np.mean(a_lat_pikart[a_lat_pikart<-23.5])],
                        'nlf': [np.nan, np.nan, nlf_target, nlf_pikart],
                        'nlf_frac': [np.nan, np.nan, nlf_target/a_lat_target[~np.isnan(a_lat_target)].size, nlf_pikart/Nids]
                        },
                       index=(['IPART', 'AR-CONNECT', 'tARget', 'PIKART']))
d_stats.to_csv(output_folder + 'validation_stats.csv', index=True)

# %% SAVE DATA NECESSARY FOR THE FIGURE 

# tARget V4
np.save(output_folder + 'mean_ivt_target.npy', a_ivt_target)
np.save(output_folder + 'lifetime_target.npy', a_lifetime_target)
np.save(output_folder + 'lat_target.npy', a_lat_target)

# AR-CONNECT
np.save(output_folder + 'mean_ivt_arconn.npy', a_ivt_arconn)
np.save(output_folder + 'lifetime_arconn.npy', a_lifetime_arconn)
np.save(output_folder + 'lat_arconn.npy', a_lat_arconn)

# IPART
np.save(output_folder + 'mean_ivt_ipart.npy', a_ivt_ipart)
np.save(output_folder + 'lifetime_ipart.npy', a_lifetime_ipart)
np.save(output_folder + 'lat_ipart.npy', a_lat_ipart)

# PIKART
np.save(output_folder + 'mean_ivt_pikart.npy', a_ivt_pikart)
np.save(output_folder + 'lifetime_pikart.npy', a_lifetime_pikart)
np.save(output_folder + 'lat_pikart.npy', a_lat_pikart)
