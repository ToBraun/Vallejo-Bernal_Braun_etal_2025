# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script computes the inland penetration statistics necessary to plot Fig. 11 of 
Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python analysis_figure11.py
'''

# %% IMPORT MODULES

import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import scipy.stats as st

# %% PARAMETERS

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Years 
years = range(1940, 2024)

# Folders and files
output_base_name = config['output_base_name']
footprints_file = os.path.join(
    config['output_folder'], 'continuous_trackids', 
    output_base_name %('AR_footprints', '%s', 'nc'))
lagrangian_file = '1940-2023_PIKARTV1_Vallejo-Bernal_Braun_etal_2025.csv'

continents_folder = ...
output_folder = './output/figure11/'
os.makedirs(output_folder, exist_ok=True)

# %% LOAD DATA

# Open csv file with AR records
ars = pd.read_csv(lagrangian_file)
ars['time'] = pd.to_datetime(ars['time'])
lf_continents = ['na', 'sa', 'eu', 'af', 'as', 'au', 'oc', 'an']

# %% FIGURE 11a

# %% LIFE CYCLES

# Unique AR trajectories
tracks = ars.drop_duplicates(['trackid'])
nid = np.size(tracks,0)

# How many ARs are IP?
nip = np.sum(tracks['inland_pen'])

# How many ARs are LF but NOT IP?
nlf = np.sum(tracks['land_falling'] & ~tracks['inland_pen'])

# PARAMETERS
## INTERPOLATION GRID
nbins = 12
bins = np.linspace(0,1,nbins)
# SIZE OF POPULATIONS
print('Number of landfalling ARs that do NOT inland penetrate: ', nlf)
print('Number of inland penetrating ARs: ', nip)

# LOOP
ivtlc, ivtlc_lf, ivtlc_ip = np.nan*np.ones((nid, nbins)), np.nan*np.ones((nlf, nbins)), np.nan*np.ones((nip, nbins))
lftime, l_iptime = np.nan*np.ones(nlf), []
allcounts, lfcounts, IPcounts = np.zeros(nbins), np.zeros(nbins), np.zeros(nbins)
klf, kip = 0, 0
for i in range(nid):
    
    # Extract AR
    AR = ars[ars['trackid'] == tracks['trackid'].iloc[i]]
    L = AR.shape[0]
    
    # There's one AR with only 6h for some reason (probably due to year clipping), kick it out...
    if L == 1:
        print(tracks['trackid'].iloc[i])
        continue
    
    # Lifecycle axis with birth-death approach
    tmp_lc = np.array(AR['life_cycle'])  
    # IVT
    tmp_ivt = np.array(AR['mean_ivt'])

    # Assign to bins
    tmp_assign = np.digitize(tmp_lc, bins)
    tmp_bins = np.unique(tmp_assign)
    
    # Compute binwise means
    tmp_ivtlc = np.zeros(nbins)
    tmp_count = np.zeros(nbins)
    for nbin in tmp_bins:
        tmp_ivtlc[nbin-1] += np.mean(tmp_ivt[tmp_assign == nbin]) 
        tmp_count[nbin-1] += 1

    # Interpolate IVT values to common time steps
    interpolation_func = interp1d(tmp_lc, tmp_ivt, kind='linear')
    tmp_ivtlc = interpolation_func(bins)
    
    ivtlc[i,] = tmp_ivtlc
    allcounts += np.ones(nbins)
    
    # Only inland penetrating ARs
    if np.any(AR['inland_pen']):
        ivtlc_ip[kip,] = tmp_ivtlc
        IPcounts += np.where(tmp_ivtlc>0, 1, 0)

        # IP times
        ## extract longest continuous interval of IP (i.e., more than 30% of cont. intersection)
        tmp_contiidx = np.where(AR['land'] > 30)[0]
        tmp_seq = max(np.split(tmp_contiidx, np.where(np.diff(tmp_contiidx) != 1)[0] + 1), key=len).tolist()
        ## locate IP point at the center of that interval
        ipIDX = tmp_lc[tmp_seq[int(len(tmp_seq)/2)]]
        l_iptime.append(ipIDX)
        kip += 1
        
    # Only landfalling ARs
    elif np.any(~np.isnan(AR['lf_lat'])) and not np.any(AR['inland_pen']):
        ivtlc_lf[klf,] = tmp_ivtlc
        lfcounts += np.where(tmp_ivtlc>0, 1, 0)

        # LF times
        idx = np.where(~np.isnan(AR['lf_lat']))[0][0]
        lftime[klf] = tmp_lc[idx]
        klf += 1

# %% RESULTS

## TIMES
## exclude landfalls and inland penetrations at lifecycle instances 0 and 1 for better visibility
lifetimes_disp = lftime[(lftime > 0) & (lftime < 1)]
iptimes = np.hstack(l_iptime)
iptimes_disp = iptimes[(iptimes > 0) & (iptimes <1)]

# MEAN IVT CURVES
med_ivt_lc0 = np.nansum(ivtlc,0)/allcounts
med_ivt_lf0 = np.nansum(ivtlc_lf,0)/lfcounts
med_ivt_ip0 = np.nansum(ivtlc_ip,0)/IPcounts

# normalise by first element
med_ivt_lc = med_ivt_lc0/med_ivt_lc0[0]
med_ivt_lf = med_ivt_lf0/med_ivt_lf0[0]
med_ivt_ip = med_ivt_ip0/med_ivt_ip0[0]

# RATES
rate_all = np.diff(med_ivt_lc)/med_ivt_lc[:-1]
rate_lf = np.diff(med_ivt_lf)/med_ivt_lf[:-1]
rate_ip = np.diff(med_ivt_ip)/med_ivt_ip[:-1]

# %% SAVE OUTPUT

np.save(output_folder + 'bins.npy', bins)

np.save(output_folder + 'lifetimes_disp.npy', lifetimes_disp)
np.save(output_folder + 'iptimes_disp.npy', iptimes_disp)

np.save(output_folder + 'med_ivt_lc.npy', med_ivt_lc)
np.save(output_folder + 'med_ivt_lf.npy', med_ivt_lf)
np.save(output_folder + 'med_ivt_ip.npy', med_ivt_ip)

np.save(output_folder + 'rate_all.npy', rate_all) 
np.save(output_folder + 'rate_lf.npy', rate_lf)
np.save(output_folder + 'rate_ip.npy', rate_ip)

print('Results of the life cycle analysis saved')

# %% FIGURE 11b

# %%

def kernel_map(X, Y, q, k=100):
    """
    Generates a kernel density estimation map.
    
    Parameters
    ----------
    X : array-like
        Array of x-coordinates for data points.
    Y : array-like
        Array of y-coordinates for data points.
    q : float
        Quantile value for pruning the kernel density estimation.
    k : int, optional
        Number of grid points along each axis, default is 100.
    
    Returns
    ------- 
    tx : ndarray
        Meshgrid of x-coordinates for the generated map.
    ty : ndarray
        Meshgrid of y-coordinates for the generated map.
    a_ckernel : ndarray
        Pruned kernel density estimation values corresponding to the grid points.
    """
    #src: https://ipython-books.github.io/76-estimating-a-probability-distribution-nonparametrically-with-a-kernel-density-estimation/
    # projection
    crs = ccrs.PlateCarree(central_longitude=0)
    geo = ccrs.Geodetic()
    # Coordinates of the four corners of the map.
    ax = plt.axes(projection=crs)
    x0, x1, y0, y1 = ax.get_extent()
    #k = 100
    
    ### GENESIS CENTROIDS
    h = crs.transform_points(geo, X, Y)[:, :2].T
    hnonan = h[:,np.where(np.all(~np.apply_along_axis(np.isnan, 0, h),0))[0]]
    # kernel
    kde = st.gaussian_kde(hnonan)
    # We create the grid.
    tx, ty = np.meshgrid(np.linspace(x0, x1, 2 * k),
                         np.linspace(y0, y1, k))
    # We reshape the grid for the kde() function.
    mesh = np.vstack((tx.ravel(), ty.ravel()))
    # We evaluate the kde() function on the grid.
    v = kde(mesh).reshape((k, 2 * k))
    # prune
    a_ckernel = np.where(v>np.quantile(v, q), v, np.nan)
    
    return tx, ty, a_ckernel

# %% GENESIS AND TERMINATION PATCHES

# COLLECT THE INTERSECTION COLUMNS
li_cols = []
for cont in lf_continents:
    for comp in ['conti', 'insu']:
        li_cols.append(['%s_%s_lon' %(comp, cont), '%s_%s_lat' %(comp, cont)])

### GENESIS AND TERMINATION
l_genesis_locs, l_term_locs = [], []
for i in range(nid):
    
    # Extract AR
    AR = ars[ars['trackid'] == tracks['trackid'].iloc[i]]

    # Genesis centroids
    l_genesis_locs.append((AR['centroid_lon'].iloc[0], AR['centroid_lat'].iloc[0]))

    # Termination land intersection
    ARterm = AR.tail(1)
    tmp_LIs = np.vstack([ARterm[li_cols[j]] for j in range(np.size(li_cols,0))])
    l_term_locs.append(tmp_LIs[np.where(np.sum(np.isnan(tmp_LIs), 1)==0)[0],:])

# KERNEL DENSITY ESTIMATES
txg, tyg, g_kernel = kernel_map(np.vstack(l_genesis_locs)[:,0], np.vstack(l_genesis_locs)[:,1], 0.8)
txt, tyt, t_kernel = kernel_map(np.vstack(l_term_locs)[:,0], np.vstack(l_term_locs)[:,1], 0.8)

# %% SAVE OUTPUT

np.save(output_folder + 'txg.npy', txg)
np.save(output_folder + 'tyg.npy', tyg)
np.save(output_folder + 'g_kernel.npy', g_kernel)

np.save(output_folder + 'txt.npy', txt)
np.save(output_folder + 'tyt.npy', tyt)
np.save(output_folder + 't_kernel.npy', t_kernel)

print('Results of the genesis-termination analysis saved')

# %% INLAND PENETRATION 

# %% LOAD DATA
   
# Load latitude and longitude
footprints_dataset = xr.open_dataset(footprints_file %years[0])  
lat = footprints_dataset['latitude'].values  # extract/copy the data
lon = footprints_dataset['longitude'].values
footprints_dataset.close()

# Load continents mask
cont_mask = rioxarray.open_rasterio(continents_folder + 'world_continents_3000sqkm.tif')

# Check grids and extract data
if np.array_equal(cont_mask['x'].values, lon) and np.array_equal(np.flip(cont_mask['y'].values, axis=0), lat):
    cont_mask = np.flip(np.squeeze(cont_mask.values), axis=0)
else:
   print('Continents mask is on an incorrect grid')
   
# Mask water and land
land_mask = np.where(cont_mask > 0, 1, 0)

# Filter inland penetrating ARs
ip_ars = ars[(ars['inland_pen'])]
   
# %% FREQUENCY OF INLAND-PENETRATING ARs

# Initialize variables
ip_count = []
total_time = 0

# Loop through years
for year in years:
    
    print(year)
    
    # Open dataset
    footprints_dataset = xr.open_dataset(footprints_file %year)
    
    # Load time
    time = pd.to_datetime(footprints_dataset['time'].values)
    
    # Add the total number of time steps in the specific year
    total_time += len(time)

    # Initialize variable
    ip_count_year = np.zeros((len(time), len(lat), len(lon)))
    
    # Loop through time
    for t in range(len(time)):
        
        # Extract active ARs
        ar_record = ip_ars[ip_ars['time'] == time[t]]
        
        # If there are inland-penetrating ARs active
        if np.size(ar_record,0) > 0:
            
            # Extract footprints
            footprints_t = footprints_dataset.sel(time=time[t])['footprint'].values
            
            # Get footprints of inland-penetrating ARs
            ip_trackids = np.unique(ar_record['trackid'])
            
            # Put a 1 in the grid cells that are over land and under the footprint of an inland-penetrating AR
            ip_count_year[t,:,:] = land_mask * np.isin(footprints_t, ip_trackids) 
            
    # Calculate the total number of time steps with inland-penetrating ARs
    ip_count.append(np.sum(ip_count_year, axis=0))
    
# Concatenate the years along the time axis
ip_count = np.stack(ip_count, axis=0)

# Calculate inland-penetration frequency
ip_freq = 100*np.sum(ip_count, axis=0)/total_time
print('Total number of time steps: %s' %total_time)

# %% SAVE OUTPUT

np.save(output_folder + 'ip_freq.npy', ip_freq)
np.save(output_folder + 'lon_pikart.npy', lon)
np.save(output_folder + 'lat_pikart.npy', lat)

print('Results of the inland penetration analysis saved')
