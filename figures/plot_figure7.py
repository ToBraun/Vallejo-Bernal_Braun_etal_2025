# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. 7 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure7.py
'''

# %% IMPORT MODULES

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

# Update Matplotlib parameters
plt.rcParams.update({'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'font.family': 'Arial'})

mpl.rcParams.update({'axes.linewidth': 0.75,
                     'axes.edgecolor': '#262626',
                     'xtick.color':'#262626',
                     'ytick.color':'#262626',
                     'font.size': 8,
                     'axes.titlesize': 8.5})

root_folder = ...

# %% INITIAL SETTINGS

# Folders and files
analysis_folder = root_folder + '/data/figure07/'
figure_file = root_folder + 'manuscript/PIKART_Figure7.png'

# %% LOAD DATA NECESSARY FOR THE FIGURE 

# tARget V4
ivt_target = np.load(analysis_folder + 'mean_ivt_target.npy')
lifetime_target = np.load(analysis_folder + 'lifetime_target.npy')
lat_target = np.load(analysis_folder + 'lat_target.npy')

# AR-CONNECT
ivt_arconn = np.load(analysis_folder + 'mean_ivt_arconn.npy')
lifetime_arconn = np.load(analysis_folder + 'lifetime_arconn.npy')
lat_arconn = np.load(analysis_folder + 'lat_arconn.npy')

# IPART
ivt_ipart = np.load(analysis_folder + 'mean_ivt_ipart.npy')
lifetime_ipart = np.load(analysis_folder + 'lifetime_ipart.npy')
lat_ipart = np.load(analysis_folder + 'lat_ipart.npy')

# PIKART
ivt_pikart = np.load(analysis_folder + 'mean_ivt_pikart.npy')
lifetime_pikart = np.load(analysis_folder + 'lifetime_pikart.npy')
lat_pikart = np.load(analysis_folder + 'lat_pikart.npy')

# %% HISTOGRAMS

## Percentage of 6h-tracks in tARget:
print(np.nansum(lifetime_target == 6)/lifetime_target[~np.isnan(lifetime_target)].size)
# how many tracks are there with dur>6h?
print(lifetime_target[~np.isnan(lifetime_target)].size - np.nansum(lifetime_target == 6))

# Lifetime
lifetimehist_target, tbins = np.histogram(lifetime_target, bins=np.arange(0, 400, 6), density=False)
lifetimehist_arconn, abins = np.histogram(lifetime_arconn, bins=np.arange(0, 400, 6), density=False)
lifetimehist_ipart, ipbins = np.histogram(lifetime_ipart, bins=np.arange(0, 400, 6), density=False)
lifetimehist_pikart, pbins = np.histogram(lifetime_pikart, bins=np.arange(0, 400, 6), density=False)

# Convert to float. Necessary to mask values
lifetimehist_target = lifetimehist_target.astype('float32')
lifetimehist_arconn = lifetimehist_arconn.astype('float32')
lifetimehist_ipart = lifetimehist_ipart.astype('float32')
lifetimehist_pikart = lifetimehist_pikart.astype('float32')

# Mask values
lifetimehist_pikart[:2] = np.nan # due to cropping
lifetimehist_ipart[:2] = np.nan
lifetimehist_target[0] = np.nan
lifetimehist_arconn[3] = np.nan

# Mean IVT
kbins = 50
ivthist_target, ibins = np.histogram(ivt_target, bins=np.linspace(0, 1000, kbins), density=False)
ivthist_target_long, ibins = np.histogram(ivt_target[lifetime_target>6], bins=np.linspace(0, 1000, kbins), density=False)
ivthist_arconn, _ = np.histogram(ivt_arconn, bins=np.linspace(0, 1000, kbins), density=False)
ivthist_ipart, _ = np.histogram(ivt_ipart, bins=np.linspace(0, 1000, kbins), density=False)
ivthist_pikart, _ = np.histogram(ivt_pikart[ivt_pikart>0], bins=np.linspace(0, 1000, kbins), density=False)
ibins = ibins[1:] - np.diff(ibins)[0]/2

# Centroid latitude
kbins = 75
lathist_target, lbins = np.histogram(lat_target, bins=np.linspace(-90, 90, kbins), density=False)
lathist_target_long, lbins = np.histogram(lat_target[lifetime_target>6], bins=np.linspace(-90, 90, kbins), density=False)
lathist_arconn, _ = np.histogram(lat_arconn, bins=np.linspace(-90, 90, kbins), density=False)
lathist_ipart, _ = np.histogram(lat_ipart, bins=np.linspace(-90, 90, kbins), density=False)
lathist_pikart, _ = np.histogram(lat_pikart, bins=np.linspace(-90, 90, kbins), density=False)
lbins = lbins[1:] - np.diff(lbins)[0]/2

# Mask out tropics for catalogs that don't have them
lathist_arconn = np.where(lathist_arconn == 0, np.nan, lathist_arconn)
lathist_ipart = np.where(lathist_ipart == 0, np.nan, lathist_ipart)

# %% PLOT

# Create the figure
figure = plt.figure(figsize=(120/25.4,200/25.4), dpi=500)

# Create the subplot axis
ax1 = figure.add_subplot(3,1,1)

# Lifecycle panel
ax1.plot(tbins.astype(int)[:-1], lifetimehist_target, color='goldenrod', linewidth=0.75, label='tARget-4'); 
ax1.plot(abins.astype(int)[:-1], lifetimehist_arconn, color='mediumseagreen', linewidth=0.75, label='AR-CONNECT'); 
ax1.plot(ipbins.astype(int)[:-1], lifetimehist_ipart, color='royalblue', linewidth=0.75, label='IPART'); 
ax1.plot(pbins.astype(int)[:-1], lifetimehist_pikart, color='black', linewidth=1.5, label='PIKART'); 
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(top=70000)
ax1.set_xticks([6,12,24,48,96,192,384])
ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax1.xaxis.set_minor_locator(mpl.ticker.NullLocator())
ax1.set_xlabel('Lifetime (hours)')
ax1.set_ylabel('Number of AR trajectories')
ax1.legend(loc='lower left', edgecolor='none', framealpha=1)

# Create the subplot axis
ax2 = figure.add_subplot(3,1,2)

ax2.plot(ibins.astype(int), ivthist_target, color='goldenrod', linewidth=0.75, label='tARget-4'); 
ax2.plot(ibins.astype(int), ivthist_target_long, linestyle="dotted", color='goldenrod', linewidth=0.75, label='tARget V3'); 
ax2.plot(ibins.astype(int), ivthist_arconn, color='mediumseagreen', linewidth=0.75, label='AR-CONNECT'); 
ax2.plot(ibins.astype(int), ivthist_ipart, color='royalblue', linewidth=0.75, label='IPART'); 
ax2.plot(ibins.astype(int), ivthist_pikart, color='black', linewidth=1.5, label='PIKART'); 
ax2.set_ylim(top=7000)
ax2.set_xlabel('Mean IVT (kg m$^{-1}$ s$^{-1}$)')
ax2.set_ylabel('Number of AR trajectories')

# Create the subplot axis
ax3 = figure.add_subplot(3,1,3)

ax3.plot(lbins.astype(int), lathist_target, color='goldenrod', linewidth=0.75, label='tARget-4'); 
ax3.plot(lbins.astype(int), lathist_target_long, linestyle="dotted", color='goldenrod', linewidth=0.75, label='tARget V3'); 
ax3.plot(lbins.astype(int), lathist_arconn, color='mediumseagreen', linewidth=0.75, label='AR-CONNECT'); 
ax3.plot(lbins.astype(int), lathist_ipart, color='royalblue', linewidth=0.75, label='IPART'); 
ax3.plot(lbins.astype(int), lathist_pikart, color='black', linewidth=1.5, label='PIKART'); 
ax3.set_xlim(-90, 90); ax3.set_xticks(np.arange(-90,120,30)); ax3.set_xticklabels(['90°S', '60°S', '30°S', '0', '30°N', '60°N', '90°N'])
ax3.set_ylim(top=3500)
ax3.set_xlabel('Average centroid latitude')
ax3.set_ylabel('Number of AR trajectories')

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
