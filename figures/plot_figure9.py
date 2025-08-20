# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. 9 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure9.py
'''

# %% IMPORT MODULES

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import rioxarray
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cmcrameri import cm
from matplotlib.patches import Patch

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

# Parameters
years = 2023 - 1979 + 1

# Files
continents_file = root_folder + 'data/world_continents/world_continents_3000sqkm.tif'
cont_rank_file = root_folder + 'data/figure09/cont_rank_count.npy'
figure_file = root_folder + 'manuscript/PIKART_Figure9.png'

# %% LOAD DATA

# Load continents mask
cont_mask = rioxarray.open_rasterio(continents_file)

# Calculate the number of pixels per continent
continents = np.hstack(['North America', 'South America', 'Europe', 'Africa', 
                        'Asia', 'Australia', 'Oceania', 'Antarctica'])
cont_pixels = np.hstack(
    [np.sum(cont_mask == n) for n in np.arange(1, continents.size+1)])

# Load histogram
cont_rank_count = np.load(cont_rank_file)

# %% COMPUTE HITOGRAMS

# Mean number of time steps with AR conditions for each level and each continent
cont_rank_count = np.transpose(np.sum(cont_rank_count, axis=(2,3)).T/cont_pixels)

# INSET: Mean frequency of ARs for each continent (no level distinction)
cont_freq = np.sum(cont_rank_count, 1)/years

# RIGHT AXIS: Mean frequency of ARs for each level (no continental distinction)
rank_freq = np.sum(cont_rank_count, 0)/years

# LEFT AXIS: Mean frequency of each rank for each continent
cont_rank_freq = cont_rank_count/years

# Order by total continental count
permut = np.flip(np.argsort(cont_freq))
cont_labels = continents[permut]
cont_rank_freq = cont_rank_freq[permut,]
cont_freq = cont_freq[permut]

# %% PLOT FIGURE 9

# LOCAL plotting parameters
colors = cm.bukavu(np.linspace(0,1, 8+1))
cont_abbrevs = np.array(['NA', 'SA', 'EU', 'AF', 'AS', 'AU', 'OC', 'AT'])
ranks = np.arange(1,6)
n_cont = 8
x_axis = np.arange(5)
x_shifts = np.linspace(-.35, .35, n_cont)

# Create the figure
figure = plt.figure(figsize=(180/25.4,110/25.4), dpi=500)

# Create the subplot axis
ax = figure.add_subplot(1,1,1)

# Plot the bars
ax.bar(x_axis, rank_freq**(1/2), fill=False, zorder=6, edgecolor='black', alpha=.8)
ax.bar(x_axis, rank_freq**(1/2), color='white', alpha=.25, zorder=1)

# Maxfactor the plot
ax.set_xticks(np.arange(5))
ax.set_xticklabels(np.arange(1,6))
ax.set_ylim([0,5])
ax.set_xlabel('AR rank')
ax.set_ylabel('Rank-specific frequency (AR events/year)');
ax.set_xlim(-0.6, 4.6)
ax.set_yticklabels(np.arange(0,6,1)**2)

# Legend elements
legend_elements = []
for i in range(np.size(continents)):
    legend_elements.append(Patch(
        edgecolor=colors[i,], facecolor=colors[i,]))
lgd = ax.legend(legend_elements, cont_labels,
            bbox_to_anchor=(0.5,-0.26), loc="lower center",
            ncol=4, edgecolor='w', framealpha=1) 

## INSET
axins = inset_axes(ax, width="40%", height="40%", loc=1, borderpad=2)

for i in range(8):
    axins.bar([i], cont_freq[i], color=colors[i,],)
    
axins.set_xticks(np.arange(8))
axins.set_xticklabels(cont_abbrevs[permut])
axins.set_ylim([0,11])
axins.set_yticks(np.arange(0,11,2))
axins.set_ylabel('Continent-specific\nfrequency (AR events/year)')

# SECOND AXIS
secx = ax.twinx()
for r in range(5):
    
    rank_cont = cont_rank_freq[:,r]
    
    for i in range(n_cont):    
    
        secx.bar(x_axis[r] + x_shifts[i], rank_cont[i]**(1/2), color=colors[i,], 
                 label = cont_labels[i], width=.5/8)

secx.set_ylim([0,3])
secx.set_ylabel('Continent-specific frequency (AR events/year)')
secx.set_yticks(np.arange(0,4,1))
secx.set_yticklabels(np.arange(0,4,1)**2)

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
