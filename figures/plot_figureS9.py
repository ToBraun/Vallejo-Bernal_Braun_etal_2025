# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# Tobias Braun <tobraun@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. S9 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figureS9.py
'''

# %% IMPORT MODULES

import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cmcrameri import cm
from plot_functions import plot_landmasses, maxfactor_map
from cartopy.util import add_cyclic_point

import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

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
year0 = 1940
yearf = 2023
years = yearf - year0 + 1
central_lon = -150
color_map = cm.bilbao_r

# Folders and files
input_folder = root_folder + 'data/figureS9/'
figure_file = root_folder + 'manuscript/PIKART_FigureS9.png'

# %% LOAD DATA

# Coordinates
lat = np.arange(-89.75, 89.75 + 0.5, 0.5)
lon = np.arange(-180, 180, 0.5)

# AR conditions frequency
ar_freq = np.load(input_folder + '%s-%s_ar_freq.npy' %(year0, yearf))

# AR landfall
unique_lons = np.load(input_folder + '%s-%s_unique_lons.npy' %(year0, yearf))
unique_lats = np.load(input_folder + '%s-%s_unique_lats.npy' %(year0, yearf))
lf_map = np.load(input_folder + '%s-%s_lf_map.npy' %(year0, yearf))/years

# Add cyclic point
ar_freq, lon_freq = add_cyclic_point(ar_freq, coord=lon, axis=1)

# Calculate number of days per year
ar_freq = 365*ar_freq/100

# %% PLOT

# Create the figure
figure = plt.figure(figsize=(180/25.4,93/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# Create the subplot axis
ax = figure.add_subplot(1,1,1, projection=proj)
ax.set_global()

# Plot the landmasses
ax = plot_landmasses(ax)

# PLOT CONTOURS
map_freq = ax.pcolormesh(lon_freq, lat, ar_freq, shading='nearest',
                         cmap=color_map, vmin=0, vmax=70, 
                         transform=ccrs.PlateCarree(), zorder=0)

#PLOT BUBBLES
pltlons, pltlats = np.meshgrid(unique_lons, unique_lats)
nonzero_mask = lf_map > 0
lf_scatter = ax.scatter(pltlons[nonzero_mask], pltlats[nonzero_mask], 
                        transform=ccrs.PlateCarree(), zorder=5,
                        s=10*lf_map[nonzero_mask], alpha=.7, color='steelblue')

# Maxfactor the map
ax = maxfactor_map(ax, draw_grid=False)
ax.set_aspect('equal')

# Color bar
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size='1.8%', pad=0.1, axes_class=plt.Axes)
figure.add_axes(ax_cb)
cbar = plt.colorbar(map_freq, cax=ax_cb, extend='max')
cbar.ax.tick_params(length=2.5, width=0.6)
cbar.set_label(label='AR conditions (days/year)', labelpad=8)

# Add landfall legend
handles, labels = lf_scatter.legend_elements(
    prop='sizes', num=[10, 50, 100, 150], alpha=0.7, color='steelblue')
legend = ax.legend(handles, ['1', '5', '10', '15'], ncol=4, loc='upper left', 
                   bbox_to_anchor=(0.6, 1.1), edgecolor='none', framealpha=1)

ax.annotate('Land-falling ARs/year:', xy=(35, 96), xycoords='data', 
            ha='right', va='bottom')

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
