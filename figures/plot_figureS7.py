# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. S7 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figureS7.py
'''

# %% IMPORT MODULES

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
from cmcrameri import cm
from plot_functions import plot_landmasses, maxfactor_map_polar
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

# %% PARAMETERS

# Threshold for the Antarctic panel
threshold_sp = 2  

# Threshold for the Artic panel
threshold_np = 5  

# Color map
color_map = cm.batlowW_r
coast_color = '#505050'

# Files and folders
input_folder = root_folder + 'data/figureS7/'
figure_file = root_folder + 'manuscript/PIKART_FigureS7.png'

# %% LOAD AND PREPARE DATA

# Load PIKART
ar_freq_pikart = np.load(input_folder + '1980-2020_ar_freq_pikart.npy')
lon_pikart = np.load(input_folder + 'lon_pikart.npy')
lat_pikart = np.load(input_folder + 'lat_pikart.npy')

# Add cyclic point
ar_freq_pikart, lon_pikart = add_cyclic_point(
    ar_freq_pikart, coord=lon_pikart, axis=1)

# Load tARget
ar_freq_target = np.load(input_folder + '1980-2020_ar_freq_target.npy')
lon_target = np.load(input_folder + 'lon_target.npy')
lat_target = np.load(input_folder + 'lat_target.npy')

# Add cyclic point
ar_freq_target, lon_target = add_cyclic_point(
    ar_freq_target, coord=lon_target, axis=1)

# %% FIGURE LAYOUT

# Create the figure
figure = plt.figure(figsize=(180/25.4,180/25.4), dpi=500)

# Define the CartoPy projections
proj_np = ccrs.NorthPolarStereo()
proj_sp = ccrs.SouthPolarStereo()

# Define a circular boundary to crop the map within a circular region
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = Path(verts * radius + center)

# Create the layout
gs = gridspec.GridSpec(2, 3, width_ratios=[1,1,0.05], wspace=0.4)

# %% ANTARCTICA

##### PIKART

# Create the subplot axis
ax = figure.add_subplot(gs[0], projection=proj_sp)

# Set the extension of the map
ax.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())

# Plot the landmasses
ax = plot_landmasses(ax, ec=coast_color, hatch=False)

# Apply the circular boundary to the polar map
ax.set_boundary(circle, transform=ax.transAxes)

# Plot the AR frequency
cs = ax.pcolormesh(lon_pikart, lat_pikart, ar_freq_pikart, shading='nearest',
                   cmap=color_map, vmax=threshold_sp, zorder=0,
                   transform=ccrs.PlateCarree())

# Maxfactor the map
ax = maxfactor_map_polar(ax)
ax.set_title('PIKART', fontweight='bold')
ax.text(-0.10, 1, '(a)', ha='left', transform=ax.transAxes)

##### tARget V4

# Create the subplot axis
ax = figure.add_subplot(gs[1], projection=proj_sp)

# Set the extension of the map
ax.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())

# Plot the landmasses
ax = plot_landmasses(ax, ec=coast_color, hatch=False)

# Apply the circular boundary to the polar map
ax.set_boundary(circle, transform=ax.transAxes)

# Plot the AR frequency
cs = ax.pcolormesh(lon_target, lat_target, ar_freq_target, shading='nearest',
                   cmap=color_map, vmax=threshold_sp, zorder=0,
                   transform=ccrs.PlateCarree())

# Maxfactor the map
ax = maxfactor_map_polar(ax)
ax.set_title('tARget-4', fontweight='bold')
ax.text(-0.10, 1, '(b)', ha='left', transform=ax.transAxes)

# Add subplot to place the colorbar
ax = figure.add_subplot(gs[2])
ax.set_axis_off()
cbaxes = inset_axes(ax, width=0.12, height=2.1, loc='center left',
                    bbox_to_anchor=(-0.5, 0.5),
                    bbox_transform=ax.transData,
                    borderpad=0)

cbar = plt.colorbar(
    cs, cax=cbaxes, orientation='vertical', spacing='proportional',
    extend='max')
cbar.set_label('AR frequency (%)', labelpad=6)
cbar.ax.tick_params(length=2.5, width=0.6)
cbar.ax.set_yticks(np.arange(0,2.5,0.5))

# %% ARCTIC

##### PIKART

# Create the subplot axis
ax = figure.add_subplot(gs[3], projection=proj_np)

# Set the extension of the map
ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())

# Plot the landmasses
ax = plot_landmasses(ax, ec=coast_color, hatch=False)

# Apply the circular boundary to the polar map
ax.set_boundary(circle, transform=ax.transAxes)

# Plot the AR frequency
cs = ax.pcolormesh(lon_pikart, lat_pikart, ar_freq_pikart, shading='nearest',
                   cmap=color_map, vmax=threshold_np, zorder=0, 
                   transform=ccrs.PlateCarree())

# Maxfactor the map
ax = maxfactor_map_polar(ax)
ax.text(-0.10, 1, '(c)', ha='left', transform=ax.transAxes)

##### tARget V4

# Create the subplot axis
ax = figure.add_subplot(gs[4], projection=proj_np)

# Set the extension of the map
ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())

# Plot the landmasses
ax = plot_landmasses(ax, ec=coast_color, hatch=False)

# Apply the circular boundary to the polar map
ax.set_boundary(circle, transform=ax.transAxes)

# Plot the AR frequency
cs = ax.pcolormesh(lon_target, lat_target, ar_freq_target, shading='nearest',
                   cmap=color_map, vmax=threshold_np, zorder=0, 
                   transform=ccrs.PlateCarree())

# Maxfactor the map
ax = maxfactor_map_polar(ax)
ax.text(-0.10, 1, '(d)', ha='left', transform=ax.transAxes)

# Add subplot to place the colorbar
ax = figure.add_subplot(gs[5])
ax.set_axis_off()
cbaxes = inset_axes(ax, width=0.12, height=2.1, loc='center left',
                    bbox_to_anchor=(-0.5, 0.5),
                    bbox_transform=ax.transData,
                    borderpad=0)

# Add colorbar
cbar = plt.colorbar(
    cs, cax=cbaxes, orientation='vertical', spacing='proportional',
    extend='max')
cbar.set_label('AR frequency (%)', labelpad=6)
cbar.ax.tick_params(length=2.5, width=0.6)
cbar.ax.set_yticks(np.arange(0,6,1))

# Save figure
figure.tight_layout()
plt.savefig(figure_file, bbox_inches='tight')