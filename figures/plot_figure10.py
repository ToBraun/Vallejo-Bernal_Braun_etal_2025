# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. 10 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure10.py
'''

# %% IMPORT MODULES

import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
from cmcrameri import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_functions import plot_landmasses, maxfactor_fig10
import matplotlib.gridspec as gridspec
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
years = 2023 - 1979 + 1
central_lon = -150
color_map = cm.bilbao_r
contours=[13, 10, 12, 8, 7]
ax_labels = [False, False, False, False, True] 

# Files
rank_file = root_folder + 'data/figure10/rank_count.npy'
figure_file = root_folder + 'manuscript/PIKART_Figure10.png'

# %% LOAD DATA

rank_count = np.load(rank_file)/years
longitude = np.arange(-180, 180, .5)
latitude = np.arange(-89.75, 90, .5)

# %% PLOT

# Create the figure
figure = plt.figure(figsize=(115/25.4,240/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# Create the layout
gs = gridspec.GridSpec(5, 1)

literal = ['(a)', '(b)', '(c)', '(d)', '(e)']

for rank in range(5):
    
    # Extract the data    
    plot_map = np.where(rank_count[rank,:,:]>0, rank_count[rank,:,:], np.nan)#/24
        
    # Create the subplot axis
    ax = figure.add_subplot(gs[rank], projection=proj)
    ax.set_global()
    
    # Plot the landmasses
    ax = plot_landmasses(ax, lw=0.5, ec='gray')
    
    # Plot the AR rank frequency
    plot_map, longitude_plot = add_cyclic_point(plot_map, coord=longitude, axis=1)
    
    map_rank = ax.contourf(longitude_plot, latitude, plot_map, cmap=color_map, 
                           vmin=0, vmax=np.ceil(np.nanmax(plot_map)), 
                           levels=contours[rank], 
                           transform=ccrs.PlateCarree(), zorder=0)
        
    # Maxfactor the map
    ax = maxfactor_fig10(ax, ll=True, bl=ax_labels[rank])
    ax.set_aspect('equal')
    ax.set_ylabel('Rank %s' %(rank+1))
        
    # Color bar
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size='3%', pad=0.1, axes_class=plt.Axes)
    figure.add_axes(ax_cb)
    cbar = plt.colorbar(map_rank, cax=ax_cb)
    cbar.ax.tick_params(length=2, width=0.6)
    ax.text(1.15, 0.5, 'AR events/year', rotation=90, va='center', ha='center',
            transform=ax.transAxes)
    
    # Literal
    ax.text(-0.015, 0.95, literal[rank], ha='right', transform=ax.transAxes)
    
# Save figure
gs.tight_layout(figure, h_pad=0.5)
plt.savefig(figure_file, bbox_inches='tight')
