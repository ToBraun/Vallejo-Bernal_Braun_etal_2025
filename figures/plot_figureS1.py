# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. S1 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figureS1.py
'''

# %% IMPORT MODULES

import geopandas as gpd
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from cmcrameri import cm
from plot_functions import maxfactor_map

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

# Parameters
central_lon = -150
color_map = cm.managua
colors = [29,86,199,114,58,143,228,171]
continents = ['North America', 'South America', 'Europe', 'Africa',
               'Asia', 'Australia', 'Oceania', 'Antarctica']

# Files
shape_cont_file = root_folder + 'data/world_continents/world_continental_3000sqkm.shp'
shape_ins_file = root_folder + 'data/world_continents/world_insular_3000sqkm.shp'
figure_file = root_folder + 'manuscript/PIKART_FigureS1.png'

# Legend elements
insular = Patch(edgecolor='dimgrey', facecolor='none', hatch='\\\\\\\\\\',
                linewidth=0.6, linestyle='solid')

continental = Patch(edgecolor='dimgrey', facecolor='none', linewidth=0.6,
                    linestyle='solid')

# %% LOAD DATA

# Load continental and insular shapefiles
shape_cont = gpd.read_file(shape_cont_file)
shape_ins = gpd.read_file(shape_ins_file)

# %% PLOT LANDMASSES

# Create the figure
figure = plt.figure(figsize=(180/25.4,102/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# Create the subplot axis
ax = figure.add_subplot(1,1,1, projection=proj)
ax.set_global()

# Change hatch linewidth
mpl.rcParams.update({'hatch.linewidth': 0.6})

# Plot the landmasses
handles = []
for i in range(len(continents)):

    # Plot the islands
    insdf = shape_ins[shape_ins['CONTINENT'] == continents[i]]
    ax.add_geometries(insdf['geometry'], crs=ccrs.PlateCarree(), linewidth=0.6,
                      linestyle='solid', edgecolor=color_map(colors[i]), 
                      facecolor='none', hatch='\\\\\\\\\\')

    # Plot the continents
    contdf = shape_cont[shape_cont['CONTINENT'] == continents[i]]
    ax.add_geometries(contdf['geometry'], crs=ccrs.PlateCarree(), linewidth=0.6,
                      linestyle='solid', edgecolor=color_map(colors[i]), 
                      facecolor=color_map(colors[i]))

    handles.append(Patch(edgecolor=color_map(colors[i]), 
                         facecolor=color_map(colors[i]),
                         linewidth=0.6, linestyle='solid'))

# Maxfactor the map
ax = maxfactor_map(ax, draw_grid=False)
ax.set_aspect('equal')

# Legend
handles.append(continental)
handles.append(insular)
continents.append('Continental landmass')
continents.append('Insular landmass')
plt.legend(handles, continents, bbox_to_anchor =(0.5,-0.2), loc="lower center",
           ncol=5, edgecolor='none', framealpha=1)

# Save figure
figure.tight_layout()
plt.savefig(figure_file, bbox_inches='tight')