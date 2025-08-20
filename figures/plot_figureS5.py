# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. S5 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figureS5.py
'''

# %% IMPORT MODULES

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import datetime as dt
import xarray as xr
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cmcrameri import cm
from plot_functions import read_csv_record, plot_landmasses, maxfactor_map_figS5

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

# Parameters
central_lon = 90
color_ar = cm.roma(220)

# Files
ivt_file = root_folder + 'data/figureS5/ivt_ERA5_0p75_6hr_2016_SOND.nc'
ipart_file = root_folder + 'data/figureS5/IPART_Lagrangian_ERA5_0p75_6hr_2016.csv'
figure_file = root_folder + 'manuscript/PIKART_FigureS5.png'

# %% LOAD DATA

# Open csv file with AR records
ardf = read_csv_record(ipart_file)
ardf['time'] = pd.to_datetime(ardf['time'])

# Load IVT dataset
ivt_dataset = xr.open_dataset(ivt_file)

# Load geographic coordinates
lat = ivt_dataset.coords['lat'].values
lon = ivt_dataset.coords['lon'].values
time = pd.to_datetime(ivt_dataset.coords['time'].values)

# IVT colormap
step = int(np.floor(120/10))
num = int(9)
colors = np.stack(cm.grayC_r.colors, axis=0)
colors = colors[np.arange(0, num) * step, :]
cmap = (mpl.colors.ListedColormap(colors[:-1,:]).with_extremes(over=colors[-1,:]))
bounds = np.arange(0, 510, 125/2)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# %% PLOT DETECTION

# Create the figure
figure = plt.figure(figsize=(180/25.4,60/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# Create the layout
gs = gridspec.GridSpec(1, 3, width_ratios=[0.6,0.6,0.1]) # height_ratios=[1,0.17])

# %% PANEL A

# Panel time
plot_time = dt.datetime(2016, 9, 8, 9)
ar_id = 14

# Create the subplot axis
ax = figure.add_subplot(gs[0], projection=proj)

# Set the extension of the map
ax.set_extent([50, 150, 30, 90], crs=ccrs.PlateCarree())

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.6)

# Extract the IVT
ivt = ivt_dataset.sel(time=plot_time)['ivt'].values

# Plot the IVT field
cs = ax.pcolormesh(lon, lat, ivt, shading='nearest', cmap=cmap, norm=norm, 
                   transform=ccrs.PlateCarree(), zorder=0)

#Extract the ARs
ar = ardf[(ardf['time'] == plot_time) & (ardf['id'] == ar_id)].iloc[0]

# Plot AR
linestyle=(0,(3, 2)) if ar['is_relaxed'] else '-'

# Plot contour
polygon = mpatches.Polygon(
    list(zip(ar['contour_x'], ar['contour_y'])), closed=True, ec=color_ar,
    fill=False, lw=1, ls=linestyle, transform=ccrs.Geodetic())
ax.add_patch(polygon)

# Plot axis
ax.plot(ar['axis_x'], ar['axis_y'], color=color_ar, linewidth=0.65,
        alpha=0.75, ls=linestyle, transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map_figS5(ax)
ax.set_aspect('equal')
title = '%s' %(plot_time.strftime('%d %b %Y, %H UTC'))
ax.set_title(title)
ax.text(0, 1.045, '(a)', ha='left', transform=ax.transAxes)

# %% PANEL B

# Panel time
plot_time = dt.datetime(2016, 12, 12, 21)
ar_id = 12

# Create the subplot axis
ax = figure.add_subplot(gs[1], projection=proj)

# Set the extension of the map
ax.set_extent([50, 150, 30, 90], crs=ccrs.PlateCarree())

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.6)

# Extract the IVT
ivt = ivt_dataset.sel(time=plot_time)['ivt'].values

# Plot the IVT field
cs = ax.pcolormesh(lon, lat, ivt, shading='nearest', cmap=cmap, norm=norm, 
                   transform=ccrs.PlateCarree(), zorder=0)

#Extract the ARs
ar = ardf[(ardf['time'] == plot_time) & (ardf['id'] == ar_id)].iloc[0]

# Plot AR
linestyle=(0,(3, 2)) if ar['is_relaxed'] else '-'

# Plot contour
polygon = mpatches.Polygon(
    list(zip(ar['contour_x'], ar['contour_y'])), closed=True, ec=color_ar,
    fill=False, lw=1, ls=linestyle, transform=ccrs.Geodetic())
ax.add_patch(polygon)

# Plot axis
ax.plot(ar['axis_x'], ar['axis_y'], color=color_ar, linewidth=0.65,
        alpha=0.75, ls=linestyle, transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map_figS5(ax, ll=False)
ax.set_aspect('equal')
title = '%s' %(plot_time.strftime('%d %b %Y, %H UTC'))
ax.set_title(title)
ax.text(0, 1.045, '(b)', ha='left', transform=ax.transAxes)

# %% ADD COLORBAR

# Add subplot to place the colorbar
ax = figure.add_subplot(gs[2])
ax.set_axis_off()
cbaxes = inset_axes(ax, width=0.11, height=1.77, loc='center left',
                    bbox_to_anchor=(0, 0.5),
                    bbox_transform=ax.transData,
                    borderpad=0) #

cbar = plt.colorbar(
    cs, cax=cbaxes, extend='max', orientation='vertical', spacing='proportional',
    ticks=bounds)
cbar.set_label('IVT (kg m$^{-1}$ s$^{-1}$)', labelpad=4)
cbar.ax.tick_params(length=2.5, width=0.6)
cbar.ax.set_yticklabels(['0', '', '125', '', '250', '', '375', '', '500'])

# Save figure
figure.tight_layout()
#gs.tight_layout(figure, w_pad=1)
plt.savefig(figure_file)
