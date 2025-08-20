# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. S6 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figureS6.py
'''

# %% IMPORT MODULES

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import datetime as dt
import xarray as xr
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.patches import Patch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cmcrameri import cm
from plot_functions import read_csv_record, plot_landmasses, maxfactor_map
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

# Parameters
plot_time = dt.datetime(2015, 11, 27, 21)
central_lon = -150
color_map = cm.batlow

# Files
ivt_file = root_folder + 'data/figureS6/ivt_ERA5_0p75_6hr_2015_Nov.nc'
ipart_file = root_folder + 'data/figureS6/IPART_Lagrangian_ERA5_0p75_6hr_2015.csv'
figure_file = root_folder + 'manuscript/PIKART_FigureS6.png'

# Markers
cent_marker = 'o'
cent_marker_style = dict(color='r', markersize=2.5, fillstyle='full', markeredgewidth=0)

miss_cent_marker = 'X'
miss_cent_marker_style = dict(color='m', markersize=5, fillstyle='full', markeredgewidth=0)

# %% LOAD DATA

# Open csv file with AR records
ardf = read_csv_record(ipart_file)
ardf['time'] = pd.to_datetime(ardf['time'])

# Load IVT dataset
ivt_dataset = xr.open_dataset(ivt_file)

# Load geographic coordinates
lat = ivt_dataset['lat'].values
lon = ivt_dataset['lon'].values
time = pd.to_datetime(ivt_dataset['time'].values)

# Extract the IVT
ivt = ivt_dataset.sel(time=plot_time)['ivt'].values

# Add cyclic point to IVT
ivt, lon = add_cyclic_point(ivt, coord=lon, axis=1)

# %% PLOT DETECTION

# Create the figure
figure = plt.figure(figsize=(180/25.4,107/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# Create the subplot axis
ax = figure.add_subplot(1,1,1, projection=proj)
ax.set_global()

# Plot the landmasses
ax = plot_landmasses(ax)

# IVT colormap
step = int(np.floor(120/14))
num = int(13)
colors = np.stack(cm.grayC_r.colors, axis=0)
colors = colors[np.arange(0, num) * step, :]
cmap = (mpl.colors.ListedColormap(colors[:-1,:]).with_extremes(over=colors[-1,:]))
bounds = np.arange(0, 780, 125/2)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# Plot the IVT field
cs = ax.pcolormesh(lon, lat, ivt, shading='nearest', cmap=cmap, norm=norm, 
                   transform=ccrs.PlateCarree(), zorder=0)

# Extract the ARs
plot_ars = ardf[ardf.time == plot_time]

# Plot AR contours
for j in range(len(plot_ars)):

    # Extract AR
    ar = plot_ars.iloc[j]
    color_ar = color_map(35*j)
    color_face = color_ar[:-1] + (0.3,)

    # Mask IVT
    polygon = mpatches.Polygon(list(zip(ar['contour_x'], ar['contour_y'])),
                                closed=True, ec='w', fill=True, lw=0,
                                fc='w', transform=ccrs.Geodetic())
    ax.add_patch(polygon)

    # Plot AR
    linestyle=(0,(3, 2)) if ar['is_relaxed'] else '-'

    # Plot contour
    polygon = mpatches.Polygon(
        list(zip(ar['contour_x'], ar['contour_y'])), closed=True, ec=color_ar,
        fill=True, lw=1, fc=color_face, ls=linestyle, transform=ccrs.Geodetic())
    ax.add_patch(polygon)

    # Plot axis
    ax.plot(ar['axis_x'], ar['axis_y'], color=color_ar, linewidth=0.65,
            alpha=0.75, ls=linestyle, transform=ccrs.Geodetic())

    # Check if the AR crosses the antimeridian
    if np.any(np.abs(np.diff(ar['contour_x'])) > 350):

        # Plot missplaced centroid
        ax.plot(ar['centroid_x'], ar['centroid_y'], marker=miss_cent_marker,
                transform=ccrs.Geodetic(), **miss_cent_marker_style)

    else:

        # Plot missplaced centroid
        ax.plot(ar['centroid_x'], ar['centroid_y'], marker=cent_marker,
                transform=ccrs.Geodetic(), **cent_marker_style)

# Plot antimeridian
ax.plot([180, 180], [-90, 90], ls='--', lw=0.5, color='k', dashes=(5, 5),
        transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map(ax, draw_grid=False)
ax.set_aspect('equal')
title = 'AR identification — %s' %(plot_time.strftime('%d %b %Y, %H UTC'))
ax.set_title(title)

# %% ADD LEGEND

# Insular landmass
insular = Patch(edgecolor='dimgrey', facecolor='none', hatch='\\\\\\\\',
                linewidth=0.6, linestyle='solid')

# Continental landmass
continental = Patch(edgecolor='dimgrey', facecolor='none', linewidth=0.6,
                    linestyle='solid')

# Strict AR
color_edge = color_map(0)
color_face = color_edge[:-1] + (0.2,)
strict_ar = Patch(edgecolor=color_edge, facecolor=color_face, linewidth=1,
                    linestyle='solid')
strict_ar_line = mlines.Line2D([], [], linewidth=0.65, color=color_edge, alpha=0.75)

# Markers
centroid = mlines.Line2D([], [], marker=cent_marker, linestyle='None',
                         **cent_marker_style)
miss_centroid = mlines.Line2D([], [], marker=miss_cent_marker, linestyle='None',
                              **miss_cent_marker_style)

# Add legend
plt.legend([(strict_ar, strict_ar_line), continental, centroid, insular, 
            miss_centroid],
           ['AR shape', 'Continental landmass', 'Correctly positioned centroid',
            'Insular landmass', 'Incorrectly positioned centroid'],
           bbox_to_anchor=(0.5,-0.2), loc="lower center", ncol=3,
           edgecolor='none', framealpha=1)

figure.tight_layout()

# Add colorbar
cbaxis = inset_axes(ax, width=1.2, height=0.08, loc=3,
                    bbox_to_anchor=(-173, -86), bbox_transform=ax.transData)
cbar = figure.colorbar(cs, cax=cbaxis, orientation='horizontal',
                    spacing='proportional', extend='max')
cbar.set_ticks([0, 125, 250, 375, 500, 625, 750])
cbar.ax.set_xticklabels(['0', '', '250', '', '500', '', '750'])
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.tick_params(length=2.5, width=0.6, pad=1.5)
cbar.set_label('IVT (kg m$^{-1}$ s$^{-1}$)', labelpad=4)
cbar.ax.xaxis.set_label_position('top')

# Create a rectangle patch
rect = Rectangle((46, -75), 52, 12, edgecolor='none', facecolor='w', 
                 transform=ccrs.PlateCarree(), zorder=5)

# Add the patch to the Axes
ax.add_patch(rect)

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
