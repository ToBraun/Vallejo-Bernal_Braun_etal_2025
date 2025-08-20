# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. 2 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure2.py
'''

# %% IMPORT MODULES

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import datetime as dt
import xarray as xr
from cmcrameri import cm
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.patches import Patch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from plot_functions import plot_landmasses, maxfactor_map, read_csv_record
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
plot_time = dt.datetime(2019, 1, 23, 21)
central_lon = -150
color_map = cm.batlow
continents = ['na', 'sa', 'eu', 'af', 'as', 'au', 'oc', 'an']

# Files
ivt_file = root_folder + 'data/figure02/ivt_ERA5_0p5deg_6hr_2019-01-23.nc'
pikart_file = root_folder + 'data/figure02/AR_land_intersection_ERA5_0p5deg_6hr_2019.csv'
figure_file = root_folder + 'manuscript/PIKART_Figure2.png'

# Markers
cent_marker = 'o'
cent_marker_style = dict(
    color='r', markersize=2.5, fillstyle='full', markeredgewidth=0)

lfc_marker = 's'
lfc_marker_style = dict(
    color='m', markersize=3, fillstyle='full', markeredgewidth=0)

lfi_marker = '^'
lfi_marker_style = dict(
    color='aqua', markersize=4, fillstyle='full', markeredgewidth=0)

# %% LOAD DATA

# Open csv file with AR records
ar_records = read_csv_record(pikart_file)
ar_records['time'] = pd.to_datetime(ar_records['time'])

# Load IVT data
ivt = xr.open_dataset(ivt_file)

# Load geographic coordinates
lat = ivt['lat'].values  
lon = ivt['lon'].values

# Extract the IVT
ivt = ivt.sel(time=plot_time)['ivt'].values

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
cs = ax.pcolormesh(lon, lat, ivt, shading='nearest',
                   cmap=cmap, norm=norm, transform=ccrs.PlateCarree(),
                   zorder=0)

# Extract the ARs
ardf = ar_records[ar_records.time == plot_time]

# Plot AR contours
for j in range(len(ardf)):

    # Extract AR
    ar = ardf.iloc[j]
    color_ar = color_map(16*j)
    color_face = color_ar[:-1] + (0.3,)

    # Mask IVT
    polygon = mpatches.Polygon(list(zip(ar['contour_lon'], ar['contour_lat'])),
                                closed=True, ec='w', fill=True, lw=0,
                                fc='w', transform=ccrs.Geodetic())
    ax.add_patch(polygon)

    # Plot AR
    linestyle=(0,(3, 2)) if ar['is_relaxed'] else '-'

    # Plot contour
    polygon = mpatches.Polygon(
        list(zip(ar['contour_lon'], ar['contour_lat'])), closed=True, ec=color_ar,
        fill=True, lw=1, fc=color_face, ls=linestyle, transform=ccrs.Geodetic())
    ax.add_patch(polygon)

    # Plot axis
    ax.plot(ar['axis_lon'], ar['axis_lat'], color=color_ar, linewidth=0.65,
            alpha=0.75, ls=linestyle, transform=ccrs.Geodetic())

# Plot centroids
centroid_coord = np.array(ardf[['centroid_lon', 'centroid_lat']])

for cent in centroid_coord:
    ax.plot(cent[0], cent[1], marker=cent_marker, transform=ccrs.Geodetic(),
            **cent_marker_style)

# Get landfall locations
lf_coord_cont = []
lf_coord_ins = []
for k in continents:
    lf_coord_cont.append(
        (ardf[ardf[['conti_%s_lon' %(k), 'conti_%s_lat' %(k)]].notnull().all(1)][
            ['conti_%s_lon' %(k), 'conti_%s_lat' %(k)]]).values.tolist())
    lf_coord_ins.append(
        (ardf[ardf[['insu_%s_lon' %(k), 'insu_%s_lat' %(k)]].notnull().all(1)][
            ['insu_%s_lon' %(k), 'insu_%s_lat' %(k)]]).values.tolist())

lf_coord_cont = [item for sublist in lf_coord_cont for item in sublist]
lf_coord_ins = [item for sublist in lf_coord_ins for item in sublist]

# Plot landfall locations
for lf in lf_coord_cont:
    ax.plot(lf[0], lf[1], marker=lfc_marker, transform=ccrs.Geodetic(), 
            **lfc_marker_style)

for lf in lf_coord_ins:
    ax.plot(lf[0], lf[1], marker=lfi_marker, transform=ccrs.Geodetic(), 
            **lfi_marker_style)

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
strict_ar_line = mlines.Line2D(
    [], [], linewidth=0.65, color=color_edge, alpha=0.75)

# Relaxed AR
color_edge = color_map(16)
color_face = color_edge[:-1] + (0.2,)
relaxed_ar = Patch(edgecolor=color_edge, facecolor=color_face, linewidth=1,
                    linestyle=(0,(3, 2)))
relaxed_ar_line = mlines.Line2D([], [], linewidth=0.65, color=color_edge,
                                alpha=0.75, ls=(0,(3, 2)))

# Markers
centroid = mlines.Line2D(
    [], [], marker=cent_marker, linestyle='None', **cent_marker_style)
lf_cont = mlines.Line2D(
    [], [], marker=lfc_marker, linestyle='None', **lfc_marker_style)
lf_ins = mlines.Line2D(
    [], [], marker=lfi_marker, linestyle='None', **lfi_marker_style)

# Add legend
plt.legend([continental, insular, (strict_ar, strict_ar_line),
            (relaxed_ar, relaxed_ar_line), lf_cont, lf_ins, centroid],
            ['Continental landmass', 'Insular landmass', 'AR shape',
            'Relaxed AR shape', 'Continental intersection',
            'Insular intersection', 'AR centroid'], bbox_to_anchor=(0.5,-0.21),
            loc="lower center", ncol=4, edgecolor='none', framealpha=1)

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
rect = Rectangle((46, -73), 52, 10, edgecolor='none', facecolor='w', 
                 transform=ccrs.PlateCarree(), zorder=5)

# Add the patch to the Axes
ax.add_patch(rect)

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
