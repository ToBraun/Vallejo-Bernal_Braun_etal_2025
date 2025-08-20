# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. S3 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figureS3.py
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
import matplotlib.lines as mlines
from matplotlib.patches import Patch
from cmcrameri import cm
from cartopy.util import add_cyclic_point
from plot_functions import read_csv_record, plot_landmasses, maxfactor_map

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
plot_time = [dt.datetime(2015, 11, 13, 15),
             dt.datetime(2015, 11, 14, 15),
             dt.datetime(2015, 11, 15, 15)]
central_lon = -150
color_map_ars = cm.batlow

#Files
pikart_file = root_folder + 'data/figureS2/AR_trackids_ERA5_0p5deg_6hr_2015.csv'
connect_file = root_folder + 'data/figureS3/ARCONNECT_Eulerian_MERRA2-3hr_2015.nc4'
colors_file = root_folder + 'data/figureS3/colors_PIKART.npy'
figure_file = root_folder + 'manuscript/PIKART_FigureS3.png'

# %% LOAD PIKART CATALOG

# Open csv file with AR records
ar_records = read_csv_record(pikart_file)
ar_records['time'] = pd.to_datetime(ar_records['time'])

# Extract the ARs
ars_pikart = []
ids_pikart = []
for k in range(np.size(plot_time)):
    ars = ar_records[ar_records.time == plot_time[k]]
    ids_pikart.append(np.unique(ars['trackid']))
    ars_pikart.append(ars)
    
ids_pikart = np.unique(np.concatenate(ids_pikart))

# Load AR colors
colors_pikart = np.load(colors_file)

# %% MANUALLY DEFINE THE POSITIONS OF THE LABELS

pikart_label = [
    np.array([[121808, -118, -42], [121845, -140, 32], [121846, 87, -50],
              [121847, -25, -16], [121855, -63, -12], [121856, -40, 20],
              [121858, 105, 33], [121859, -170, -53], [121860, 50, -45],
              [121861, 166, -28], [121862, 160, 22]]),

    np.array([[121846, 132, -44], [121847, 0, -30], [121856, -29, 27],
              [121858, 135, 20], [121860, 68, -55], [121862, 170, 22],
              [121864, 113, 52], [121865, -158, -14], [121866, -140, 25]]),

    np.array([[121847, 0, -28], [121856, -21, 34], [121858, 140, 20], 
              [121862, -152, 32], [121865, -150, -15], [121867, -176, -27],
              [121868, 73, -35], [121869, -175, 10], [121870, -120, 13], 
              [121871, -50, 25]])]

# %% EXTRACT AR-CONNECT CATALOG

# Load AR-CONNECT catalog
connect_dataset = xr.open_dataset(connect_file)

# Load geographic coordinates
lat_connect = connect_dataset['lat'].values
lon_connect = connect_dataset['lon'].values
time_connect = pd.to_datetime(connect_dataset['time'].values)

# Extract AR footprints
shape_connect = np.zeros((np.size(plot_time), np.size(lat_connect), np.size(lon_connect)))

for k in range(np.size(plot_time)):
    shape_connect[k,:,:] = connect_dataset.sel(
        time=plot_time[k])['ar_object_labels'].values

# Add cyclic point
shape_connect, lon_connect = add_cyclic_point(shape_connect, coord=lon_connect, axis=2)

# Mask zero values
shape_connect = np.ma.filled(shape_connect, np.nan)
shape_connect[shape_connect == 0] = np.nan

# %% PREPARE THE FIGURE LAYOUT

# Create the figure
figure = plt.figure(figsize=(170/25.4,240/25.4), dpi=500)

# Create the layout
gs = gridspec.GridSpec(3, 1) #width_ratios=[1,0.1] height_ratios=[1,0.17])

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

for n in range(np.size(plot_time)):

    # Create the subplot axis
    ax = figure.add_subplot(gs[n,0], projection=proj)
    ax.set_global()

    # Plot the landmasses
    ax = plot_landmasses(ax, lw=0.6)

    ##### PLOT PIKART

    # Plot AR contours
    for j in range(len(ars_pikart[n])):
        
        # Extract AR
        ar = ars_pikart[n].iloc[j]
        color_ar = colors_pikart[ids_pikart == ar['trackid'],:]

        # Plot AR
        linestyle=(0,(3, 2)) if ar['is_relaxed'] else '-'

        # Plot contour
        polygon = mpatches.Polygon(
            list(zip(ar['contour_lon'], ar['contour_lat'])), closed=True, ec=color_ar,
            fill=False, lw=0.75, ls=linestyle, transform=ccrs.Geodetic(), zorder=10)
        ax.add_patch(polygon)

        # Plot axis
        ax.plot(ar['axis_lon'], ar['axis_lat'], color=color_ar, linewidth=0.5,
                ls=linestyle, transform=ccrs.Geodetic(), zorder=4)

        # Label with track ID
        x_lab = pikart_label[n][pikart_label[n][:,0] == ar['trackid'], 1]
        y_lab = pikart_label[n][pikart_label[n][:,0] == ar['trackid'], 2]
        ax.annotate(str(int(ar['trackid'])), xy=(x_lab[0], y_lab[0]), xytext=(0, 0),
                    color=color_ar, textcoords='offset points',
                    ha='center', va='center', transform=ccrs.Geodetic(), zorder=12)

    # Plot AR-CONNECT
    ax.pcolormesh(lon_connect, lat_connect, shape_connect[n],
                  shading='nearest', transform=ccrs.PlateCarree(), zorder=2,
                  cmap=mpl.colors.ListedColormap([190/256, 190/256, 190/256, 1]))

    # Maxfactor the map
    if n < (np.size(plot_time)-1):
        ax = maxfactor_map(ax, draw_grid=False, ll=True, bl=False)
    else:
        ax = maxfactor_map(ax, draw_grid=False, ll=True, bl=True)

    ax.set_aspect('equal')

    # Label with date
    ax.annotate(plot_time[n].strftime('%d %b %Y, %H UTC'), xy=(33,87),
                xytext=(0, 0), color='w', textcoords='offset points', ha='left',
                va='top', transform=ccrs.Geodetic(),
                bbox = {'facecolor': 'k', 'edgecolor':'none', 'pad': 3})

# Tight layout
gs.tight_layout(figure, h_pad=0, w_pad=-50)

# %% ADD LEGEND

# Insular landmass
insular = Patch(edgecolor='dimgrey', facecolor='none', hatch='\\\\\\\\',
                linewidth=0.6, linestyle='solid')

# Continental landmass
continental = Patch(edgecolor='dimgrey', facecolor='none', linewidth=0.6,
                    linestyle='solid')

# Strict AR
color_edge = colors_pikart[0,:]
strict_ar = Patch(edgecolor=color_edge, facecolor='none', linewidth=0.75,
                    linestyle='solid')
strict_ar_line = mlines.Line2D([], [], linewidth=0.5, color=color_edge)

# Relaxed AR
color_edge = colors_pikart[6,:]
relaxed_ar = Patch(edgecolor=color_edge, facecolor='none', linewidth=0.75,
                    linestyle=(0,(3, 2)))
relaxed_ar_line = mlines.Line2D([], [], linewidth=0.5, color=color_edge,
                                linestyle=(0,(3, 2)))

# AR-CONNECT
color_face = [190/256, 190/256, 190/256, 1]
connect_ar = Patch(edgecolor='none', facecolor=color_face)

# Add legend
plt.legend([continental, (strict_ar, strict_ar_line), (relaxed_ar, relaxed_ar_line),
            insular, connect_ar],
            ['Continental landmass', 'PIKART: AR shape', 'PIKART: relaxed AR shape',
             'Insular landmass', 'AR-CONNECT'], fontsize=7,
            bbox_to_anchor=(0.23,0.02), loc="lower center", ncol=2,
            columnspacing=1, edgecolor='none', framealpha=1)

# Save figure
figure.tight_layout()
plt.savefig(figure_file, bbox_inches='tight')
