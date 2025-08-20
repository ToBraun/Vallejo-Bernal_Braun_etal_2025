# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. S4 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figureS4.py
'''

# %% IMPORT MODULES

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import datetime as dt
import xarray as xr
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from cmcrameri import cm
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
central_lon = -150
color_map = cm.batlow

# Files
eulerian_file = root_folder + 'data/figureS4/PIKARTV1_eulerian_ERA5_0p5deg_6hr_1979.nc'
lagrangian_file = root_folder + 'data/figureS4/AR_footprints_ERA5_0p5deg_6hr_1979.nc'
relaxed_file = root_folder + 'data/figureS4/AR_trackids_ERA5_0p5deg_6hr_1979.csv'
figure_file = root_folder + 'manuscript/PIKART_FigureS4.png'

# Track IDs of ARs to plot
# Short, long, big, round, incoherent IVT direction
relaxed_ars = [[63515, dt.datetime(1979, 1, 18, 3), 7, -9, -45],
               [63440, dt.datetime(1979, 1, 22, 9), 10, -110, -25],
               [64596, dt.datetime(1979, 10, 4, 21), 10, -170, 28],
               [63622, dt.datetime(1979, 2, 12, 3), 10, 75, -50],
               [64101, dt.datetime(1979, 6, 2, 15), 10, -43, 23]]

labels = ['Short', 'Long', 'Big', 'Round', 'Incoherent IVT direction']

# %% LOAD DATA

# Open csv file with AR records
ar_records = read_csv_record(relaxed_file)
ar_records['time'] = pd.to_datetime(ar_records['time'])

# %% PLOT RELAXED ARs

# Create the figure
figure = plt.figure(figsize=(180/25.4,107/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# Create the subplot axis
ax = figure.add_subplot(1,1,1, projection=proj)
ax.set_global()

# Plot the landmasses
ax = plot_landmasses(ax)

# Plot AR contours
for j in range(len(relaxed_ars)):

    # Extract AR
    ar = ar_records.loc[
        (ar_records['trackid'] == relaxed_ars[j][0]) &
        (ar_records['time'] == relaxed_ars[j][1])].iloc[0]
    color_ar = color_map(52*j)
    color_face = color_ar[:-1] + (0.2,)

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

    # LOAD AND EXTRACT IVT

    # Load PIKART datasets
    eulerian_dataset = xr.open_dataset(eulerian_file)
    lagrangian_dataset = xr.open_dataset(lagrangian_file)

    # Load geographic coordinates
    lat = eulerian_dataset.coords['latitude'].values
    lon = eulerian_dataset.coords['longitude'].values
    time = pd.to_datetime(eulerian_dataset.coords['time'].values)

    # Get IVT data
    ivtu = eulerian_dataset.sel(time=relaxed_ars[j][1])['ivtu'].values
    ivtv = eulerian_dataset.sel(time=relaxed_ars[j][1])['ivtv'].values

    # Get AR footprints
    foot_pikart = lagrangian_dataset.sel(time=relaxed_ars[j][1])['footprint'].values

    # Mask IVT fields
    ivtu[foot_pikart != ar['trackid']] = np.nan
    ivtv[foot_pikart != ar['trackid']] = np.nan

    # Plot IVT field
    arrow_jump = relaxed_ars[j][2]
    q = ax.quiver(
        lon[::arrow_jump], lat[::arrow_jump],
        ivtu[::arrow_jump, ::arrow_jump], ivtv[::arrow_jump, ::arrow_jump],
        pivot='mid', color='#262626', transform=ccrs.PlateCarree(),
        scale_units='xy', scale=150, width=.0013,
        headwidth=4, headlength=3.5, headaxislength=3, zorder=5)
    
    # Add the track ID
    ax.annotate(relaxed_ars[j][1].strftime('%d %b %Y, %H UTC'),
                xy=(relaxed_ars[j][3], relaxed_ars[j][4]), xytext=(0, 0),
                color='k', textcoords='offset points',
                ha='center', va='center', fontsize=7, transform=ccrs.Geodetic(),
                bbox = {'facecolor': color_face, 'edgecolor': 'none', 'pad': 1})

# Maxfactor the map
ax = maxfactor_map(ax, draw_grid=False)
ax.set_aspect('equal')

# Legend elements
def leg_elem(color):

    color_edge = color
    color_face = color_edge[:-1] + (0.2,)
    patch = Patch(edgecolor=color_edge, facecolor=color_face, linewidth=1,
                  ls=(0,(3, 2)))
    line = mlines.Line2D([], [], linewidth=0.65, color=color_edge, alpha=0.75,
                         ls=(0,(3, 2)))

    return (patch, line)

# Add legend
plt.legend([leg_elem(color_map(52*0)), leg_elem(color_map(52*1)),
            leg_elem(color_map(52*2)), leg_elem(color_map(52*3)),
            leg_elem(color_map(52*4))],
           labels, bbox_to_anchor=(0.41,-0.16),
           loc="lower center", ncol=5, edgecolor='none', framealpha=1)

ax.quiverkey(q, X=0.81, Y=-0.112, U=1000, label = '1000 kg m$^{-1}$ s$^{-1}$',
             labelpos='E', coordinates='axes', labelsep=0.05)

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
