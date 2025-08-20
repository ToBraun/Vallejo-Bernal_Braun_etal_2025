# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. 6 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure6.py
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
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib import image
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from plot_functions import plot_landmasses, maxfactor_map_fig6, read_csv_record

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
trackid = 124915
plot_time = dt.datetime(2017, 10, 14, 15)
plot_time = pd.to_datetime(plot_time)
central_lon = 185
arrow_jump = 6
color_map = cm.grayC_r

# Files and folders
photo_file = root_folder + 'data/figure06/ar_viirs_14-10-2017.jpg'
ivt_file = root_folder + 'data/figure06/ivt_ERA5_0p5deg_6hr_2017-10-14.nc'
ivtu_file = root_folder + 'data/figure06/ivtu_ERA5_0p5deg_6hr_2017-10-14.nc'
ivtv_file = root_folder + 'data/figure06/ivtv_ERA5_0p5deg_6hr_2017-10-14.nc'
pikart_file = root_folder + 'data/figure06/AR_land_intersection_ERA5_0p5deg_6hr_2017.csv'
pikart_euler_file = root_folder + 'data/figure06/AR_footprints_ERA5_0p5deg_6hr_2017-10-14.nc'
target_file_shape = root_folder + 'data/figure06/tARgetv4_shapemap_2017.nc'
target_file_axis = root_folder + 'data/figure06/tARgetv4_axismap_2017.nc'
figure_file = root_folder + 'manuscript/PIKART_Figure6.png'

# %% LOAD AND EXTRACT IVT

# Load ivt data
ivt = xr.open_dataset(ivt_file) 
ivtu = xr.open_dataset(ivtu_file) 
ivtv = xr.open_dataset(ivtv_file) 

# Load geographic coordinates
lat = ivt['lat'].values  
lon = ivt['lon'].values

# Extract the IVT
ivt = ivt.sel(time=plot_time)['ivt'].values
ivtu = ivtu.sel(time=plot_time)['ivtu'].values
ivtv = ivtv.sel(time=plot_time)['ivtv'].values

# %% LOAD AN EXTRACT PIKART AR

# Open csv file with PIKART AR records
ardf = read_csv_record(pikart_file)
ardf['time'] = pd.to_datetime(ardf['time'])

# Extract the ARs
ardf = ardf[ardf['trackid'] == trackid]
ar_pikart = ardf[ardf['time'] == plot_time].iloc[0]

# Load nc file with AR footprint
mask_pikart = xr.open_dataset(pikart_euler_file)

# Extract the AR footprint (with the track ID to correctly mask the IVT fields)
mask_pikart = mask_pikart.sel(time=plot_time)['footprint'].values
mask_pikart = np.where(mask_pikart == trackid, 1, np.nan)

# Mask the variables
ivtu = ivtu * mask_pikart
ivtv = ivtv * mask_pikart

# %% EXTRACT tARget AR

# Load tARget3 catalog
shape_target_dataset = xr.open_dataset(target_file_shape)
axis_target_dataset = xr.open_dataset(target_file_axis)

# Load geographic coordinates
lat_target = shape_target_dataset['lat'].values
lon_target = shape_target_dataset['lon'].values

# Extract AR footprint
shape_target = shape_target_dataset.sel(
    time=(plot_time - dt.timedelta(hours=3)))['shapemap'].values
shape_target = (shape_target == 2).astype('int')

# Extract AR contour
cs = plt.contour(lon_target, lat_target, shape_target, [0.9,1.1])
plt.close()  # Prevent plot from popping up
conts = cs.collections[0].get_paths()
conts.sort(key=lambda x:len(x.vertices))
cont_target = conts[-1]
cont_target = cont_target.vertices

# Extract the AR axis
axis_target = np.squeeze(axis_target_dataset.sel(
    time=(plot_time - dt.timedelta(hours=3)))['axismap'].values)
axis_target = np.where((axis_target >= 2) & (axis_target < 3), 1, np.nan) 

# %% PREPARE THE FIGURE LAYOUT

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude = central_lon)

# Create the figure
figure = plt.figure(figsize=(180/25.4, 128/25.4), dpi=500)

# Create the layout
gs = gridspec.GridSpec(2, 2, width_ratios=[1,0.1])

# %% PLOT PHOTO

# Create the subplot axis
ax = figure.add_subplot(gs[0,0])

# Open JPG image
photo14 = image.imread(photo_file)

# Plot the photo
plt.imshow(photo14)
plt.xticks([])
plt.yticks([])

# Maxfactor the map
ax.set_aspect('equal')
title = 'Satellite image — %s' %(plot_time.strftime('%d %b %Y'))
ax.set_title(title)
ax.text(0, 1.04, '(a)', ha='left', transform=ax.transAxes)

# %% PLOT AR DETECTION

# Create the subplot axis
ax = figure.add_subplot(gs[1,0], projection=proj)

# Set the extension of the map
ax.set_extent([105, 265, 13, 67], crs=ccrs.PlateCarree())

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.75, ec='gray')

# Colorbar settings
step = int(np.floor(250/11))
num = int(10)
colors = np.stack(color_map.colors, axis=0)
colors = colors[np.arange(0, num) * step, :]
cmap = (mpl.colors.ListedColormap(colors[:-1,:]).with_extremes(over=colors[-1,:]))
bounds = [0, 100, 250, 375, 500, 625, 750, 1000, 1250, 1500]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# Plot the IVT field
cs = ax.pcolormesh(lon, lat, ivt, shading='nearest',
                   cmap=cmap, norm=norm, transform=ccrs.PlateCarree(),
                   zorder=0)

# Plot the arrows
q = ax.quiver(
    lon[::arrow_jump], lat[::arrow_jump],
    ivtu[::arrow_jump, ::arrow_jump], ivtv[::arrow_jump, ::arrow_jump],
    pivot='mid', color='#262626', transform=ccrs.PlateCarree(),
    scale_units='xy', scale=200, width=.0013,
    headwidth=4, headlength=6.5, headaxislength=5.3)

# Plot contour tARget
color_target = [0.8899, 0.592087, 0.330454]
lw_target = 0.9
polygon = mpatches.Polygon(list(zip(cont_target[:,0], cont_target[:,1])),
                            closed=True, ec=color_target, fill=False,
                            lw=lw_target, transform=ccrs.Geodetic())
ax.add_patch(polygon)

#Plot axis
ax.pcolormesh(lon_target, lat_target, axis_target, shading='nearest',
              cmap=mpl.colors.ListedColormap(color_target),
              transform=ccrs.PlateCarree(), zorder=3)

# Plot contour PIKART
color_pikart = [0.118992, 0.362849, 0.382713]
lw_pikart = 1.25
polygon = mpatches.Polygon(
    list(zip(ar_pikart['contour_lon'], ar_pikart['contour_lat'])), closed=True,
    ec=color_pikart, fill=False, lw=lw_pikart, transform=ccrs.Geodetic())
ax.add_patch(polygon)

#Plot axis
lw_axis_pikart = 0.75
ax.plot(ar_pikart['axis_lon'], ar_pikart['axis_lat'], color=color_pikart,
        lw=lw_axis_pikart, transform=ccrs.Geodetic(), zorder=5)

# Maxfactor the map
ax = maxfactor_map_fig6(ax, draw_grid=False)
ax.set_aspect('equal')
title = 'PIKART identification — %s' %(plot_time.strftime('%d %b %Y, %H UTC'))
ax.set_title(title)
ax.text(0, 1.04, '(b)', ha='left', transform=ax.transAxes)

# Legend markers
pikart = Patch(edgecolor=color_pikart, facecolor='none', lw=lw_pikart)
pikart_line = mlines.Line2D([], [], lw=lw_axis_pikart, color=color_pikart)

target = Patch(edgecolor=color_target, facecolor='none', lw=lw_target)
target_line = mlines.Line2D([], [], lw=0.65, color=color_target)

# Add legend
plt.legend(
    [(pikart, pikart_line), (target, target_line)], ['PIKART', 'tARget-4'],
    bbox_to_anchor=(0.01, 0.895), loc='upper left', edgecolor='none', facecolor='w')

ax.quiverkey(q, X=0.05, Y=0.91, U=750, label = '750 kg m$^{-1}$ s$^{-1}$',
             labelpos='E', coordinates='axes', labelsep=0.05)

# Add subplot to place the colorbar
ax = figure.add_subplot(gs[1,1])
ax.set_axis_off()
cbaxes = inset_axes(ax, width=0.1, height=2, loc='center left',
                    bbox_to_anchor=(0, 0.5),
                    bbox_transform=ax.transData,
                    borderpad=0) #

cbar = plt.colorbar(
    cs, cax=cbaxes, extend='max', orientation='vertical', spacing='proportional',
    ticks=bounds)
cbar.set_label('IVT (kg m$^{-1}$ s$^{-1}$)', labelpad=4)
cbar.ax.tick_params(length=2.5, width=0.6)
cbar.ax.set_yticklabels(
    ['0', '100', '250', '', '500', '', '750', '1000', '1250', '1500'])

# Tight layout
gs.tight_layout(figure, h_pad=0, w_pad=1)

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
