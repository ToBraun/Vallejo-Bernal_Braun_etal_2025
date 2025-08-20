# Copyright (C) 2023 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script contains useful functions to plot.
'''

# %% IMPORT MODULES

import numpy as np
import pandas as pd
import re
import geopandas as gpd
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import ListedColormap
from cartopy.util import add_cyclic_point
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
from matplotlib.ticker import MultipleLocator

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

# %% FOLDERS

world_continents_folder = ...

# %% READ CSV RECORD

def read_csv_record(file_name):

    def conv_array(text):
        '''Convert array texts to ndarray'''
        text = text.replace('[','').replace(']','')
        array = np.array(text.split()).astype('float')
        return array

    conv_keys = ['contour_lon', 'contour_lat', 'axis_lon', 'axis_lat', 
                 'contour_x', 'contour_y', 'axis_x', 'axis_y']

    converters = dict([(key_i, conv_array) for key_i in conv_keys])

    dtypes={'id': 'int', 
            'time': 'str',
            'deleted': 'str',
            'is_relaxed': 'bool',
            'why_relaxed': 'str'}

    ardf = pd.read_csv(file_name, dtype=dtypes, converters=converters)

    return ardf

# %% PLOT LANDMASSES

# Add coast lines
def plot_landmasses(map_axes, lw=0.6, ec='darkgrey', hatch=True):
    
    # Change hatch linewidth
    mpl.rcParams.update({'hatch.linewidth': lw})

    # Load continental and insular shapefiles
    shape_cont = gpd.read_file(world_continents_folder + 'world_continental_3000sqkm.shp')
    shape_ins = gpd.read_file(world_continents_folder + 'world_insular_3000sqkm.shp')
    continents = ['North America', 'South America', 'Europe', 'Africa',
                  'Asia', 'Australia', 'Oceania', 'Antarctica']

    for i in range(len(continents)):
        
        # Plot the continents
        contdf = shape_cont[shape_cont['CONTINENT'] == continents[i]]
        map_axes.add_geometries(contdf['geometry'], crs=ccrs.Geodetic(), linewidth=lw,
                          linestyle='solid', edgecolor=ec, facecolor='none')
        
        if hatch:

            # Plot the islands
            insdf = shape_ins[shape_ins['CONTINENT'] == continents[i]]
            p = insdf['geometry'].iloc[0]
            map_axes.add_geometries(p.reverse(), crs=ccrs.Geodetic(), linewidth=lw,
                              linestyle='solid', edgecolor=ec, facecolor='none',
                              hatch='\\\\\\\\')
            
        else:
            # Plot the islands
            insdf = shape_ins[shape_ins['CONTINENT'] == continents[i]]
            p = insdf['geometry'].iloc[0]
            map_axes.add_geometries(p.reverse(), crs=ccrs.Geodetic(), linewidth=lw,
                              linestyle='solid', edgecolor=ec, facecolor='none')

    return map_axes

# %% GRIDS AND TICKS

# Add grid, ticks and tick labels
def maxfactor_map(map_axes, draw_grid=True, ll=True, bl=True):

    # Major grid
    if draw_grid:
        alpha_lines = 1
    else:
        alpha_lines = 0

    gl = map_axes.gridlines(
        draw_labels=True, crs=ccrs.PlateCarree(), x_inline=False, y_inline=False, 
        linewidth=0.5, color='darkgrey', alpha=alpha_lines, linestyle=':', zorder=0)
    gl.top_labels = False
    gl.bottom_labels = bl
    gl.right_labels = False
    gl.left_labels = ll
    gl.xlocator = mticker.FixedLocator(
        [-150,-120,-90,-60,-30,0,30,60,90,120,150,180])
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
    gl.xformatter = LongitudeFormatter(direction_label=False)
    gl.yformatter = LatitudeFormatter(direction_label=False)
    gl.rotate_labels=False
    gl.xpadding = 6
    gl.ypadding = 6

    # Major ticks
    map_axes.set_xticks([-150,-120,-90,-60,-30,0,30,60,90,120,150,180],
                        crs=ccrs.PlateCarree())
    map_axes.set_yticks([-90,-60,-30,0,30,60,90], crs=ccrs.PlateCarree())
    map_axes.set_xticklabels([])
    map_axes.set_yticklabels([])
    map_axes.tick_params('both', length=3.5, width=0.6)

    # Minor ticks
    map_axes.xaxis.set_minor_locator(MultipleLocator(10))
    map_axes.yaxis.set_minor_locator(MultipleLocator(10))
    map_axes.tick_params('both', length=2.5, which='minor')

    return map_axes

# Addapted to Figure 3
def maxfactor_map_fig3(map_axes, draw_grid=True, ll=True, bl=True):

    # Major grid
    if draw_grid:
        lw_major = 0.5
    else:
        lw_major = 0

    gl = map_axes.gridlines(
        draw_labels=True, crs=ccrs.PlateCarree(), x_inline=False, y_inline=False, 
        linewidth=lw_major, color='darkgrey', alpha=1, linestyle=':')
    gl.top_labels = False
    gl.bottom_labels = bl
    gl.right_labels = False
    gl.left_labels = ll
    gl.xlocator = mticker.FixedLocator(
        [-120,-90,-60,-30,0,30])
    gl.ylocator = mticker.FixedLocator([0,20,40,60,80])
    gl.xformatter = LongitudeFormatter(direction_label=False)
    gl.yformatter = LatitudeFormatter(direction_label=False)
    gl.rotate_labels=False
    gl.xpadding = 4
    gl.ypadding = 4

    # Major ticks
    map_axes.set_xticks([-120,-90,-60,-30,0,30],
                        crs=ccrs.PlateCarree())
    map_axes.set_yticks([0,20,40,60,80], crs=ccrs.PlateCarree())
    map_axes.set_xticklabels([])
    map_axes.set_yticklabels([])
    map_axes.tick_params('both', length=4)

    # Minor ticks
    map_axes.xaxis.set_minor_locator(MultipleLocator(10))
    map_axes.yaxis.set_minor_locator(MultipleLocator(10))
    map_axes.tick_params('both', length=2.5, which='minor')

    return map_axes

# Figure validation of AR detection
def maxfactor_map_fig6(map_axes, draw_grid=True, ll=True, bl=True):

    # Major grid
    if draw_grid:
        alpha_lines = 1
    else:
        alpha_lines = 0

    # Grid lines
    lw_major = 0.5

    # Major grid
    gl = map_axes.gridlines(
        draw_labels=True, crs=ccrs.PlateCarree(), x_inline=False, y_inline=False, 
        linewidth=lw_major, color='darkgrey', alpha=alpha_lines, linestyle=':', zorder=0)
    gl.top_labels = False
    gl.bottom_labels = bl
    gl.right_labels = False
    gl.left_labels = ll
    gl.xlocator = mticker.FixedLocator([120,140,160,180,-160,-140,-120,-100])
    gl.ylocator = mticker.FixedLocator([20,40,60])
    gl.xformatter = LongitudeFormatter(direction_label=False)
    gl.yformatter = LatitudeFormatter(direction_label=False)
    gl.rotate_labels=False
    gl.xpadding = 5
    gl.ypadding = 3

    # Major ticks
    map_axes.set_xticks([120,140,160,180,-160,-140,-120,-100], crs=ccrs.PlateCarree())
    map_axes.set_yticks([20,40,60], crs=ccrs.PlateCarree())
    map_axes.set_xticklabels([])
    map_axes.set_yticklabels([])
    map_axes.tick_params('both', length=3.5, width=0.6)

    # Minor ticks
    map_axes.xaxis.set_minor_locator(MultipleLocator(base=10, offset=5))
    map_axes.yaxis.set_minor_locator(MultipleLocator(10))
    map_axes.tick_params('both', length=2.5, which='minor')

    return map_axes

# Addapted to Figure 10
def maxfactor_fig10(map_axes, ll=True, bl=True):

    map_axes.set_xticks([-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    map_axes.set_yticks([-60,-30,0,30,60], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(direction_label=False)
    lat_formatter = LatitudeFormatter(direction_label=False)
    map_axes.xaxis.set_major_formatter(lon_formatter)
    map_axes.yaxis.set_major_formatter(lat_formatter)     
    map_axes.tick_params(axis='x', which='major', pad=6)

    # Major ticks
    if bl == False:
        map_axes.set_xticklabels([])
    if ll == False:
        map_axes.set_yticklabels([])
    
    map_axes.tick_params('both', length=2, width=0.6)

    return map_axes

#Figure IPART wrong contours
def maxfactor_map_figS5(map_axes, ll=True, bl=True):

    gl = map_axes.gridlines(
        draw_labels=True, crs=ccrs.PlateCarree(), x_inline=False,
        y_inline=False, linewidth=0.5, color='darkgrey',
        alpha=0, linestyle=':', zorder=0)
    gl.top_labels = False
    gl.bottom_labels = bl
    gl.right_labels = False
    gl.left_labels = ll
    gl.xlocator = mticker.FixedLocator([60,80,100,120,140])
    gl.ylocator = mticker.FixedLocator([40,60,80])
    gl.xformatter = LongitudeFormatter(direction_label=False)
    gl.yformatter = LatitudeFormatter(direction_label=False)
    gl.rotate_labels=False
    gl.xpadding = 6
    gl.ypadding = 4

    # Major ticks
    map_axes.set_xticks([60,80,100,120,140],
                        crs=ccrs.PlateCarree())
    map_axes.set_yticks([40,60,80], crs=ccrs.PlateCarree())
    map_axes.set_xticklabels([])
    map_axes.set_yticklabels([])
    map_axes.tick_params('both', length=3.5, width=0.6)

    # Minor ticks
    map_axes.xaxis.set_minor_locator(MultipleLocator(10))
    map_axes.yaxis.set_minor_locator(MultipleLocator(10))
    map_axes.tick_params('both', length=2.5, which='minor')

    return map_axes

def maxfactor_map_polar(map_axes):

    # Add grid
    gl = map_axes.gridlines(draw_labels=True, color='gray', linestyle='--',
                            linewidth=0.5, alpha=.5)
    
    plt.draw()
    
    # Reposition the tick labels
    # Labels at 150d meridian will be moved to 180d
    for ea in gl.label_artists:
        
        pos = ea.get_position()
        if pos[0] == 150:
            ea.set_position([90, pos[1]])
        
        text = ea.get_text()
        match = re.search(r'(\d+)', text)
        if match:
            value = int(match.group(1))
            if value % 10 != 0:
                ea.set_visible(False)
                
    return map_axes

# %% PLOT TIME SERIES WITH TREND

def plot_time_series(ax, years, time_series, trend, color, ylims, yticks, 
                     xtick_labels=False, add_title=False, title_str=''):
    
    # Plot time series
    ax.plot(years, time_series, color=color, linewidth=0.75)
    
    # Plot trend
    ax.plot(years, trend, color=color, alpha=.5, linewidth=1.5)
    
    # X axis
    ax.set_xlim(1940,2023)
    ax.set_xticks(np.arange(1940, 2040, 20))  
    
    if xtick_labels:
        ax.set_xlabel('Year')
    else:
        ax.set_xticklabels([])  # Hide x-axis labels
        
    # Y axis
    ax.set_ylim(ylims[0], ylims[-1])
    ax.set_yticks(yticks)
    
    # Axis ticks
    ax.tick_params(axis='both', length=2, width=0.6)
    
    # Grid
    ax.grid()
    
    # Title
    if add_title:
        ax.set_title(title_str, fontweight='bold')
        
# %% PLOT MAP WITH TRENDS

def plot_trends_map(figure, ax, lat, lon, map_data, signif_data, color_map, 
                    label, yticks, mask_nans=False):
        
    # Create the subplot axis
    ax.set_global()
    
    # Plot the landmasses
    ax = plot_landmasses(ax, lw=0.5, ec='gray')
    
    # Redefine non-uniform normalization to display real range of values, centered around zero
    norm = TwoSlopeNorm(vmin=np.floor(np.nanmin(map_data)), vcenter=0, vmax=np.floor(np.nanmax(map_data)))

    # Plot the shadding
    map_data, lon_plot = add_cyclic_point(map_data, coord=lon, axis=1)
    
    map_shad = ax.pcolormesh(lon_plot, lat, map_data, shading='nearest',
                             cmap=color_map, norm=norm, 
                             transform=ccrs.PlateCarree(), zorder=0)
    
    if mask_nans:
        
        # Overlay yellow on NaNs
        nan_mask = np.isnan(map_data)
        if np.any(nan_mask):
            ax.pcolormesh(lon_plot, lat, nan_mask, shading='nearest',
                          cmap=ListedColormap(['none', 'cornsilk']), vmin=0, 
                          vmax=1, transform=ccrs.PlateCarree(), zorder=1)
    
    # Plot the significance
    scatx, scaty = np.meshgrid(lon, lat)
    ax.scatter(scatx, scaty, marker='o', c='black', edgecolors='none', 
                s=0.2*signif_data, transform=ccrs.PlateCarree())
        
    # Maxfactor the map
    ax = maxfactor_fig10(ax, ll=True, bl=True)
    ax.set_aspect('equal')
        
    # Color bar
    pos = ax.get_position()
    ax_cb = figure.add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
    figure.add_axes(ax_cb)
    cbar = plt.colorbar(map_shad, cax=ax_cb, extend='max')
    cbar.ax.set_yticks(yticks)
    cbar.ax.tick_params(length=2, width=0.6)
    cbar.set_label(label=label, labelpad=4)
    