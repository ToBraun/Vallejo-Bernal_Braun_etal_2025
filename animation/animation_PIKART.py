# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script creates an animation of the atmospheric river (AR) trajectories 
from a selected year of the Lagrangian PIKART catalog.

**Inputs**
    - A csv file with AR trajectories from the Lagrangian PIKART catalog, named
      following the convention:
          AR_trackids_{dataset}_{spatialres}_{timeres}_{year}.csv

**Outputs**
    - A mp4 file visualizing the AR trajectories for the specified year.

**Usage**
    1. Configure initial settings in the file "config_PIKART.yml".
    2. In the submission script "submit_animation.sh", specify the year to 
       animate by editing the line:

           python animation_PIKART.py {year}

    3. Run the script from the command line with:

           sbatch submit_animation.sh
'''

# %% IMPORT MODULES

import os
import sys
import yaml
import numpy as np
import pandas as pd
import random
import xarray as xr
import datetime as dt
import time as tm
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from cmcrameri import cm
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))
from ardt_functions import read_csv_record

# Change default tick direction
params = {'xtick.direction': 'in',
          'ytick.direction': 'in',}
plt.rcParams.update(params)
mpl.rcParams['font.size'] = 11

# %% PARAMETERS

# Get year from the submission code
year = int(sys.argv[1])

print('PIKART MOVIE %s\n' %year)

# Starting time
start = tm.time()
print('Starting time: %s' %dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Parameters
lon_var = config['lon_var']
lat_var = config['lat_var']
var = config['var']
data_base_name = config['data_base_name']
output_base_name = config['output_base_name']

# Files and folders
ivt_file = os.path.join(config['ivt_folder'], data_base_name %(var, year))
tracks_file = os.path.join(
    config['output_folder'], 'continuous_trackids',
    output_base_name %('AR_trackids', year, 'csv'))
wc_file = os.path.join(config['wc_folder'], config['wc_file'][:-3] + 'shp')

# %% LOAD DATA

# Load IVT dataset
ivt_dataset = xr.open_dataset(ivt_file)

# Load coordinates
lat = ivt_dataset[lat_var].values
lon = ivt_dataset[lon_var].values
time = pd.to_datetime(ivt_dataset['time'].values)

# Load IVT
ivt = ivt_dataset[var].values

# Add cyclic point to IVT
ivt, lon = add_cyclic_point(ivt, coord=lon, axis=2)

# Open csv file with AR trajectories
ardf = read_csv_record(tracks_file)
ardf['time'] = pd.to_datetime(ardf['time'])

# %% SHUFFLE COLORS

'''Generate n unique shuffled numbers with no consecutive values'''

def shuffled_no_consecutives(n):
    nums = list(range(n))
    
    # Try shuffling until no adjacent values are consecutive
    for _ in range(10000):  # limit attempts to avoid infinite loop
        random.shuffle(nums)
        if all(abs(a - b) > 3 for a, b in zip(nums, nums[1:])):
            return nums
    
    raise ValueError("Couldn't find a valid ordering after many attempts.")

# %% CREATE ANIMATION

# Select portion of colormap for the IVT
colors = cm.grayC_r(np.linspace(0.05, 0.35, 236))
color_map = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
color_min, color_max = 0, 750

# AR colors
tracks = ardf.drop_duplicates('trackid')
colors = cm.batlow(np.linspace(0, 1, 50))
rows = shuffled_no_consecutives(50)
color_map_ar = colors[rows,:]

# Load the continents shapefile
world_continents = ShapelyFeature(Reader(wc_file).geometries(), 
                                  ccrs.PlateCarree(), 
                                  facecolor='none')

# Create the figure
figure = plt.figure(figsize=(320/25.4,160/25.4), dpi=300)
ax = figure.add_subplot(
    1,1,1, projection=ccrs.PlateCarree(central_longitude=-150))

def update_frame(plot_idx):

    # Define the variables
    plot_time = time[plot_idx]
    ivt_t = ivt[plot_idx,:,:]
    
    # Extract the ARs
    ardf_t = ardf[ardf['time'] == plot_time]
    
    # Clear the axis
    ax.clear()
    
    # Plot the IVT field
    ax.set_extent([0, 360, -90, 90], ccrs.PlateCarree())
    cs = ax.contourf(lon, lat, ivt_t,
                     levels = np.linspace(color_min, color_max, 15+1),
                     cmap=color_map,
                     extend='max', transform=ccrs.PlateCarree())
    
    # Plot our coastlines
    ax.add_feature(world_continents, linewidth=0.75, linestyle='solid', 
                   edgecolor='grey')
    
    # Plot ARs
    for k in range(np.size(ardf_t,0)):
    
        # Extract AR
        ar = ardf_t.iloc[k]
        trackid = int(ar['trackid'])
        relaxed = ar['is_relaxed']
        color_ar = color_map_ar[trackid%50,:]
        color_face = np.copy(color_ar)
        color_face[-1] = 0.4
    
        # Mask IVT
        polygon = mpatches.Polygon(list(zip(ar['contour_lon'], ar['contour_lat'])),
                                   closed=True, ec='w', fill=True, lw=0,
                                   fc='w', transform=ccrs.Geodetic())
        ax.add_patch(polygon)
    
        # Plot contour
        linestyle=':' if relaxed else '-'
        polygon = mpatches.Polygon(list(zip(ar['contour_lon'], ar['contour_lat'])),
                                   closed=True, ec=color_ar, fill=True, lw=1.25,
                                   fc=color_face, ls=linestyle, 
                                   transform=ccrs.Geodetic())
        ax.add_patch(polygon)
    
        # Plot axis
        ax.plot(ar['axis_lon'], ar['axis_lat'], color=color_ar, linewidth=1, 
                alpha=0.75, transform=ccrs.Geodetic())
    
        # Plot centroid
        ax.plot(ar['centroid_lon'], ar['centroid_lat'], 'o', color=color_ar, ms=3,
                fillstyle='full', markeredgewidth=0, transform=ccrs.PlateCarree())
        ax.plot(ar['centroid_lon'], ar['centroid_lat'], 'o', color=color_ar, ms=8,
                fillstyle='none', markeredgewidth=1, transform=ccrs.PlateCarree())
    
    # Major grid
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), x_inline=False,
                      y_inline=False, linewidth=0, color='darkgrey', alpha=1, 
                      linestyle=':')
    gl.top_labels = False
    gl.bottom_labels = True
    gl.right_labels = False
    gl.left_labels = True
    gl.xlocator = mticker.FixedLocator([-160,-120,-80,-40,0,40,80,120,160])
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
    gl.xformatter = LongitudeFormatter(direction_label=False)
    gl.yformatter = LatitudeFormatter(direction_label=False)
    gl.rotate_labels=False
    gl.xpadding = 10
    gl.ypadding = 10
    
    # Minor grid
    ax.gridlines(xlocs=[-180,-140,-100,-60,-20,20,60,100,140,180], 
                 ylocs=[-80,-70,-50,-40,-20,-10,10,20,40,50,70,80], 
                 draw_labels=False, crs=ccrs.PlateCarree(), linewidth=0, 
                 color='silver', alpha=1, linestyle=':')
    
    # Major ticks
    ax.set_xticks([-160,-120,-80,-40,0,40,80,120,160], crs=ccrs.PlateCarree())
    ax.set_yticks([-90,-60,-30,0,30,60,90], crs=ccrs.PlateCarree())
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params('both', length=5)
    
    # Minor ticks
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params('both', length=3.5, which='minor')
    
    # Maxfactor the map
    title='PIKART catalog \N{MINUS SIGN} %s UTC \N{MINUS SIGN} https://ar.pik-potsdam.de' %(str(plot_time)[:13])
    ax.set_title(title, fontsize=14)
    
    # Add colorbar
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size='1.4%', pad=0.2, axes_class=plt.Axes)
    figure.add_axes(ax_cb)
    cbar = plt.colorbar(cs, cax=ax_cb)
    cbar.set_label('IVT ($\mathregular{kg \ m^{-1} \ s^{-1}}$)', labelpad=8)
    
    plt.tight_layout()

# Create and save animation
anim = animation.FuncAnimation(figure, update_frame, frames=np.size(time,0), 
                               interval=100)
f = 'PIKART_%s_movie.mp4' %(year)
writervideo = animation.FFMpegWriter(fps=2)
anim.save(f, writer=writervideo)

# Ending time
end = tm.time()

# Print processing time
print('Ending time: %s' %dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print('Processing time: %s hours' %round((end-start)/3600, 2))
