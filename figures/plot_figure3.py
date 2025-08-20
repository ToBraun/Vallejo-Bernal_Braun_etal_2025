# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. 3 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure3.py
'''

# %% IMPORT MODULES

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from cmcrameri import cm
from datetime import timedelta
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_functions import plot_landmasses, maxfactor_map_fig3, read_csv_record

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
central_lon = 0
main = 127367    # Main AR
merge = 127366    # AR that attaches
split = 127382   # AR that dettaches
color_map = cm.managua

# Files
pikart_file = root_folder + 'data/figure02/AR_land_intersection_ERA5_0p5deg_6hr_2019.csv'
figure_file = root_folder + 'manuscript/PIKART_Figure3.png'

# %% LOAD DATA

# Open csv file with AR records
ardf = read_csv_record(pikart_file)
ardf['time'] = pd.to_datetime(ardf['time'])

# Extract the AR tracks
main = ardf[ardf['trackid'] == main]
merge = ardf[ardf['trackid'] == merge]
split = ardf[ardf['trackid'] == split]

# %% COLORMAP SETTINGS

# Create list of dates
initial_date = main['time'].iloc[0]
final_date = main['time'].iloc[-1]
dates = pd.DataFrame(pd.date_range(initial_date, final_date, freq='6h'), columns=['time'])
levels = np.size(dates,0)

# Define color positions
color_pos = np.arange(0, levels * round(256*0.8/(levels-1)), round(256*0.8/(levels-1)))

# Create colorbar
colors = []
for c in range(len(color_pos)):
    colors.append(color_map(color_pos[c]))
colors = mpl.colors.ListedColormap(colors)

bounds = np.arange(0,len(color_pos),1)
ticks = np.arange(0,len(color_pos),4)
norm = mpl.colors.BoundaryNorm(bounds, colors.N)

# %% PLOT MULTIPANEL FIGURE

# Create the figure
figure = plt.figure(figsize=(180/25.4,65/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# %% PANEL C

# Create the subplot axis
# add_axes([xmin,ymin,dx,dy])
ax = figure.add_axes([0.275, 0.05, 0.45, 0.9], projection=proj)

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.5)

# Extract the dates of the track
track_dates = (dates.reset_index().merge(main, how='inner').set_index('index')).index

# Plot the temporal evolution of the track
for k in range(len(track_dates)):

    color_edge = colors(track_dates[k])

    # Plot axis
    ax.plot(main['axis_lon'].iloc[k], main['axis_lat'].iloc[k],
            color=color_edge, linewidth=0.65,
            alpha=1, transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map_fig3(ax, draw_grid=False)
ax.set_extent([-90, 20, 15, 69], ccrs.PlateCarree())
ax.set_aspect('equal')
#title =
ax.set_title('Spatiotemporal evolution of the AR trajectory', pad=15, fontweight='bold')
ax.text(0.5, 1.04,
        '%s to %s' %(main['time'].iloc[0].strftime('%d %b %Y, %H UTC'),
                     main['time'].iloc[-1].strftime('%d %b %Y, %H UTC')),
        ha='center', transform=ax.transAxes)

# Instance used to place the colorbar
divider = make_axes_locatable(ax)

# %% PANEL A

# Create the subplot axis
ax = figure.add_axes([0.04, 0.51, 0.19, 0.42], projection=proj)

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.5)

# Extract the merging times
time0 = merge['time'].iloc[-1]

# Extract the tracks during the merging
track_main = main[main['time'] == time0].iloc[0]
track_merge = merge[merge['time'] == time0].iloc[0]

# Prepare color
color_edge = color_map(color_pos[dates[dates['time']==time0].index[0]])
color_face = color_edge[:-1] + (0.2,)

# Plot contours
polygon = mpatches.Polygon(list(zip(track_main['contour_lon'], track_main['contour_lat'])),
                           closed=True, ec=color_edge, fill=True, lw=1,
                           fc=color_face, transform=ccrs.Geodetic())
ax.add_patch(polygon)

polygon = mpatches.Polygon(list(zip(track_merge['contour_lon'], track_merge['contour_lat'])),
                           closed=True, ec='gray', fill=True, lw=1,
                           fc='gainsboro', transform=ccrs.Geodetic())
ax.add_patch(polygon)

# Plot axes
ax.plot(track_main['axis_lon'], track_main['axis_lat'], color=color_edge,
        linewidth=0.65, transform=ccrs.Geodetic())

ax.plot(track_merge['axis_lon'], track_merge['axis_lat'], color='gray',
        linewidth=0.65, transform=ccrs.Geodetic())

# Tag the ARs
ax.text(-46, 41, 'AR1', color=color_edge, ha='right', weight='bold', transform=ccrs.Geodetic())
ax.text(-20, 41, 'AR2', color='gray', ha='right', weight='bold', transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map_fig3(ax, draw_grid=False, ll=True, bl=False)
ax.set_extent([-70, 0, 20, 60], ccrs.PlateCarree())
ax.set_aspect('equal')

ax.set_title('AR merger', pad=15, fontweight='bold')
ax.text(0.5, 1.08, time0.strftime('%d %b %Y, %H UTC'), ha='center',
        transform=ax.transAxes)

# %% PANEL B

# Create the subplot axis
ax = figure.add_axes([0.04, 0.07, 0.19, 0.42], projection=proj)

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.5)

# Extract the merging times
timef = time0 + timedelta(hours=6)

# Extract the tracks during the merging
track_main = main[main['time'] == timef].iloc[0]

# Prepare color
color_edge = color_map(color_pos[dates[dates['time']==timef].index[0]])
color_face = color_edge[:-1] + (0.2,)

# Plot contours
polygon = mpatches.Polygon(list(zip(track_main['contour_lon'], track_main['contour_lat'])),
                           closed=True, ec=color_edge, fill=True, lw=1,
                           fc=color_face, transform=ccrs.Geodetic())
ax.add_patch(polygon)

# Plot axes
ax.plot(track_main['axis_lon'], track_main['axis_lat'], color=color_edge,
        linewidth=0.65, transform=ccrs.Geodetic())

# Tag the ARs
ax.text(-20, 41, 'AR1', color=color_edge, ha='right', weight='bold', transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map_fig3(ax, draw_grid=False, ll=True, bl=True)
ax.set_extent([-70, 0, 20, 60], ccrs.PlateCarree())
ax.set_aspect('equal')

ax.text(0.5, 1.08, timef.strftime('%d %b %Y, %H UTC'), ha='center',
        transform=ax.transAxes)

# %% PANEL D

# Create the subplot axis
ax = figure.add_axes([0.77, 0.51, 0.218, 0.42], projection=proj)

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.5)

# Extract the splitting times
time0 = split['time'].iloc[0] - timedelta(hours=6)

# Extract the tracks during the splitting
track_main = main[main['time'] == time0].iloc[0]

# Prepare color
color_edge = color_map(color_pos[dates[dates['time']==time0].index[0]])
color_face = color_edge[:-1] + (0.2,)

# Plot contours
polygon = mpatches.Polygon(list(zip(track_main['contour_lon'], track_main['contour_lat'])),
                            closed=True, ec=color_edge, fill=True, lw=1,
                            fc=color_face, transform=ccrs.Geodetic())
ax.add_patch(polygon)

# Plot axes
ax.plot(track_main['axis_lon'], track_main['axis_lat'], color=color_edge,
        linewidth=0.65, transform=ccrs.Geodetic())

# Tag the ARs
ax.text(-25, 44, 'AR1', color=color_edge, ha='right', weight='bold', transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map_fig3(ax, draw_grid=False, ll=True, bl=False)
ax.set_extent([-55, 25, 25, 65], ccrs.PlateCarree())
ax.set_aspect('equal')

ax.set_title('AR separation', pad=15, fontweight='bold')
ax.text(0.5, 1.08, time0.strftime('%d %b %Y, %H UTC'), ha='center',
        transform=ax.transAxes)

# %% PANEL E

# Create the subplot axis
ax = figure.add_axes([0.77, 0.07, 0.218, 0.42], projection=proj)

# Plot the landmasses
ax = plot_landmasses(ax, lw=0.5)

# Extract the merging times
timef = time0 + timedelta(hours=6)

# Extract the tracks during the merging
track_main = main[main['time'] == timef].iloc[0]
track_split = split[split['time'] == timef].iloc[0]

# Prepare color
color_edge = color_map(color_pos[dates[dates['time']==timef].index[0]])
color_face = color_edge[:-1] + (0.2,)

# Plot contours
polygon = mpatches.Polygon(list(zip(track_main['contour_lon'], track_main['contour_lat'])),
                            closed=True, ec=color_edge, fill=True, lw=1,
                            fc=color_face, transform=ccrs.Geodetic())
ax.add_patch(polygon)

polygon = mpatches.Polygon(list(zip(track_split['contour_lon'], track_split['contour_lat'])),
                            closed=True, ec='gray', fill=True, lw=1,
                            fc='gainsboro', transform=ccrs.Geodetic())
ax.add_patch(polygon)

# Plot axes
ax.plot(track_main['axis_lon'], track_main['axis_lat'], color=color_edge,
        linewidth=0.65, transform=ccrs.Geodetic())

ax.plot(track_split['axis_lon'], track_split['axis_lat'], color='gray',
        linewidth=0.65, transform=ccrs.Geodetic())

# Tag the ARs
ax.text(-25, 45, 'AR1', color=color_edge, ha='right', weight='bold', transform=ccrs.Geodetic())
ax.text(-8, 57, 'AR3', color='gray', ha='right', weight='bold', transform=ccrs.Geodetic())

# Maxfactor the map
ax = maxfactor_map_fig3(ax, draw_grid=False, ll=True, bl=True)
ax.set_extent([-55, 25, 25, 65], ccrs.PlateCarree())
ax.set_aspect('equal')

ax.text(0.5, 1.08, timef.strftime('%d %b %Y, %H UTC'), ha='center',
        transform=ax.transAxes)

# %% LEGEND

# Colorbar
cax = divider.new_vertical(size = "5%", pad = 0.25, axes_class=plt.Axes, pack_start=True)
figure.add_axes(cax)

colorbar = figure.colorbar(mpl.cm.ScalarMappable(cmap=colors, norm=norm), cax=cax,
                           ticks=ticks, spacing='proportional', orientation='horizontal',
                           label='Days after the genesis of AR1')
colorbar.ax.set_xticklabels((ticks/4).astype('int'))

# Literals
plt.figtext(0.032, 0.895, '(a)', ha='right')
plt.figtext(0.032, 0.456, '(b)', ha='right')
plt.figtext(0.267, 0.895, '(c)', ha='right')
plt.figtext(0.762, 0.895, '(d)', ha='right')
plt.figtext(0.762, 0.456, '(e)', ha='right')

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
