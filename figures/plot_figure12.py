# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved..
# GNU General Public License v3.0.

'''
This script produces Fig. 12 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure12.py
'''

# %% IMPORT MODULES

import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
from cmcrameri import cm
import matplotlib.lines as mlines
from plot_functions import plot_time_series, plot_trends_map

import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

# Update Matplotlib parameters
plt.rcParams.update({'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'font.family': 'Arial'})

mpl.rcParams.update({'axes.linewidth': 0.6,
                     'axes.edgecolor': '#262626',
                     'xtick.color':'#262626',
                     'ytick.color':'#262626',
                     'grid.color': 'silver',
                     'grid.linewidth': 0.5,
                     'font.size': 8,
                     'axes.titlesize': 8.5})

root_folder = ...

# %% INITIAL SETTINGS

# Parameters
central_lon = -150
FACT = 1e4

# Folders and files
input_folder = root_folder + 'data/figure12/'
figure_file = root_folder + 'manuscript/PIKART_Figure12.png'

# %% LOAD DATA

slopes = np.load(input_folder + 'slopes.npy') 
signif = np.load(input_folder + 'signif.npy')
aggr_lons = np.load(input_folder + 'aggr_lons.npy')
aggr_lats = np.load(input_folder + 'aggr_lats.npy')

nh = np.load(input_folder + 'nh.npy')
trop = np.load(input_folder + 'trop.npy')
sh = np.load(input_folder + 'sh.npy')
nh_trend = np.load(input_folder + 'nh_trend.npy')
trop_trend = np.load(input_folder + 'trop_trend.npy')
sh_trend = np.load(input_folder + 'sh_trend.npy')
years = np.load(input_folder + 'years.npy')

# %% FIGURE LAYOUT

# Create the figure
figure = plt.figure(figsize=(180/25.4,180/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# Define relative row heights
h_small = 0.04
h_big = 0.21
wpad = 0.1  # horizontal space between columns
vpad_small = 0.008
vpad_big = 0.05

# Define function to place a pair of axes:
def add_row(row_height, vpad, column, geo_ax=False):
    global y_start
    y_start -= row_height
    if (column == 0) & (not geo_ax):
        ax = figure.add_axes([0.04, y_start, 0.42 - wpad/2, row_height])
    elif (column == 0) & geo_ax:
        ax = figure.add_axes([0.04, y_start, 0.42 - wpad/2, row_height], 
                             projection=proj)
    if (column == 1) & (not geo_ax):
        ax = figure.add_axes([0.50 + wpad/2, y_start, 0.42 - wpad/2, row_height])
    elif (column == 1) & geo_ax:
        ax = figure.add_axes([0.50 + wpad/2, y_start, 0.42 - wpad/2, row_height],
                             projection=proj)
    
    y_start -= vpad  # vertical padding after row
    
    return ax

# %% AR EVENTS

n=0
y_start = 0.95

# Northern hemisphere
ax_a1 = add_row(h_small, vpad_small, 0, False)
plot_time_series(ax_a1, years, nh[n,]/FACT, nh_trend[n,]/FACT, color='navy', 
                 ylims=[10, 15], yticks=(11,14), xtick_labels=False, 
                 add_title=True, title_str=r'AR events/year (x$\mathbf{10^4}$)')

# Tropics
ax_a2 = add_row(h_small, vpad_small, 0, False)
plot_time_series(ax_a2, years, trop[n,]/FACT, trop_trend[n,]/FACT, color='hotpink', 
                 ylims=[3.5, 11], yticks=(5,10), xtick_labels=False,
                 add_title=False, title_str='')

# Southern hemisphere
ax_a3 = add_row(h_small, vpad_big, 0, False)
plot_time_series(ax_a3, years, sh[n,]/FACT, sh_trend[n,]/FACT, color='darkgreen', 
                 ylims=[15, 21], yticks=(16,20), xtick_labels=True, 
                 add_title=False, title_str='')

# Literal
ax_a1.text(-0.02, 1.3, '(a)', ha='right', transform=ax_a1.transAxes)

# Map
print(np.nanmin(slopes[n,]*10))
print(np.nanmax(slopes[n,]*10))
yticks = [-2,-1,0,1,2,3]
ax_b = add_row(h_big, vpad_big, 0, True)
plot_trends_map(figure, ax_b, aggr_lats, aggr_lons, map_data=slopes[n,]*10, 
                signif_data=signif[n,], color_map='BrBG',
                label=r'$\Delta$AR events/decade', yticks=yticks)

# Literal
ax_b.text(-0.02, 1, '(b)', ha='right', transform=ax_b.transAxes)

# %% AVERAGE DURATION OF AR EVENTS

n=1
y_start = 0.95

# Northern hemisphere
ax_c1 = add_row(h_small, vpad_small, 1)
plot_time_series(ax_c1, years, nh[n,], nh_trend[n,], color='navy', 
                 ylims=[14.5, 17.5], yticks=(15,17), xtick_labels=False, 
                 add_title=True, title_str='Mean AR event duration (h)')

# Tropics
ax_c2 = add_row(h_small, vpad_small, 1)
plot_time_series(ax_c2, years, trop[n,], trop_trend[n,], color='hotpink', 
                 ylims=[17.5, 27], yticks=(20,25), xtick_labels=False, 
                 add_title=False, title_str='')

# Southern hemisphere
ax_c3 = add_row(h_small, vpad_big, 1)
plot_time_series(ax_c3, years, sh[n,], sh_trend[n,], color='darkgreen', 
                 ylims=[14.4, 17.5], yticks=(15,17), xtick_labels=True, 
                 add_title=False, title_str='')

# Literal
ax_c1.text(-0.02, 1.3, '(c)', ha='right', transform=ax_c1.transAxes)

# Map
print(np.nanmin(slopes[n,]*10))
print(np.nanmax(slopes[n,]*10))
yticks = [-2,-1,0,1,2,3]
ax_d = add_row(h_big, vpad_big, 1, True)
plot_trends_map(figure, ax_d, aggr_lats, aggr_lons, map_data=slopes[n,]*10, 
                signif_data=signif[n,], color_map=cm.bam, 
                label=r'$\Delta \tau$/decade (h)', yticks=yticks)

# Literal
ax_d.text(-0.02, 1, '(d)', ha='right', transform=ax_d.transAxes)

# %% AR-RELATED IVT

n=2
y_start = 0.47

# Northern hemisphere
ax_e1 = add_row(h_small, vpad_small, 0)
plot_time_series(ax_e1, years, nh[n,], nh_trend[n,], color='navy', 
                 ylims=[290, 385], yticks=(300,350), xtick_labels=False, add_title=True, 
                 title_str='Mean AR-related IVT (kg m$\mathbf{^{-1}}$ s$\mathbf{^{-1}}$)')

# Tropics
ax_e2 = add_row(h_small, vpad_small, 0)
plot_time_series(ax_e2, years, trop[n,], trop_trend[n,], color='hotpink', 
                 ylims=[380, 495], yticks=(400,450), xtick_labels=False, add_title=False, 
                 title_str='')

# Southern hemisphere
ax_e3 = add_row(h_small, vpad_big, 0)
plot_time_series(ax_e3, years, sh[n,], sh_trend[n,], color='darkgreen', 
                 ylims=[310, 440], yticks=(350,400), xtick_labels=True, add_title=False, 
                 title_str='')

# Literal
ax_e1.text(-0.02, 1.3, '(e)', ha='right', transform=ax_e1.transAxes)

# Map
print(np.nanmin(slopes[n,]*10))
print(np.nanmax(slopes[n,]*10))
yticks = [-12,-6,0,14,28]
ax_f = add_row(h_big, vpad_small, 0, True)
plot_trends_map(figure, ax_f, aggr_lats, aggr_lons, map_data=slopes[n,]*10, 
                signif_data=signif[n,], color_map=cm.vik_r, 
                label=r'$\Delta$IVT/decade (kg m$^{-1}$ s$^{-1}$)', yticks=yticks)

# Literal
ax_f.text(-0.02, 1, '(f)', ha='right', transform=ax_f.transAxes)

# %% NON AR-RELATED IVT

n=3
y_start = 0.47

# Northern hemisphere
ax_g1 = add_row(h_small, vpad_small, 1)
plot_time_series(ax_g1, years, nh[n,], nh_trend[n,], color='navy', 
                 ylims=[87, 103], yticks=(90,100), xtick_labels=False, add_title=True, 
                 title_str='Mean non-AR-related IVT (kg m$\mathbf{^{-1}}$ s$\mathbf{^{-1}}$)')

# Tropics
ax_g2 = add_row(h_small, vpad_small, 1)
plot_time_series(ax_g2, years, trop[n,], trop_trend[n,], color='hotpink', 
                 ylims=[176, 205], yticks=(180,200), xtick_labels=False, add_title=False, 
                 title_str='')

# Southern hemisphere
ax_g3 = add_row(h_small, vpad_big, 1)
plot_time_series(ax_g3, years, sh[n,], sh_trend[n,], color='darkgreen', 
                 ylims=[87, 108], yticks=(90,100), xtick_labels=True, add_title=False, 
                 title_str='')

# Literal
ax_g1.text(-0.02, 1.3, '(g)', ha='right', transform=ax_g1.transAxes)

# Map
print(np.nanmin(slopes[n,]*10))
print(np.nanmax(slopes[n,]*10))

yticks = [-4,-2,0,3,6]
ax_h = add_row(h_big, vpad_small, 1, True)
plot_trends_map(figure, ax_h, aggr_lats, aggr_lons, map_data=slopes[n,]*10, 
                signif_data=signif[n,], color_map=cm.vik_r, 
                label=r'$\Delta$IVT/decade (kg m$^{-1}$ s$^{-1}$)', 
                yticks=yticks, mask_nans=True)

# Literal
ax_h.text(-0.02, 1, '(h)', ha='right', transform=ax_h.transAxes)

# %% LEGEND AND LITERALS

# Markers
ghost = mlines.Line2D([], [], color='none', linewidth=0)

nh_ts = mlines.Line2D([0,0.8], [0,0], color='navy', linewidth=0.75)
nh_trend = mlines.Line2D([0,0.8], [0,0], color='navy', alpha=.5, linewidth=1.5)

trop_ts = mlines.Line2D([0,0.8], [0,0], color='hotpink', linewidth=0.75)
trop_trend = mlines.Line2D([0,0.8], [0,0], color='hotpink', alpha=.5, linewidth=1.5)

sh_ts = mlines.Line2D([0,0.8], [0,0], color='darkgreen', linewidth=0.75)
sh_trend = mlines.Line2D([0,0.8], [0,0], color='darkgreen', alpha=.5, linewidth=1.5)

time_series = mlines.Line2D([0,0.8], [0,0], color='grey', linewidth=0.75)
trend = mlines.Line2D([0,0.8], [0,0], color='grey', alpha=.5, linewidth=1.5)

# Add legend
figure.legend([time_series, trend, nh_ts, trop_ts, sh_ts], 
              ['Time series', 'Linear trend:', 'Northern Hemisphere', 'Tropics', 'Southern Hemisphere'], 
              loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0), edgecolor='none', framealpha=1)

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
