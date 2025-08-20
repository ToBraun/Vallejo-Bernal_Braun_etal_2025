# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de>
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script produces Fig. 11 of Vallejo-Bernal & Braun et al., (2025)

Run this script from the command line with:

    python plot_figure11.py
'''

# %% IMPORT MODULES

import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cmcrameri import cm
from cartopy.util import add_cyclic_point
from matplotlib.patches import Patch
from plot_functions import plot_landmasses, maxfactor_map

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
central_lon = -150

# Folders and files
input_folder = root_folder + 'data/figure11/'
figure_file = root_folder + 'manuscript/PIKART_Figure11.png'

# %% LOAD DATA

# Coordinates
lon = np.load(input_folder + 'lon_pikart.npy')
lat = np.load(input_folder + 'lat_pikart.npy')

bins = np.load(input_folder + 'bins.npy')

lifetimes_disp = np.load(input_folder + 'lifetimes_disp.npy')
iptimes_disp = np.load(input_folder + 'iptimes_disp.npy')

med_ivt_lc = np.load(input_folder + 'med_ivt_lc.npy')
med_ivt_lf = np.load(input_folder + 'med_ivt_lf.npy')
med_ivt_ip = np.load(input_folder + 'med_ivt_ip.npy')

rate_all = np.load(input_folder + 'rate_all.npy') 
rate_lf = np.load(input_folder + 'rate_lf.npy')
rate_ip = np.load(input_folder + 'rate_ip.npy')

txg = np.load(input_folder + 'txg.npy')
tyg = np.load(input_folder + 'tyg.npy')
g_kernel = np.load(input_folder + 'g_kernel.npy')

txt = np.load(input_folder + 'txt.npy')
tyt = np.load(input_folder + 'tyt.npy')
t_kernel = np.load(input_folder + 't_kernel.npy')

ip_freq = np.load(input_folder + 'ip_freq.npy')

# Add cyclic point
ip_freq, lon = add_cyclic_point(ip_freq, coord=lon, axis=1)

# Calculate number of days per year
ip_freq = 365*ip_freq/100

# Nonlinear scale
ip_freq = np.sqrt(ip_freq)

# %% PLOT

# Create the figure
figure = plt.figure(figsize=(180/25.4,180/25.4), dpi=500)

# Define the CartoPy projection
proj = ccrs.PlateCarree(central_longitude=central_lon)

# LF SCATTER
ax = figure.add_axes([0.06, 0.89, 0.86, 0.05]) 
ax.scatter(np.sort(lifetimes_disp), 400 + np.random.normal(0, 2, lifetimes_disp.size),
           alpha=.03, color='steelblue', s=25, marker='o')
ax.plot(np.nanmean(lifetimes_disp), 400, marker='*', color='darkblue', markersize=10)
ax.set_ylim(395, 405)
ax.set_xlim(0, 1)
ax.set_yticks([])
ax.set_xticks([])

# IP SCATTER
ax = figure.add_axes([0.06, 0.82, 0.86, 0.05]) 
ax.scatter(np.sort(iptimes_disp), 400 + np.random.normal(0, 2, iptimes_disp.size),
              alpha=.03, color='orchid', s=25, marker='o')
ax.plot(np.nanmean(iptimes_disp), 400, marker='*', color='purple', markersize=10)
ax.set_ylim(395, 405)
ax.set_xlim(0, 1)
ax.set_yticks([]);
ax.set_xticks([])

# IVT CURVES
ax = figure.add_axes([0.06, 0.58, 0.86, 0.22])
ax.plot(bins, med_ivt_lc, color='black', linewidth=1.5, label='All ARs')
ax.plot(bins, med_ivt_lf, color='cornflowerblue', linewidth=1.5, label='Land-falling ARs')
ax.plot(bins, med_ivt_ip, color='plum', linewidth=1.5, label='Inland-penetrating ARs')
ax.set_ylabel('IVT/IVT$_{t_0}$')
ax.set_xlabel('Life cycle')
ax.set_xlim(0,1)
ax.set_ylim(0.73, 1.1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1])

# IVT RATES
secx = ax.twinx()
secx.plot(bins[1:], rate_all, color='dimgray', linewidth=1.5, linestyle='dashed')
secx.plot(bins[1:], rate_lf, color='steelblue', linewidth=1.5, linestyle='dashed')
secx.plot(bins[1:], rate_ip, color='indigo', linewidth=1.5, linestyle='dashed')
secx.scatter(bins[1], rate_all[0], color='dimgray', s=20)
secx.scatter(bins[1], rate_lf[0], color='steelblue', s=20)
secx.scatter(bins[1], rate_ip[0], color='indigo', s=20)
secx.hlines(0, 0, 1, color='#262626', linewidth=0.75)
secx.set_ylabel('Gain/loss rate')
secx.set_xlim(0,1)
secx.set_ylim(-0.1, 0.05)
secx.set_yticks([-0.1, -0.05, 0, 0.05])

# Global legend on top of the entire plot
figure.legend(loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=3, edgecolor='none')

##### MAP

# Create the subplot axis
ax = figure.add_axes([0.06, 0.06, 0.86, 0.43], projection=proj)
ax.set_global()

# Plot the landmasses
ax = plot_landmasses(ax)

# AR genesis
ax.contour(txg, tyg, g_kernel, transform=ccrs.PlateCarree(), linewidths=1, 
           colors='mediumaquamarine', levels=2)

# AR termination
ax.contour(txt, tyt, t_kernel, transform=ccrs.PlateCarree(), linewidths=1, 
           colors='gray', levels=2)

# Inland penetration
map_ip = ax.pcolormesh(lon, lat, np.where(ip_freq==0, np.nan, ip_freq), 
                        cmap=cm.lapaz_r, shading='nearest', 
                        vmin=0, vmax=np.ceil(np.nanmax(ip_freq)),
                        transform=ccrs.PlateCarree(), 
                        zorder=0)

# Maxfactor the map
ax = maxfactor_map(ax, draw_grid=False)
ax.set_aspect('equal')

# Color bar
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size='1.8%', pad=0.1, axes_class=plt.Axes)
figure.add_axes(ax_cb)
cbar = plt.colorbar(map_ip, cax=ax_cb)
cbar.ax.tick_params(length=2.5, width=0.6)
cbar.ax.set_yticklabels((np.arange(0, np.ceil(np.nanmax(ip_freq) + 1)).astype(int))**2)
cbar.set_label(label='Inland-penetrating AR conditions (days/year)', labelpad=8)

# Legend elements
genesis = Patch(edgecolor='mediumaquamarine', facecolor='none', linewidth=1)
termination = Patch(edgecolor='gray', facecolor='none', linewidth=1)

# Add legend
plt.legend([genesis, termination], ['AR genesis', 'AR termination'],
            loc='upper left', bbox_to_anchor=(-22.5, 1.1), ncol=2, edgecolor='none', framealpha=1)

# Literals
plt.figtext(0.05, 0.964, '(a)', ha='right')
plt.figtext(0.05, 0.50, '(b)', ha='right')

# Save figure
figure.tight_layout()
plt.savefig(figure_file)
