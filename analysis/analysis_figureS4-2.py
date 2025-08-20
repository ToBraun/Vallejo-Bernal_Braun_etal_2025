# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script allows the selection of the relaxed AR shapes to be plotted in 
Figure S4 of Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python analysis_figureS4-2.py
'''

# %% IMPORT MODULES

import sys
import pandas as pd
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))
from ardt_functions import read_csv_record

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
year = 1979
central_lon = -150

# File
relaxed_file = root_folder + 'data/figureS4/1940-2023_PIKARTV1_relaxed_ARs.csv'

# %% PLOT COLLECTION OF ARs

def plot_ar_collection(ardf):

    for j in range(len(ardf)):

        # Create the figure
        figure = plt.figure(figsize=(180/25.4,110/25.4), dpi=300)

        # Define the CartoPy projection
        proj = ccrs.PlateCarree(central_longitude=central_lon)

        # Create the subplot axis
        ax = figure.add_subplot(1,1,1, projection=proj)
        ax.set_global()

        # Plot the landmasses
        ax.coastlines()

        # Extract AR
        ar = ardf.iloc[j]
        color_ar = 'steelblue'
        is_relax = (ar['why_relaxed'] != 'not_relaxed')

        # Plot contour
        linestyle=':' if is_relax else '-'
        polygon = mpatches.Polygon(list(zip(ar['contour_lon'], ar['contour_lat'])),
                                   closed=True, ec=color_ar, fill=True, lw=1.75,
                                   fc='none', ls=linestyle, transform=ccrs.Geodetic())
        ax.add_patch(polygon)

        # Plot axis
        ax.plot(ar['axis_lon'], ar['axis_lat'], color=color_ar, linewidth=1.25, alpha=0.75,
                transform=ccrs.Geodetic())

        # Add the track ID and the time
        trackid = ar['trackid']
        tail_lat = ar['axis_lat'][-1]
        tail_lon = ar['axis_lon'][-1]
        ax.annotate(str(int(trackid)) + ', ' + str(ar['time'])[:-6], xy=(tail_lon, tail_lat), xytext=(0, 0),
                          color='white', textcoords='offset points',
                          ha='center', va='center', fontsize=6, transform=ccrs.Geodetic(),
                          bbox = {'facecolor': 'dimgrey', 'edgecolor': 'dimgrey', 'pad': 1})

# %% LOAD DATA

# Open csv file with AR records
relaxed_ars = read_csv_record(relaxed_file)
relaxed_ars['time'] = pd.to_datetime(relaxed_ars['time'])

# %% EXTRACT RELAXED ARs

# Extract the year
relaxed_ars = relaxed_ars[relaxed_ars['time'].dt.year == year]

# Extract each category
short_ars = relaxed_ars[relaxed_ars['why_relaxed'] == 'short']
long_ars = relaxed_ars[relaxed_ars['why_relaxed'] == 'long']
big_ars = relaxed_ars[relaxed_ars['why_relaxed'] == 'big']
round_ars = relaxed_ars.loc[relaxed_ars['why_relaxed'] == 'round']
incoh_ars = relaxed_ars.loc[relaxed_ars['why_relaxed'] == 'incoherent_ivt_direction']

# %% PLOT ARs

# Sort big ARs
big_ars = big_ars.sort_values(by=['area'], ascending=False)

# Plot big ARs
plot_ar_collection(big_ars)

# Sort long ARs
long_ars = long_ars.sort_values(by=['length'], ascending=False)

# Plot long ARs
plot_ar_collection(long_ars.iloc[:10])

# Sort short ARs
short_ars = short_ars.sort_values(by=['length'], ascending=False)

# Plot short ARs
plot_ar_collection(short_ars.iloc[:20])

# Sort round ARs
round_ars = round_ars.sort_values(by=['isop_quotient'], ascending=False)

# Plot round ARs
plot_ar_collection(round_ars.iloc[:20])

# Sort incoherent ARs
incoh_ars = incoh_ars.sort_values(by=['length'], ascending=False)

# Plot round ARs
plot_ar_collection(incoh_ars.iloc[:20])

