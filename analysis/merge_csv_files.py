# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script compiles a reduced version of the PIKART lagrangian catalog, 
containing the full period from 1940 to 2023, but only the variables necessary
for the analysis presented in Vallejo-Bernal & Braun et al., (2025).

Run this script from the command line with:

    python merge_csv_files.py
'''

# %% IMPORT MODULES

import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))
from ardt_functions import read_csv_record

# %% PARAMETERS

# Load configuration file
config = yaml.safe_load(open('../config_PIKART.yml'))

# Years
year_start = config['year_start']
year_end = config['year_end']

# Input files
output_base_name = config['output_base_name']
lagrangian_file = os.path.join(
    config['output_folder'], 'PIKARTV1_lagrangian', 
    output_base_name %('PIKARTV1_lagrangian', '%s', 'csv'))

# Output file
output_file = '%s-%s_PIKARTV1_Vallejo-Bernal_Braun_etal_2025.csv' %(year_start, year_end)

# %% COLUMNS TO DROP

# Continents
continents = np.array(['north_america', 'south_america', 'europe', 'africa', 
                       'asia', 'australia', 'oceania', 'antarctica'])
conti_str = ['na', 'sa', 'eu', 'af', 'as', 'au', 'oc', 'an']

# Land intersection columns
li_cols = []
for i in range(len(continents)):
    for comp in ['conti', 'insu']:
        for v in ['ivt', 'ivtu', 'ivtv']:
            li_cols.append('%s_%s_%s' %(comp, conti_str[i], v))

# Columns to drop
drop_cols = [##### Shape
             'contour_lon', 'contour_lat',
             'axis_lon', 'axis_lat',
             'head_lon', 'head_lat',
             'tail_lon', 'tail_lat',
             ##### Geometric properties
             'area', 'length', 'width', 'perimeter', 
             'isop_quotient', 'lw_ratio',
             ##### Physical properties
             'core_lon', 'core_lat',
             'core_ivt', 'core_ivtu', 'core_ivtv', 
             'cent_speed',
             ##### Intersection with landmasses
             'water', 'north_america', 'south_america', 'europe', 
             'africa', 'asia', 'australia', 'oceania', 'antarctica',
             'continent']

drop_cols = drop_cols + li_cols

# AR classification
drop_cols = drop_cols + ['is_relaxed', 'land_inters'] 

# %% MERGE CSV FILES

ardf_list = []
for year in range(year_start, year_end + 1):
    
    # Load lagrangian version of PIKART
    ardf = read_csv_record(lagrangian_file %year)
    ardf['time'] = pd.to_datetime(ardf['time'])
    
    ardf = ardf.drop(columns=drop_cols)
    
    # Append AR track
    ardf_list.append(ardf)
    
    print(year)
    
# Merge dataframes
ardf_full = pd.concat(ardf_list, axis=0, ignore_index=True)

# Necessary: to remove summarization in csv file
if sys.version_info.major == 2:
    np.set_printoptions(threshold = np.inf)
elif sys.version_info.major == 3:
    np.set_printoptions(threshold = sys.maxsize)

# Save catalog
ardf_full.to_csv(output_file, index=False)
