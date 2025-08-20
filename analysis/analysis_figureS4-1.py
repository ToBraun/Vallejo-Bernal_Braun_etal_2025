# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# All rights reserved.
# GNU General Public License v3.0.

'''
This script subsets the PIKART catalog to its relaxed AR shapes.

Run this script from the command line with:

    python analysis_figureS4-1.py
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
relaxed_file = os.path.join(
    config['output_folder'], 'continuous_trackids', 
    output_base_name %('AR_trackids', '%s', 'csv'))

# Output file
output_file = '%s-%s_PIKARTV1_relaxed_ARs.csv' %(year_start, year_end)

# %% MERGE CSV FILES

ardf_list = []
for year in range(year_start, year_end + 1):
    
    # Load lagrangian version of PIKART
    ardf = read_csv_record(relaxed_file %year)
    ardf['time'] = pd.to_datetime(ardf['time'])
    
    ardf = ardf[ardf['is_relaxed']]
    
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
