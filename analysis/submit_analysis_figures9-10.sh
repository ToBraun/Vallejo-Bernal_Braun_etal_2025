#!/bin/bash

### Copyright (C) 2025 by
### Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
### All rights reserved.
### GNU General Public License v3.0.
###
### This script submits the Python script "analysis_figures9-10.py" to Foote
### It is parallelized to perform 36 tasks.
###
### Run this script from the command line with: 
### sbatch submit_analysis_figures9-10.sh
###
### SUBMISSION CODE BEGIN ###

#SBATCH --qos=standby
#SBATCH --time=2:00:00 
#SBATCH --job-name=PIKART_job
#SBATCH --array=0-35
#SBATCH --output=control_files/analysis_figures9-10_%j.out
#SBATCH --error=control_files/analysis_figures9-10_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

python analysis_figures9-10.py