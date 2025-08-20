#!/bin/bash

### Copyright (C) 2025 by
### Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
### All rights reserved.
### GNU General Public License v3.0.
###
### This script submits the Python script "animation_PIKART.py" to Foote.
###
### Run this script from the command line with: sbatch submit_animation.sh
###
### SUBMISSION CODE BEGIN ###

#SBATCH --qos=priority 
#SBATCH --job-name=PIKART_job
#SBATCH --output=animation_PIKART_%j.out
#SBATCH --error=animation_PIKART_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

python animation_PIKART.py 2019