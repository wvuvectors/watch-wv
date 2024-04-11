#!/bin/bash

#SBATCH -J analysis_test
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -p standby
#SBATCH -t 4:00:00
source /shared/software/conda/etc/profile.d/conda.sh

#example commandline usage : ./analysis_script_testing.sh RUNID path/to/rawdata 

cd /scratch/viv0001/COVID_ANALYSIS/SCRIPTS/

./analysis_script_testing.sh SARS_20240326 /scratch/viv0001/COVID/SARS_20240326

