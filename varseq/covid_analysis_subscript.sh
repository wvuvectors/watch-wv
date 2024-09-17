#!/bin/bash

#SBATCH -J analysis_test
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -p standby
#SBATCH -t 4:00:00
source /shared/software/conda/etc/profile.d/conda.sh

#example commandline usage : ./analysis_script_testing.sh RUNID path/to/rawdata freyja_barcode_file

cd /scratch/viv0001/COVID_ANALYSIS/SCRIPTS/

./covid_analysis.sh SARS_20240326 /scratch/viv0001/COVID/SARS_20240326 usher_barcodes09_05_2024-00-47.feather

